### Performs alignment of channels to areas that the probe goes through for a given mouse id

#from this import d
from tkinter import ANCHOR, N
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pathlib
import os
import SimpleITK as sitk
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, QProgressBar,\
    QVBoxLayout, QWidget, QPushButton, QGridLayout, QCheckBox, QFormLayout, QFileDialog, QComboBox, QMessageBox
from PyQt5 import QtCore
from pyqtgraph.Qt import QtGui
from PyQt5.QtGui import QColor
from PyQt5 import QtWidgets
from qtrangeslider import QRangeSlider
import sys
from get_tissuecyte_info import get_tc_info
from warp_image import warp_channels
import xmltodict
from generate_metrics_paths import generate_metrics_path_days, generate_templeton_metric_path_days, generate_metrics_path_ephys
import visvis as vis
from warp_image import warp_execute, cluster_annotations
import pickle
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import pathlib
import threading
from qc_check import qcChecker
import colorbar as cb

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)
parser.add_argument('-t', '--templeton', help='Templeton experiment', required=False, default = "")
parser.add_argument('-r', '--record', help='Record Node number', required=False, default = "")
parser.add_argument('-o', '--oldDirStruct', help='Old Directory Structure', required=False, default="")
parser.add_argument('-s', '--synology', help='On Synology drive or not', required=False, default="")
parser.add_argument('--dr', help='Dynamic Routing experiment', required=False, default='')

SCALING_FACTOR = 1.5
DEFAULT_COLOR_VALUES = [[0, 256], [0, 256], [0, 1000]]
DEFAULT_COLOR_VALUES_CCF = [[0, 3000], [0, 3000], [0, 1000]]

# graph represents a metric plot that can be clicked on
# pyqtgraph documentation
class Graph(pg.GraphItem):
    def __init__(self, label=None):
        self.dragPoint = None
        self.dragOffset = None
        self.textItems = []
        self.pointPos = {}
        self.clicked = False
        pg.GraphItem.__init__(self)
        self.scatter.sigClicked.connect(self.onclick)
        self.label = label
        
    def setData(self, **kwds):
        self.text = kwds.pop('text', [])
        self.data = kwds

        if 'pos' in self.data:
            npts = self.data['pos'].shape[0]
            self.data['data'] = np.empty(npts, dtype=[('index', int)])
            self.data['data']['index'] = np.arange(npts)
            self.setTexts(self.text)
            self.updateGraph()
        
    def setTexts(self, text):
        for i in self.textItems:
            i.scene().removeItem(i)
        self.textItems = []
        for t in text:
            item = pg.TextItem(t)
            self.textItems.append(item)
            item.setParentItem(self)
        
    def updateGraph(self):
        pg.GraphItem.setData(self, **self.data)
        for i,item in enumerate(self.textItems):
            item.setPos(*self.data['pos'][i])
       
    """
    def mouseDragEvent(self, ev):
        if ev.button() != QtCore.Qt.LeftButton:
            ev.ignore()
            return
        
        if ev.isStart():
            # We are already one step into the drag.
            # Find the point(s) at the mouse cursor when the button was first 
            # pressed:
            pos = ev.buttonDownPos()
            print('Pos', pos)
            pts = self.scatter.pointsAt(pos)
            if len(pts) == 0:
                ev.ignore()
                return
            self.dragPoint = pts[0]
            ind = pts[0].data()[0]
            self.dragOffset = self.data['pos'][ind] - pos
        elif ev.isFinish():
            self.dragPoint = None
            return
        else:
            if self.dragPoint is None:
                ev.ignore()
                return
        
        ind = self.dragPoint.data()[0]
        print(ind)
        self.data['pos'][ind] = ev.pos() + self.dragOffset
        self.updateGraph()
        ev.accept()
        self.pointPos[ind] = ev.pos()
    
    def clicked(self, pts):
        print("clicked: %s" % pts)
    """
    def onclick(self, plot, points):
        if hasattr(self, 'lastPoint'):
            self.lastPoint.setBrush(QtGui.QBrush(QColor('grey')))
            self.lastPoint.setPen(QtGui.QPen(QColor('grey')))

        self.clicked = True
        points[0].setBrush(QtGui.QBrush(QColor('green')))
        points[0].setPen(QtGui.QPen(QColor('green')))
        self.lastPoint = points[0]
        self.scatterPoint = points[0].pos()
        
        channel = self.scatter.pointsAt(self.scatterPoint)[0].data()[0]
        if self.label is not None:
            self.label.setText('Channel %d' %(channel))

# class for metric plot
# each metric plot is represented by a plot display item which uses the graph class above
class PlotDisplayItem():
    # create default image
    def __init__(self, measurement, waveform_metrics, probe_annotations, mouse_id, metrics_list, label=None):
        self.width = int(4000) #default = coronal view
        self.height = int(4000)
        self.remove = True
        self.show = True

        self.probeAnnotations = probe_annotations
        self.mouseID = mouse_id
        self.metricsList = metrics_list
        self.otherPlots = []
        self.oldChannels = [] # stack for undoing lines
        self.texts = ['Channel %d' %(i) for i in range(383, -1, -1)]
        self.waveform_metrics = waveform_metrics
        self.measurement = measurement
        #self.generateMetrics(measurement)

        if label is not None:
            self.channelsPlot = Graph(label)
        else:
            self.channelsPlot = Graph()

        self.ogChannelsPlot = Graph()
        #self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
        self.textItem = pg.TextItem(measurement.upper(), anchor=(1, 1))

        #self.channels = self.channelsOriginal
    
    # Updates the metrics with the given path to the metrics csv
    # path: string, the path to the csv
    def updateMetrics(self, path, templeton=False):
        self.waveform_metrics = pd.read_csv(path)
        print('Metrics Path', path)
        self.generateMetrics(self.measurement, templeton=templeton)

    # generates the metrics based on the measurement passed in
    # measurement: string, the metric
    def generateMetrics(self, measurement, templeton=False, channels=None):
        if self.measurement == 'unit_density':
            self.processMetrics(templeton=templeton)
            self.frequencyCounts = self.waveform_metrics['peak_channel'].value_counts().sort_index().reset_index().to_numpy()
            #x = np.linspace(-10, 100, num=384)
            self.channelsRecorded = self.frequencyCounts[:, 0].tolist()
            self.channelsOriginal = []

            for i in range(384):
                if i in self.channelsRecorded:
                    ind = self.channelsRecorded.index(i)
                    self.channelsOriginal.append([i, self.frequencyCounts[ind, 1]])
                else:
                    self.channelsOriginal.append([i, 0])
           
            self.channelsOriginal = [[self.channelsOriginal[i][1], 384 - i - 1 + 256] for i in range(len(self.channelsOriginal))]
            
            x_val = [p[0] * 10 for p in self.channelsOriginal]
            conv = np.ones(10)

            smoothed = np.convolve(x_val, conv, mode='same')
            smoothed = smoothed / np.sum(conv)
            self.channelsOriginal = [[smoothed[i] - 100, self.channelsOriginal[i][1]] for i in range(384)]

            self.ogChannelsShift = [[p[0] - 200, p[1]] for p in self.channelsOriginal]

            self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
            self.ogChannelsPlot.setData(pos=np.array(self.ogChannelsShift, dtype=float), adj=np.array(self.adj, dtype=int))
        else:
            self.processMetrics(templeton=templeton)
            self.generateMetricChannels(measurement, scale_value=1/100, shift_value=250)

    def updatePlotChannelPositions(self, channel_positions:list):
        self.ogChannelsShift = [[p[0] - 200, p[1]] for p in channel_positions]
        self.ogChannelsShift = [[self.ogChannelsShift[i][0], 384 - i - 1 + 256] for i in range(384)]

        self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
        self.ogChannelsPlot.setData(pos=np.array(self.ogChannelsShift, dtype=float), adj=np.array(self.adj, dtype=int))

    # helper function to generate the metric channels
    # metric: string
    # scale_value: float, value to scale the plot by this amount when displayed
    # shit_value: int, value to shift the plot by this amount when displayed
    def generateMetricChannels(self, metric, scale_value, shift_value):
        #print('Metric', metric)
        peak_values = self.averageMetricsChannels['peak_channel'].values.tolist()
        values = self.averageMetricsChannels[metric]

        if 'velocity' in metric:
            values = (2 * (values - values.min()) / (values.max() - values.min())) - 1
        else:
            values = (values - values.min()) / (values.max() - values.min())

        conv = np.ones(10)
        self.channelsOriginal = []

        for i in range(384):
            if i in peak_values:
                """
                index = peak_values.index(i)
                if not pd.isna(values[index]):
                    self.channelsOriginal.append([values[index], 384 - i - 1 + 256])
                else:
                    self.channelsOriginal.append([np.nan, 384 - i - 1 + 256])
                #conv[index] = 1
                """
                index = peak_values.index(i)
                self.channelsOriginal.append([values[index], 384 - i - 1 + 256])
            else:
                self.channelsOriginal.append([0, 384 - i - 1 + 256])
        
        x_val = [p[0] for p in self.channelsOriginal]
        print('X_val', len(x_val))
        smoothed = np.convolve(x_val, conv, mode='same')
        smoothed = smoothed / np.sum(conv)
        #print(smoothed.shape)
        for i in range(384):
            if scale_value != 0:
                self.channelsOriginal[i] = [(smoothed[i] / scale_value) - shift_value, self.channelsOriginal[i][1]]
            else:
                self.channelsOriginal[i] = [(smoothed[i]) - shift_value, self.channelsOriginal[i][1]]
        
    def processMetrics(self, templeton=False):
        if not templeton:
            if 'quality' in self.waveform_metrics.columns:
                self.waveform_metrics = self.waveform_metrics.drop(columns=['epoch_name_quality_metrics', 'epoch_name_waveform_metrics', 'quality'])
            elif 'epoch_name_quality_metrics' in self.waveform_metrics.columns:
                self.waveform_metrics = self.waveform_metrics.drop(columns=['epoch_name_quality_metrics', 'epoch_name_waveform_metrics'])
        #self.wnorm = (2 * (self.waveform_metrics - self.waveform_metrics.min()) / (self.waveform_metrics.max() - self.waveform_metrics.min())) - 1

        #self.wnorm['peak_channel'] = self.waveform_metrics['peak_channel']
        self.averageMetricsChannels = (self.waveform_metrics.groupby('peak_channel').mean().reset_index())

        #print(self.averageMetricsChannels)
        #self.averageMetricsChannels = (self.averageMetricsChannels - self.averageMetricsChannels.mean()) / self.averageMetricsChannels.std()
        #self.averageMetricsChannels = self.averageMetricsChannels.rolling(rolling_value, win_type='boxcar').mean()

    # when a new probe is displayed
    def resetPlot(self, remove_probe=False):
        self.channels = self.channelsOriginal
        self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
        self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

    # helper function for shfiting points when new alignment is added
    # returns closest value to K
    # lst: list to look for closest
    # K: closest number to look for
    def closest(self, lst, K):
        if len(lst) > 0:
            return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

    # helper function to update channel positions for linear interpolation between 2 alignments
    # lp: list of points 
    # points_between: points that need to be updated with points from lp, same length as lp
    def replaceValues(self, lp, points_between):
        #print(points_between)
        for i in range(len(lp)):
            ind = self.channels.index(points_between[len(lp) - i - 1])
            #print(self.channels[ind])
            self.channels[ind][1] = lp[i]
            #print(self.channels[ind])

    # linearly space points between 2 alignments
    # pointsAdded: the list of alignments that have been made so far
    def linearSpacePoints(self, pointsAdded):
        if len(pointsAdded) > 1:
            sorted_linepts = sorted(pointsAdded)

            for i in range(len(sorted_linepts) - 1): # each time new line is added
                anchor_top = sorted_linepts[i]
                anchor_bottom = sorted_linepts[i + 1]

                #print(anchor_top, anchor_bottom)

                points_between = [p for p in self.channels if p[1] > anchor_top and p[1] < anchor_bottom] # current points between
                #print(len(points_between))
                lp = np.linspace(anchor_top + 1, anchor_bottom - 1, num=len(points_between)).tolist() # new linear interpolation
                #print(len(lp))
                #print(lp)
                #lp = [[p[0], int(p[1])] for p in lp]
                self.replaceValues(lp, points_between) # replace values with new interpolation
                #print(len(self.channels))

            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

    # for alignments above current channel
    # shift up, i.e. subtract
    # srt: list with sorted alignments already done
    # channel: y coord of point clicked on
    # diff: difference between channel and alignment 
    # flag: boolean, returns False if interpolation is done
    def linspaceAbove(self, srt, channel, diff, flag):
        for point in srt:
            if self.channelsPlot.scatterPoint[1] > point:
                points_between = [p for p in self.channels if p[1] < self.channelsPlot.scatterPoint[1] and p[1] > point]
                self.channels[channel][1] -= diff
                lp = np.linspace(point + 1, self.channels[channel][1] - 1, len(points_between))
                self.replaceValues(lp, points_between)
                flag = False
                break

        return flag

    # for alignments below current channel
    # shift down, i.e. add
    # srt: list with sorted alignments already done
    # channel: y coord of point clicked on
    # diff: difference between channel and alignment 
    # flag: boolean, returns False if interpolation is done
    def linspaceBelow(self, srt, channel, diff, flag):
        for point in srt:
            if self.channelsPlot.scatterPoint[1] < point:
                points_between = [p for p in self.channels if p[1] > self.channelsPlot.scatterPoint[1] and p[1] < point]
                self.channels[channel][1] += diff
                lp = np.linspace(self.channels[channel][1] + 1, point - 1, len(points_between))
                self.replaceValues(lp, points_between)
                flag = False
                break

        return flag

    # scales the plots using the given scale factors below
    def updateDisplay(self, probe, linepts, intensity_values, keep_y=False, old_y=None, points_added=None):
        if not keep_y:
            ap_scale = 1
            dv_scale = 1/0.94
            lr_scale = 1
 
            # get vectors in xyz
            self.dz = (linepts[-1, 0] - linepts[0, 0])
            self.dy = (linepts[-1, 1] - linepts[0, 1])
            self.dx = (linepts[-1, 2] - linepts[0, 2])
            self.vector = np.array([self.dz, self.dy, self.dx])
            
            # apply scale
            self.dzScale = self.dz * ap_scale
            self.dyScale = self.dy * dv_scale
            self.dxScale = self.dx * lr_scale
            self.vectorScale = np.array([self.dzScale, self.dyScale, self.dxScale])

            # get scale factor
            self.scale = np.linalg.norm(self.vectorScale) / np.linalg.norm(self.vector)
            #print(self.scale)
            #channels_scale = np.linspace(0, self.scale, 384)

            newPoints = []
            for i in range(len(self.channelsOriginal)):
                newPoints.append([self.channelsOriginal[i][0], self.channelsOriginal[i][1] * self.scale]) # apply scale factor
        
            #print('Initial Probe length in pixels', newPoints[0][1] - newPoints[-1][1])
            self.ogPoints = newPoints
            self.channels = newPoints
            self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
        else: # some alignment has been done already, so use existing coordinates
            newPoints = [] 
            for i in range(len(self.channelsOriginal)):
                newPoints.append([self.channelsOriginal[i][0], old_y[i][1]])

            self.ogPoints = newPoints
            self.channels = newPoints
            self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
            self.linearSpacePoints(points_added)

# class to do the alignment between channels and regions
class VolumeAlignment(QWidget):
    def __init__(self, mouse_id, dr=False, templeton=False, old_struct=False, synology=False, base_path=None, record_node=None):
        super().__init__()
        self.mouseID = mouse_id
        self.title = 'Volume Alignment for Mouse {}'.format(self.mouseID)
        self.probe = 'ProbeA'
        self.templeton = templeton
        #self.waveform = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.{}-AP/waveform_metrics.csv'.format(self.probe)))
        #self.metrics = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.{}-AP/metrics_test.csv'.format(self.probe)))
        
        self.basePath = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
        self.templeBasePath = pathlib.Path('//allen/programs/mindscope/workgroups/templeton/TTOC/pilot recordings')
        self.noRecord = False

        if not templeton:
            if not dr:
                self.waveMetricsPath = generate_metrics_path_ephys(self.basePath, self.mouseID)
            else:
                self.waveMetricsPath = generate_metrics_path_days(self.basePath, self.mouseID)

            #self.days = sorted(list(self.waveMetricsPath.keys()))
            self.waveform_metrics = pd.read_csv(os.path.join(self.basePath, '1178173272_608671_20220518/1178173272_608671_20220518_probeB_sorted/continuous/Neuropix-PXI-100.0', 
                                                             'metrics.csv'))
        else:
            if record_node is None:
                self.waveMetricsPath = generate_metrics_path_days(base_path, self.mouseID)
                self.days = sorted(list(self.waveMetricsPath.keys()))
                self.waveform_metrics = pd.read_csv(os.path.join(self.basePath, '1178173272_608671_20220518/1178173272_608671_20220518_probeB_sorted/continuous/Neuropix-PXI-100.0', 
                                                                 'metrics.csv'))
                self.noRecord = True
            else:
                self.waveMetricsPath = generate_templeton_metric_path_days(self.mouseID, base_path, record_node, old_struct=old_struct, synology=synology)
                self.days = sorted(list(self.waveMetricsPath.keys()))
                self.waveform_metrics = pd.read_csv(os.path.join(self.templeBasePath,
                                                                '2022-07-26_14-09-36_620263/Record Node 101/experiment1/recording1/continuous/Neuropix-PXI-100.ProbeB-AP', 
                                                                 'waveform_metrics.csv'))

        self.initUI()
        self.displayRegion()
    
    def initUI(self):
        self.setWindowTitle(self.title)
        self.width = int(4000) #default = coronal view
        self.height = int(4000)
        self.remove = True
        self.showProbe = True
        self.showMask = True
        self.prevProbe = ''

        self.im8 = np.ones((self.height,self.width),dtype='uint8')*255
        self.image = pg.image(self.im8)
        self.image.ui.histogram.hide()
        self.image.ui.roiBtn.hide()
        self.image.ui.menuBtn.hide()
        self.image.setObjectName('image')

        self.imageMask = pg.image(self.im8)
        self.imageMask.ui.histogram.hide()
        self.imageMask.ui.roiBtn.hide()
        self.imageMask.ui.menuBtn.hide()
        self.imageMask.setObjectName('image')
        #self.image.setImage(im8.transpose())


        self.pointsAdded = [] # y coords of alignments done
        self.channelCoords = {}
        self.lineItems = [] # pyqtgraph item for alignments
        self.oldChannels = [] # stack for undoing lines
        self.distPoints = []
        self.anchorPts = []
        self.s = []
        self.textItems = []
        self.ccfTextItems = []
        self.ccfPlotItems = []
        self.probeItems = []
        self.anchorItems = []
        self.anchorPos = []
        self.plotImageItems = []
        self.plotImageTextItems = []
        self.plotImageCCFItems = []
        self.ccfAreaItems = []
        self.anchorText = []
        self.pos = []
        self.myFont = ImageFont.load_default()
        self.alignments = {}

        # important directories
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
        self.modelDirectory = pathlib.Path('//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56/atlasVolume')
        self.storageDirectory = pathlib.Path(os.path.join(self.workingDirectory, self.mouseID))

        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))
        print(self.field_file)
        self.reference_file = os.path.join( self.modelDirectory, 'average_template_25.nrrd')
        
        print('Loading app...')
        if not os.path.exists(os.path.join(self.storageDirectory, 'anchors')):
            os.mkdir(os.path.join(self.storageDirectory, 'anchors'))
        
        with open(os.path.join(self.workingDirectory, 'field_reference', 'name_map.pkl'), 'rb') as f:
            self.name_map = pickle.load(f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'acrnm_map.pkl'), 'rb') as f:
            self.acrnm_map = pickle.load(f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'color_map.pkl'), 'rb') as f:
            self.colormap = pickle.load(f)

        self.anno = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.workingDirectory, 'field_reference', 'ccf_ano.mhd')))
        
        # read reference image, deformation field, probe annotations
        self.reference = sitk.ReadImage( self.reference_file )
        self.field = sitk.ReadImage( self.field_file)

        if os.path.exists(os.path.join(self.storageDirectory, 'reassigned', 'probe_annotations_{}_reassigned.csv'.format(self.mouseID))):
            self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'reassigned', 'probe_annotations_{}_reassigned.csv'.format(self.mouseID)))
        else:
            self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))

        self.metricsList = ['Spread', 'Amplitude', 'Isi_Viol', 'NN_Hit_Rate', 'Firing_Rate', 'Presence_Ratio', 'Velocity_Above', 'Velocity_Below']
        self.label = QLabel('Channel Clicked')
        self.label.setFocusPolicy(QtCore.Qt.NoFocus)
        self.label.setStyleSheet('border: 1px solid black;')

        self.plots = {}

        if self.templeton and not self.noRecord:
            print(self.waveform_metrics.columns)
            self.waveform_metrics.drop(columns=(['Unnamed: 0', 'cluster_id', 'epoch_name']), inplace=True)
        else:
            self.waveform_metrics.drop(columns=['Unnamed: 0', 'cluster_id', 'epoch_name_quality_metrics', 'epoch_name_waveform_metrics', 'quality'], inplace=True)

        for column in self.waveform_metrics.columns:
            if column == 'peak_channel':
                self.plots['unit_density'] = PlotDisplayItem('unit_density', self.waveform_metrics, self.probeAnnotations, self.mouseID, self.metricsList, label=self.label)
            else:
                self.plots[column] = PlotDisplayItem(column, self.waveform_metrics, self.probeAnnotations, self.mouseID, self.metricsList, label=self.label)

        for plot in self.plots:
            others = [k for k in self.plots if k != plot]
            for other in others:
                self.plots[plot].otherPlots.append(self.plots[other]) # add other plots to list for updating 

        #self.getColorVolume()

        # main layout with region of interest
        self.mainLayout = QVBoxLayout()
        self.imageLayout = QHBoxLayout()

        self.imageLayout.addWidget(self.image)
        self.imageLayout.addWidget(self.imageMask)

        # ui features: probe/metric drop downs, red/green toggle, toggle probe, toggle mask, warp to ccf
        self.probes = self.probeAnnotations['probe_name'].unique()
        self.probeLetters = [s[s.index(' ')+1:][0] for s in self.probes]
        self.probeLetters = sorted(list(set(self.probeLetters)))
        print(self.probeLetters)
        self.probeDropDown = QComboBox()
        for probe in sorted(self.probes):
            self.probeDropDown.addItem(probe)

        self.probeDropDown.setFocusPolicy(QtCore.Qt.NoFocus)
        self.probeDropDown.currentTextChanged.connect(self.displayRegion)
        self.redCheck = QCheckBox('Toggle Red')
        self.redCheck.setFocusPolicy(QtCore.Qt.NoFocus)
        self.redOld = None
        self.isRedChecked = False
        self.redCheck.clicked.connect(self.toggleRed)

        self.greenCheck = QCheckBox('Toggle Green')
        self.greenCheck.setFocusPolicy(QtCore.Qt.NoFocus)
        self.greenOld = None
        self.isGreenChecked = False
        self.greenCheck.clicked.connect(self.toggleGreen)
        
        # rgb sliders to toggle value
        self.redSlider = QRangeSlider(Qt.Horizontal)
        self.redSlider.setFocusPolicy(QtCore.Qt.NoFocus)
        self.redSlider.setMinimum(DEFAULT_COLOR_VALUES[0][0])
        self.redSlider.setMaximum(DEFAULT_COLOR_VALUES[0][1])
        self.redSlider.setValue((DEFAULT_COLOR_VALUES[0][0], DEFAULT_COLOR_VALUES[0][1]))
        self.redSlider.setTickPosition(QSlider.TicksBelow)
        self.redSlider.setTickInterval(10)
        self.redSlider.valueChanged.connect(self.redSliderMoved)

        self.greenSlider = QRangeSlider(Qt.Horizontal)
        self.greenSlider.setFocusPolicy(QtCore.Qt.NoFocus)
        self.greenSlider.setMinimum(DEFAULT_COLOR_VALUES[1][0])
        self.greenSlider.setMaximum(DEFAULT_COLOR_VALUES[1][1])
        self.greenSlider.setValue((DEFAULT_COLOR_VALUES[1][0], DEFAULT_COLOR_VALUES[1][1]))
        self.greenSlider.setTickPosition(QSlider.TicksBelow)
        self.greenSlider.setTickInterval(10)
        self.greenSlider.valueChanged.connect(self.greenSliderMoved)

        # layout for sliders
        self.sliderLayout = QFormLayout()
        self.sliderLayout.addRow('Red Slider', self.redSlider)
        self.sliderLayout.addRow('Green Slider', self.greenSlider)

        """
        self.metrics = QComboBox()
        self.metrics.setFocusPolicy(QtCore.Qt.NoFocus)
        self.metrics.currentTextChanged.connect(self.metricChanged)
        for metric in self.waveform_metrics.columns:
            if metric != 'peak_channel':
                self.metrics.addItem(metric)
        """

        self.toggleProbeButton = QPushButton('Toggle Probe')
        self.toggleProbeButton.clicked.connect(self.toggleProbe)
        self.toggleProbeButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.resetPlotButton = QPushButton('Reset Metric Plot')
        self.resetPlotButton.clicked.connect(self.resetPlotImage)
        self.resetPlotButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.viewButton = QPushButton('View Probe with Selected Metric')
        self.viewButton.clicked.connect(self.displayRegion)
        self.viewButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.calcDistanceButton = QPushButton('Calculate Probe Distance')
        self.calcDistanceButton.clicked.connect(self.calcDistance)
        self.calcDistanceButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.viewCCFButton = QPushButton('Toggle Mask')
        self.viewCCFButton.clicked.connect(self.toggleCCFRegions)
        self.viewCCFButton.setFocusPolicy(QtCore.Qt.NoFocus)
        self.showCCF = False

        self.warpAnchorButton = QPushButton('Warp Selected Anchor')
        self.warpAnchorButton.clicked.connect(self.warpSelectedAnchor)
        self.warpAnchorButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.warpButton = QPushButton('Warp Selected Probe to CCF')
        self.warpButton.clicked.connect(self.warpProbe)
        self.warpButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.warpAllButton = QPushButton('Warp All Probes to CCF')
        self.warpAllButton.clicked.connect(self.warpAllProbes)
        self.warpAllButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.viewWarpedButton = QPushButton('View Warped Channels for Selected Probe')
        self.viewWarpedButton.setFocusPolicy(QtCore.Qt.NoFocus)
        self.viewWarpedButton.clicked.connect(self.viewWarpedChannels)
        self.saveCCFButton = QPushButton('Save channel alignments')
        self.saveCCFButton.setFocusPolicy(QtCore.Qt.NoFocus)
        self.saveCCFButton.clicked.connect(self.saveCCFChannels)

        self.anchorRefineLayout = QVBoxLayout()
        self.anchorRefineLayout.addWidget(self.saveCCFButton)
        """
        self.anchorRefineLayout.addWidget(self.warpAnchorButton)
        self.anchorRefineLayout.addWidget(self.warpButton)
        self.anchorRefineLayout.addWidget(self.warpAllButton)
        self.anchorRefineLayout.addWidget(self.viewWarpedButton)
        """

        self.probeViewLayout = QHBoxLayout()
        self.probeViewLayout.addWidget(self.probeDropDown)
        #self.probeViewLayout.addWidget(self.metrics)
        #self.probeViewLayout.addWidget(self.toggleProbeButton)
        self.probeViewLayout.addWidget(self.viewButton)
        self.probeViewLayout.addWidget(self.resetPlotButton)
        self.probeViewLayout.addWidget(self.toggleProbeButton)
        #self.probeViewLayout.addWidget(self.calcDistanceButton)
        self.probeViewLayout.addWidget(self.viewCCFButton)
        #self.probeViewLayout.addWidget(self.label)
        self.probeViewLayout.addWidget(self.redCheck)
        self.probeViewLayout.addWidget(self.greenCheck)
        self.probeViewLayout.addLayout(self.sliderLayout)

        #self.probeViewLayout.addWidget(self.warpButton)
        self.probeViewLayout.addLayout(self.anchorRefineLayout)
        self.probeViewLayout.setAlignment(QtCore.Qt.AlignTop)
        self.mainLayout.addLayout(self.imageLayout)
        self.mainLayout.addLayout(self.probeViewLayout)

        self.warped = False
        self.readAnchors()

        self.setLayout(self.mainLayout)
        self.showMaximized()

    # helper function to threshold image
    def imageClip(self, color, val):
        if color == 'red':
            array = self.int_arrays['red']
        elif color == 'green':
            array = self.int_arrays['green']
        elif color == 'blue':
            array = self.int_arrays['blue']

        clip = np.clip(array, a_min=val[0], a_max=val[1]) - val[0]
        clip_8 = (clip * 255. / (val[1] - val[0])).astype('uint8')

        return clip_8

    def redSliderMoved(self):
        if not self.isRedChecked:
            self.volArray[:, :, 0] = self.imageClip('red', self.redSlider.value())

            if not self.showMask:
                rotate = np.rot90(self.volArray)
                flip = np.flipud(rotate)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else:
                self.showMaskHelper()

    def greenSliderMoved(self):
        if not self.isGreenChecked:
            self.volArray[:, :, 1] = self.imageClip('green', self.greenSlider.value())

            if not self.showMask:
                rotate = np.rot90(self.volArray)
                flip = np.flipud(rotate)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else:
                self.showMaskHelper()

    # function to bring mask back 
    def showMaskHelper(self):
        im = Image.fromarray(self.volArray)
        overlay = Image.blend(im, self.mask, 0.5)
        #overlay.show()

        draw = ImageDraw.Draw(overlay)

        for label in self.labels_pos:
            positions = self.labels_pos[label]
            centroid = np.array(positions).mean(axis=0)
            draw.text((int(np.round(centroid[1])), int(np.round(centroid[0]))), label, font=self.myFont, fill=(255, 255, 255))

        self.blended = np.array(overlay)
        rot = np.rot90(self.blended)
        flip = np.flipud(rot)
        self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)

    # hides/shows the red channel of the volume
    def toggleRed(self):
        if not self.isRedChecked: # hide red channel
            self.redOld = self.volArray[:, :, 0].copy()
            self.volArray[:, :, 0] = 0
            self.isRedChecked = True

            if not self.showMask: # mask is currently hidden
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else: # bring back mask on to volume slice image
                self.showMaskHelper()
        else: # show red channel
            self.volArray[:, :, 0] = self.redOld.copy()
            self.volArray[:, :, 0] = self.imageClip('red', self.redSlider.value())
            self.isRedChecked = False

            if not self.showMask:
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else:
                self.showMaskHelper()

    # hides/shows the green channel of the volume
    def toggleGreen(self):
        if not self.isGreenChecked: # hide green
            self.greenOld = self.volArray[:, :, 1].copy()
            self.volArray[:, :, 1] = 0
            self.isGreenChecked = True

            if not self.showMask: # mask is currently hidden
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else: # bring back mask on to volume slice image
                self.showMaskHelper()
        else:
            self.volArray[:, :, 1] = self.greenOld.copy()
            self.isGreenChecked = False

            if not self.showMask:
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else:
                self.showMaskHelper()

    # calculates the euclidean distance between 2 consecutive points annotated on the probe
    def calcDistance(self):
        sorted_dist = sorted(self.distPoints, key=lambda x: x[1])
        for i in range(len(sorted_dist) - 1):
            popup = QMessageBox()
            popup.setText('{} Distance in mm {:.2f}'.format(self.probeDropDown.currentText(), np.linalg.norm(sorted_dist[i] - sorted_dist[i + 1]) / 100))
            popup.exec_()

    def resetPlotImage(self):
        self.resetPlot()
        anchor = [self.plots['unit_density'].channels, None, self.anchorPts.copy(), 
                                                  self.pointsAdded.copy()] 

        self.alignments[self.probeDropDown.currentText()] = anchor 

        self.saveAnchor(anchor)
    # resets the plots
    def resetPlot(self):
        if hasattr(self.plots['unit_density'], 'channelsOriginal'):
            self.plots['unit_density'].resetPlot()
        
        #self.plots[self.metrics.currentText().lower()].resetPlot()
        self.showProbe = True
        view = self.image.getView()

        if hasattr(self, 'plotImage') and self.prevProbe != self.probeDropDown.currentText():
            view.removeItem(self.plotImage)

        for item in self.plotImageItems:
            view.removeItem(item)

        for item in self.plotImageTextItems:
            view.removeItem(item)
        
        for item in self.plotImageCCFItems:
            view.removeItem(item)

        for item in self.lineItems:
            view.removeItem(item)

        for item in self.textItems:
            view.removeItem(item)

        for item in self.ccfTextItems:
            view.removeItem(item)

        for item in self.ccfPlotItems:
            view.removeItem(item)

        if len(self.probeItems) > 0 and self.prevProbe != self.probeDropDown.currentText():
            if self.showProbe:
                for item in self.ccfAreaItems:
                    view.removeItem(item)

                for item in self.probeItems:
                    view.removeItem(item)

            self.ccfAreaItems.clear()
            self.probeItems.clear()
        
        self.distPoints.clear()
        self.anchorPts.clear()
        self.lineItems.clear()
        self.textItems.clear()
        self.ccfTextItems.clear()
        self.plotImageItems.clear()
        self.plotImageTextItems.clear()
        self.plotImageCCFItems.clear()
        self.ccfPlotItems.clear()
        self.anchorItems.clear()
        self.anchorPos.clear()
        self.pointsAdded.clear()
        self.oldChannels.clear()

    # hides/shows the probe
    def toggleProbe(self):
        print(self.showProbe)
        view = self.image.getView()
        if self.showProbe:
            for item in self.ccfAreaItems:
                view.removeItem(item)

            for item in self.probeItems:
                view.removeItem(item)

            self.showProbe = False
        else:
            for item in self.ccfAreaItems:
                view.addItem(item)

            for item in self.probeItems:
                view.addItem(item)

            self.showProbe = True
    
    # warps the selected anchor point
    def warpAnchorPoint(self, point):
        # get 3d point in 25 resolution
        ap = int(np.round(point[0] / 2.5))
        dv = int(np.round((point[1]) / 2.5)) # offset for scaling and taking into account image generation 
        ml = int(np.round(point[2] / 2.5))

        df = pd.DataFrame({'AP': [ap], 'DV': [dv], 'ML': [ml]})
        df_result = warp_channels(df, self.field, self.reference, self.probeDropDown.currentText()) # warp anchor

        # display ccf slice
        plane = sitk.GetArrayFromImage(self.reference).T[df_result.AP.values[0], :, :]
        rotate = np.rot90(plane)
        flip = np.flipud(rotate)
        self.imageRefine.setImage(flip, levels=(0, 255), autoRange=False)

        view_ccf = self.imageRefine.getView()

        if hasattr(self, 'lastAnchor'):
            view_ccf.addItem(self.lastAnchor)
        if len(self.anchorItems) > 1:
            view_ccf.removeItem(self.anchorItems.pop(-1))

        pts = [[df_result.ML.values[0], df_result.DV.values[0]]] # plot anchor on ccf

        if len(self.anchorItems) == 0:
            anchor_item = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('blue')), brush=QtGui.QBrush(QColor('blue')), size=5)
        else:
            anchor_item = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('red')), brush=QtGui.QBrush(QColor('red')), size=5)

        view_ccf.addItem(anchor_item)
        self.anchorItems.append(anchor_item)
        self.anchorPos.append([df_result.ML.values[0], df_result.DV.values[0]])
        
        self.warped = True
        self.anchorInd = 0

    # warps the anchor point after it is clicked and displays it
    def clickLine(self, plots, points):
        y = [p[1] for p in self.plots['unit_density'].channels]

        if hasattr(self, 'y_coord'):
            self.oldy = self.y_coord

        for item in self.lineItems: # set all colors to yellow
            item.setPen(QtGui.QPen(QColor('yellow')))
            item.setBrush(QtGui.QBrush(QColor('yellow')))

        self.y_coord = points[0].pos()[1]
        if hasattr(self, 'channel'):
            self.oldChannel = self.channel
        self.channel = y.index(self.y_coord)
        point_ind = self.pointsAdded.index(self.y_coord) # get anchor that was selected
        # changes color to show its been selected
        chosen_item = self.lineItems[point_ind]
        chosen_item.setPen(QtGui.QPen(QColor('white')))
        chosen_item.setBrush(QtGui.QBrush(QColor('white')))

    # helper function to update the anchor
    def updateAnchor(self, anchor, view_ccf):
        # update image
        plane = self.volume[anchor.AP.values[0], :, :]
        rotate = np.rot90(plane)
        flip = np.flipud(rotate)
        self.imageRefine.setImage(flip, levels=(0, 255), autoRange=False)

        # add anchor point
        pts = [[anchor.ML.values[0], anchor.DV.values[0]]]
        anchor_item = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=5)
        view_ccf.addItem(anchor_item)
        self.anchorItems.append(anchor_item)
        self.anchorPos.append([anchor.ML.values[0], anchor.DV.values[0]])

    # moves to prev anchor
    def viewPrevAnchor(self):
        if self.anchorInd > 0:
            view_ccf = self.imageRefine.getView()

            for item in self.anchorItems: # clear items from anchor
                view_ccf.removeItem(item)

            self.anchorItems.clear()
            self.anchorPos.clear()

            y = [p[1] for p in self.channels]
            self.anchorInd -= 1
            anchor = self.df_ccf.loc[self.df_ccf['channel'] == y.index(sorted(self.pointsAdded)[self.anchorInd])] # prevs anchor 

            self.updateAnchor(anchor, view_ccf)

    # moves to the next anchor
    def viewNextAnchor(self):
        if self.anchorInd < len(self.pointsAdded) - 1:
            view_ccf = self.imageRefine.getView()

            for item in self.anchorItems: # clear items from anchor
                view_ccf.removeItem(item)

            self.anchorItems.clear()
            self.anchorPos.clear()

            y = [p[1] for p in self.channels]
            self.anchorInd += 1
            anchor = self.df_ccf.loc[self.df_ccf['channel'] == y.index(sorted(self.pointsAdded)[self.anchorInd])] # next anchor 
            
            self.updateAnchor(anchor, view_ccf)

    # warps all the probes
    def warpAllProbes(self):
        print('Warping Probes', list(self.alignments.keys()))
        for probe in self.alignments:
            coords = self.generateLine(probe)
            channels = self.alignments[probe][0].copy()

            thread = threading.Thread(target=self.warpChannels, args=[channels, probe, coords])
            thread.daemon = True # kill thread once the main is stopped
            thread.start()

        popup = QMessageBox()
        popup.setText('Warping Started. Check console in 2-3 hrs')
        popup.exec_()

    # warps a single probe
    def warpProbe(self):
        print('Warping', self.probeDropDown.currentText())
        coords = self.generateLine(self.probeDropDown.currentText())
        channels = self.alignments[self.probeDropDown.currentText()][0].copy()

        thread = threading.Thread(target=self.warpChannels, args=[channels, self.probeDropDown.currentText(), coords])
        thread.daemon = True # kill thread once the main is stopped
        thread.start()

        popup = QMessageBox()
        popup.setText('Warping Started. Check console in about 5-10 mins')
        popup.exec_()

    # shows the warping of the channels
    def viewWarpedChannels(self):
        probe_name = self.probeDropDown.currentText()
        self.df_ccf = pd.read_csv(os.path.join(self.storageDirectory, '{}_channels_{}_warped.csv'.format(probe_name.replace(' ', '_'), self.mouseID)))
        #self.channels = self.plots['unit_density'].channels
        y = [p[1] for p in self.plots['unit_density'].channels]
        for point in self.pointsAdded: # plot anchors
            ind = y.index(point)
            row = self.df_ccf.loc[self.df_ccf['channel'] == ind]
           
            plt.scatter(row.ML, row.DV, s=15, c='orange')

        grouped = self.df_ccf.groupby('channel_areas').mean()
        
        for index, row in grouped.iterrows(): # plot areas
            if index != 'N/A':
                plt.text(row.ML, row.DV, index, color='white')

        # plot color scheme of areas on image
        prev_area = ''
        view = self.image.getView()
        color = 'light blue'
        areas_seen = set()
        na = 0
        for index, row in self.df_ccf.reindex().sort_index(ascending=False).iterrows():
            area = row.channel_areas

            if not pd.isna(area):
                if prev_area != area: # if moved to a new area
                    if color == 'cyan':
                        color = 'pink'
                    else:
                        color = 'cyan'

                item = pg.ScatterPlotItem(pos=[[80, self.plots['unit_density'].channels[index][1]]], pen=QtGui.QPen(QColor(color)), brush=QtGui.QBrush(QColor(color)), size=5)
                self.ccfPlotItems.append(item)
                view.addItem(item)
                na += 1
                if area not in areas_seen: 
                    areas_seen.add(area)
                    text = pg.TextItem(area, color=color)
                    text.setPos(100, self.plots['unit_density'].channels[index][1])
                    self.ccfTextItems.append(text)
                    view.addItem(text)
                elif area in areas_seen and prev_area != area: # new area but has been seen already
                    text = pg.TextItem(area, color=color)
                    text.setPos(100, self.plots['unit_density'].channels[index][1])
                    self.ccfTextItems.append(text)
                    view.addItem(text)

                prev_area = area
            elif pd.isna(area): # N/A area
                if '%s%d' %('N/A', na) not in areas_seen:
                    areas_seen.add('%s%d' %('N/A', na))
                    text = pg.TextItem('N/A')
                    text.setPos(100, self.plots['unit_density'].channels[index][1])
                    self.ccfTextItems.append(text)
                    view.addItem(text)

                item = pg.ScatterPlotItem(pos=[[80, self.plots['unit_density'].channels[index][1]]], pen=QtGui.QPen(QColor('white')), brush=QtGui.QBrush(QColor('white')), size=5)
                view.addItem(item)
                self.ccfPlotItems.append(item)
                prev_area = 'N/A'

        view.removeItem(self.plItem)
        self.showProbe = False
        plt.imshow(sitk.GetArrayFromImage(self.reference).T[int(grouped.AP.mean()), :, :])
        plt.show()

    # closest value helper
    def closestValue(self, input_list, input_value):
      arr = np.asarray(input_list)
      i = (np.abs(arr - input_value)).argmin()
 
      return arr[i]

    # saves the ccf channels
    def saveCCFChannels(self):
        probe_name = self.probeDropDown.currentText()
        y = [p[1] for p in self.plots['unit_density'].channels]
        dict_final = {'AP': [], 'DV': [], 'ML': [], 'region': [], 'channel': []}
        for i in range(384):
            y_coord = int(np.round(y[i])) + 85
            key = self.closestValue(list(self.ccfAreas.keys()), y_coord)

            ccf_region_coord = self.ccfAreas[key]
            ap = ccf_region_coord[0][0]
            dv = ccf_region_coord[0][1]
            ml = ccf_region_coord[0][2]
            region = ccf_region_coord[1]

            dict_final['AP'].append(ap)
            dict_final['DV'].append(dv)
            dict_final['ML'].append(ml)
            dict_final['channel'].append(i)
            dict_final['region'].append(region)

        df = pd.DataFrame(dict_final)
        df.to_csv(os.path.join(self.storageDirectory, '{}_channels_{}_warped.csv'.format(probe_name.replace(' ', '_'), self.mouseID)), index=False)
        popup = QMessageBox()
        popup.setText('Channel regions saved')
        popup.exec_()

    # warps the channels to the ccf
    def warpChannels(self, channels, probe_name, coords, view=False):
        print('Arjun you are doing great!!!')
        values = list(self.acrnm_map.values())
        keys = list(self.acrnm_map.keys())
       
        unique = set()

        channel_map = {}
        out_of_brain = set()
        for i in range(384):
            channel_dict = {'AP': [], 'DV': [], 'ML': []} # dataframe for warping
            channel = channels[i]
            y_coord = int(np.round(channel[1])) + 85
            if y_coord in coords:
                coord = coords[y_coord] # get the 3d coordinate at that point on the probe track
                ap = int(np.round(coord[0] / 2.5))
                dv = int(np.round((coord[1]) / 2.5)) # offset for scaling and taking into account image generation 
                ml = int(np.round((coord[2]) / 2.5))

                if (ap, dv, ml) not in unique and dv > 0:
                    unique.add((ap, dv, ml))

                    channel_dict['AP'].append(ap)
                    channel_dict['DV'].append(dv)
                    channel_dict['ML'].append(ml)
        
                    df_channel = pd.DataFrame(channel_dict)
        
                    # warp using deformation field and reference
                    df_final = warp_channels(df_channel, self.field, self.reference, probe_name)
                    channel_map[(ap, dv, ml)] = (df_final['AP'].values[0], df_final['DV'].values[0], df_final['ML'].values[0])
                elif dv < 0:
                    out_of_brain.add((ap, dv, ml))

        ccf = {'AP': [], 'DV': [], 'ML': []} # dictionary for ccf coords

        for i in range(0, 384):
            channel = channels[i]
            y_coord = int(np.round(channel[1])) + 85

            if y_coord in coords:
                coord = coords[y_coord] # get the 3d coordinate at that point on the probe track
                ap = int(np.round(coord[0] / 2.5))
                dv = int(np.round((coord[1]) / 2.5)) # offset for scaling and taking into account image generation 
                ml = int(np.round((coord[2]) / 2.5))

                if dv > 0:
                    points = channel_map[(ap, dv, ml)] # get corresponding ccf coord
                    #print(points)
                    ccf['AP'].append(points[0])
                    ccf['DV'].append(points[1])
                    ccf['ML'].append(points[2])
                else: # out of brain
                    ccf['AP'].append(-1)
                    ccf['DV'].append(-1)
                    ccf['ML'].append(-1)

        df_ccf = pd.DataFrame(ccf)
        channels_num = list(range(384))
        probe = [probe_name for i in range(len(ccf['AP']))]
        df_ccf['probe_name'] = probe
        # get areas from ccf
        channel_areas = []

        # get brain areas corresponding to coords
        for index, row in  df_ccf.iterrows():
            if row.AP > 0:
                struct = self.anno[row.AP, row.DV, row.ML]
                if struct in values:
                    ind = values.index(struct)
                    key = keys[ind]

                    if not key[0].islower():
                        channel_areas.append(key)
                    else:
                        channel_areas.append('N/A')
                else:
                    channel_areas.append('N/A')
            else:
                channel_areas.append('out of brain')

        df_ccf['channel_areas'] = channel_areas
        df_ccf['channel'] = channels_num
        #print(df_channel)
        
        df_ccf.to_csv(os.path.join(self.storageDirectory, '{}_channels_{}_warped.csv'.format(probe_name.replace(' ', '_'), self.mouseID)), index=False)

        if view:
            self.viewWarpedChannels()

        print('Finsihed Warping {}.'.format(probe_name))
        print()

    # function called when the drop down for the metric changes
    # updates the plots to get the metrics file for the corresponding probe and day
    def metricChanged(self):
        print('Metric changed')
        probe = self.probeDropDown.currentText()
        metric = self.metrics.currentText()
        view = self.image.getView()

        probe_let_num = probe[probe.index(' ')+1:]

        #key = self.days[int(probe_let_num[1])]
        paths = self.waveMetricsPath[int(probe_let_num[1])]

        self.path = ''
        for p in paths:
            if 'probe' + probe_let_num[0] in p:
                self.path = p
                break
            elif self.templeton and 'Probe' + probe_let_num[0] in p:
                self.path = p
                break

        if self.path != '':
            self.plots[metric].updateMetrics(self.path, self.templeton)
        #self.plots[metric].generateMetrics(metric.lower())

        if hasattr(self, 'linepts'):
            self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues, keep_y=True, old_y=self.plots[self.oldMetric].channels, points_added=self.pointsAdded)

        if hasattr(self, 'oldMetric'):
            view.removeItem(self.plots[self.oldMetric].channelsPlot)
            view.addItem(self.plots[metric].channelsPlot)

        self.oldMetric = metric

    def readAnchors(self):
        for probe in self.probes:
            if os.path.exists(os.path.join(self.storageDirectory, 'anchors', '{}_anchors.pickle'.format(probe.replace(' ', '_')))):
                if probe not in self.alignments:
                    with open(os.path.join(self.storageDirectory, 'anchors', '{}_anchors.pickle'.format(probe.replace(' ', '_'))), 'rb') as f:
                        anchor = pickle.load(f)

                    self.alignments[probe] = anchor
    # displays the region of interest surrounding the probe track
    def displayRegion(self):
        probe = self.probeDropDown.currentText()
        #self.qc_plot = qcChecker(self.mouseID, probe[probe.index(' ') + 1:])

        #metric = self.metrics.currentText()
        self.path = ''
        view = self.image.getView()

        if os.path.exists(os.path.join(self.storageDirectory, 'anchors', '{}_anchors.pickle'.format(probe.replace(' ', '_')))):
            if probe not in self.alignments:
                with open(os.path.join(self.storageDirectory, 'anchors', '{}_anchors.pickle'.format(probe.replace(' ', '_'))), 'rb') as f:
                    anchor = pickle.load(f)

                self.alignments[probe] = anchor

            # get metrics path for this probe and day
            probe_let_num = probe[probe.index(' ')+1:]

            #key = self.days[int(probe_let_num[1])]
            if self.templeton:
                key = self.days[int(probe_let_num[1]) - 1]
            try:
                if self.templeton:
                    paths = self.waveMetricsPath[key]
                else:
                    paths = self.waveMetricsPath[int(probe_let_num[1])]

                for p in paths:
                    if 'probe' + probe_let_num[0] in p:
                        self.path = p
                        break
                    elif self.templeton and 'Probe' + probe_let_num[0] in p:
                        self.path = p
                        break

                self.plots['unit_density'].updateMetrics(self.path, self.templeton)
                #self.plots['unit_density'].updateDisplay(probe, self.linepts, self.intensityValues) # update display since no existing alignment has been done so far
                #self.plots[metric].updateMetrics(self.path, self.templeton)
                self.resetPlot()
                self.updateDisplay(probe, restore=True)
                #self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues, keep_y=True, old_y=self.alignments[probe][1], points_added=self.alignments[probe][-1])

                self.prevProbe = probe
                #self.oldMetric = metric
            except (FileNotFoundError, KeyError):
                view.removeItem(self.plots['unit_density'].channelsPlot)
                #view.removeItem(self.plots[metric].channelsPlot)
                self.resetPlot()
                self.updateDisplay(probe)
                popup = QMessageBox()
                popup.setText('Couldn\'t find metrics.csv for {}'.format(probe))
                popup.exec_()
        else:
            # get metrics path for this probe and day
            probe_let_num = probe[probe.index(' ')+1:]

            try:
                #key = self.days[int(probe_let_num[1])]
                if self.templeton:
                    key = self.days[int(probe_let_num[1]) - 1]
           
                if self.templeton:
                    paths = self.waveMetricsPath[key]
                else:
                    paths = self.waveMetricsPath[int(probe_let_num[1])]

                for p in paths:
                    if 'probe' + probe_let_num[0] in p:
                        self.path = p
                        break
                    elif self.templeton and 'Probe' + probe_let_num[0] in p:
                        self.path = p
                        break

                if self.prevProbe == '' or self.prevProbe != probe: # new alignment
                    if self.prevProbe != '' and self.prevProbe != probe:
                        if self.path == '':
                            view.removeItem(self.plots['unit_density'].channelsPlot)
                            #view.removeItem(self.plots[metric].channelsPlot)

                        self.resetPlot()
                        self.updateDisplay(probe)
                    else: # initial display when nothing has been done
                        self.plots['unit_density'].updateMetrics(self.path, self.templeton)
                        self.updateDisplay(probe)

                    if self.path != '':
                        self.plots['unit_density'].updateMetrics(self.path, self.templeton)
                        self.plots['unit_density'].updateDisplay(probe, self.linepts, self.intensityValues) # update display since no existing alignment has been done so far

                        #self.plots[metric].updateMetrics(self.path, self.templeton)
                        #self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues)

                        self.prevProbe = probe
                        #self.oldMetric = metric
                    else:
                        popup = QMessageBox()
                        popup.setText('Couldn\'t find metrics.csv for {}'.format(probe))
                        popup.exec_()
            except (FileNotFoundError, IndexError, KeyError): # metrics file not found
                view.removeItem(self.plots['unit_density'].channelsPlot)
                #view.removeItem(self.plots[metric].channelsPlot)
                self.resetPlot()
                self.updateDisplay(probe)
                popup = QMessageBox()
                popup.setText('Couldn\'t find metrics.csv for {}'.format(probe))
                popup.exec_()

    # helper function to remove the alignment line
    def removeLineHelper(self, y_coord, ind):
        view = self.image.getView()
        item = self.lineItems[ind]
        view.removeItem(item)
        self.lineItems.remove(item)
        self.pointsAdded.remove(y_coord)

        self.plots['unit_density'].channelsPlot.setData(pos=np.array(self.plots['unit_density'].channels, dtype=float), adj=np.array(self.plots['unit_density'].adj, dtype=int))
        self.plots[self.metrics.currentText().lower()].channelsPlot.setData(pos=np.array(self.plots[self.metrics.currentText().lower()].channels, dtype=float), 
                                                                            adj=np.array(self.plots[self.metrics.currentText().lower()].adj, dtype=int))
        self.plots['unit_density'].linearSpacePoints(self.pointsAdded)
        self.plots[self.metrics.currentText().lower()].linearSpacePoints(self.pointsAdded)


    # remove the line clicked on from the refinement display
    def removeLine(self, plots, points):
        y_coord = points[0].pos()[1]
        ind = self.pointsAdded.index(y_coord)
        self.removeLineHelper(y_coord, ind)

    # updates the other plots to align it with the plot that was clicked on (the unit density plot in this case)
    def updateAfterClick(self, other_plot, line_point, channel, pointsAdded):
        #print('Other plot clicked', other_plot.measurement)
        other_plot.channelsPlot.scatterPoint = other_plot.channels[channel]
        other_plot.oldChannels.append(other_plot.channels)
        flag = True

        if line_point[1] != other_plot.channelsPlot.scatterPoint[1]:
            if line_point[1] < other_plot.channelsPlot.scatterPoint[1]:
                #print('lower')
                diff = abs(line_point[1] - other_plot.channelsPlot.scatterPoint[1])
                newPoints = []

                if len(pointsAdded) == 0:
                    for p in other_plot.channels:
                        newPoints.append([p[0], p[1] - diff])
                else:
                    srt = sorted(pointsAdded, reverse=True)
                    greater = [t for t in pointsAdded if t > line_point[1]]
                    less = [t for t in pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 exisiting, don't shift, just linearly interpolate
                        flag = other_plot.linspaceAbove(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0: # new line is below exiting, so shift
                        for p in other_plot.channels:
                            if p[1] <= other_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else: # new alignment is above, so interpolate and shift
                        flag = other_plot.linspaceAbove(srt, channel, diff, flag)

                        for i in range(len(other_plot.channels)):
                            p = other_plot.channels[i]
                            if p[1] > other_plot.channelsPlot.scatterPoint[1]:
                                other_plot.channels[i] = [p[0], p[1] - diff]
            elif line_point[1] > other_plot.channelsPlot.scatterPoint[1]:
                diff = line_point[1] - other_plot.channelsPlot.scatterPoint[1]
                newPoints = []

                if len(pointsAdded) == 0:
                    for p in other_plot.channels:
                        newPoints.append([p[0], p[1] + diff])
                else:
                    srt = sorted(pointsAdded)
                    greater = [t for t in pointsAdded if t > line_point[1]]
                    less = [t for t in pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0:
                        flag = other_plot.linspaceBelow(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0:
                        flag = other_plot.linspaceBelow(srt, channel, diff, flag)

                        for i in range(len(other_plot.channels)):
                            p = other_plot.channels[i]
                            if p[1] < other_plot.channelsPlot.scatterPoint[1]:
                                other_plot.channels[i] = [p[0], p[1] + diff]
                    else:
                        for p in other_plot.channels:
                            if p[1] >= other_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
            
            if flag:
                other_plot.channels = newPoints
            #if len(self.pointsAdded) == 1:
            other_plot.channelsPlot.setData(pos=np.array(other_plot.channels, dtype=float), adj=np.array(other_plot.adj, dtype=int))

        #other_plot.linearSpacePoints(pointsAdded)

    # pefroms alignment on the unit density plot
    def onClickProbeHelper(self, click_plot, line_point, is_unit=True, channel=-1, scatter_point=None, color='yellow', arrow_key=False):
        if channel == -1:
            channel = click_plot.channels.index([click_plot.channelsPlot.scatterPoint[0], click_plot.channelsPlot.scatterPoint[1]]) # get index of point clicked on in unit density
        else:
            click_plot.channelsPlot.scatterPoint = scatter_point

        click_plot.oldChannels.append(click_plot.channels)
        flag = True
        print('Clicked', click_plot.channelsPlot.scatterPoint[1])
        if line_point[1] != click_plot.channelsPlot.scatterPoint[1]:
            if line_point[1] < click_plot.channelsPlot.scatterPoint[1]: # alignment is above currently clicked channel
                #print('lower')
                diff = abs(line_point[1] - click_plot.channelsPlot.scatterPoint[1])
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in click_plot.channels:
                        newPoints.append([p[0], p[1] - diff])
                else:
                    srt = sorted(self.pointsAdded, reverse=True)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 exisiting, don't shift, just linearly interpolate
                        flag = click_plot.linspaceAbove(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0: # new line is below existing, so shift
                        for p in click_plot.channels:
                            if p[1] <= click_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else: # new line is above so linearly interpolate and shift
                        flag = click_plot.linspaceAbove(srt, channel, diff, flag)

                        for i in range(len(click_plot.channels)):
                            p = click_plot.channels[i]
                            if p[1] > click_plot.channelsPlot.scatterPoint[1]:
                                click_plot.channels[i] = [p[0], p[1] - diff]
            elif line_point[1] > click_plot.channelsPlot.scatterPoint[1]: # alignment is below currently clicked channel
                diff = line_point[1] - click_plot.channelsPlot.scatterPoint[1]
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in click_plot.channels:
                        newPoints.append([p[0], p[1] + diff])
                else:
                    srt = sorted(self.pointsAdded)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 existing, don't shift, just linearly interpolate
                        flag = click_plot.linspaceBelow(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0: # new alignment is above, so linearly interpolate and shift
                        flag = click_plot.linspaceBelow(srt, channel, diff, flag)

                        for i in range(len(click_plot.channels)):
                            p = click_plot.channels[i]
                            if p[1] < click_plot.channelsPlot.scatterPoint[1]:
                                click_plot.channels[i] = [p[0], p[1] + diff]
                    else: # new alignment is below, so just shift
                        for p in click_plot.channels:
                            if p[1] >= click_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
                        
            if flag:
                click_plot.channels = newPoints
            #if len(self.pointsAdded) == 1:
            click_plot.channelsPlot.setData(pos=np.array(click_plot.channels, dtype=float), adj=np.array(click_plot.adj, dtype=int))
           
            pts = [[t, line_point[1]] for t in range(int(click_plot.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            self.anchorPts.append(pts)
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor(color)), brush=QtGui.QBrush(QColor(color)), size=2)
            text = pg.TextItem(str(channel))
            text.setPos(80, line_point[1])
            self.textItems.append(text)
            self.lineItems.append(h_line)
        else: # alignment is same as y coord of the channel clicked on
            pts = [[t, line_point[1]] for t in range(int(click_plot.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            self.anchorPts.append(pts)
            diff = 0
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor(color)), brush=QtGui.QBrush(QColor(color)), size=2)
            self.lineItems.append(h_line)
            text = pg.TextItem(str(channel))
            text.setPos(80, line_point[1])
            self.textItems.append(text)

        # update other plots to match alignment of unit density
        """
        if is_unit:
            self.updateAfterClick(self.plots[self.metrics.currentText().lower()], line_point, channel, self.pointsAdded)
        else:
            self.updateAfterClick(self.plots['unit_density'], line_point, channel, self.pointsAdded)
        """
        self.pointsAdded.append(line_point[1])
        click_plot.linearSpacePoints(self.pointsAdded)
        view = self.image.getView()

        if not arrow_key:
            plot_item = pg.ScatterPlotItem(pos=[[i, self.plots['unit_density'].ogChannelsShift[channel][1]] for i in range(-700, 
                                                                                                       int(self.plots['unit_density'].ogChannelsShift[channel][0]))], 
                                                      pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=3)

            view.addItem(plot_item)
            self.plotImageItems.append(plot_item)

            y = [point[1] for point in self.plots['unit_density'].channels]
            self.removeImagePlotCCFAreas()
            self.displayImageCCFAreas(y)
        """
        if is_unit:
            self.plots[self.metrics.currentText().lower()].linearSpacePoints(self.pointsAdded)
        else:
            self.plots['unit_density'].linearSpacePoints(self.pointsAdded)
        """
        view.addItem(h_line)
        view.addItem(text)
        #h_line.setClickable(True)
        h_line.sigClicked.connect(self.clickLine)

        anchor = [self.plots['unit_density'].channels, None, self.anchorPts.copy(), 
                                              self.pointsAdded.copy()] 

        self.alignments[self.probeDropDown.currentText()] = anchor
        self.saveAnchor(anchor)

    # saves the anchors for the probe
    def saveAnchor(self, anchor):
        probe = self.probeDropDown.currentText().replace(' ', '_')

        with open(os.path.join(self.storageDirectory, 'anchors', '{}_anchors.pickle'.format(probe)), 'wb') as f:
            pickle.dump(anchor, f)

    def onclickProbe(self, plot, points):
        if self.plots['unit_density'].channelsPlot.clicked:
            line_point = points[0].pos()
            self.onClickProbeHelper(self.plots['unit_density'], line_point, is_unit=True)
        else:
            line_point = points[0].pos()
            self.onClickProbeHelper(self.plots[self.metrics.currentText().lower()], line_point, is_unit=False)

        self.plots['unit_density'].channelsPlot.clicked = False
        #self.plots[self.metrics.currentText().lower()].channelsPlot.clicked = False
    
    def refineSpline(self, new_y, dist):
        y = [p[1] for p in self.channels]
        current_channel = y.index(sorted(self.pointsAdded)[self.anchorInd])

        if self.anchorInd + 1 < len(self.pointsAdded):
            next_channel = y.index(sorted(self.pointsAdded)[self.anchorInd + 1]) # next anchor channel
                        
            # get distance between old and new
            orig_dist = self.splineWarped.loc[self.splineWarped['channel'] == next_channel].DV.values[0] - \
            self.splineWarped.loc[self.splineWarped['channel'] == current_channel].DV.values[0]

            new_dist = self.splineWarped.loc[self.splineWarped['channel'] == next_channel].DV.values[0] - new_y

            scale = (orig_dist - new_dist) # scale
            scale /= 10
            print('Scale', scale)
            start = self.splineWarped.loc[self.splineWarped['channel'] == next_channel].index[0]
            end = self.splineWarped[self.splineWarped['channel'] == current_channel].index[0]
            print('Before', self.splineWarped.iloc[start + 1: end + 1])
            self.splineWarped.loc[start + 1: end + 1, 'DV'] += dist # apply shift
            #self.splineWarped.loc[start + 1: end + 1, 'DV'] *= (1 + scale) # apply scale
            result = self.spline(self.splineWarped['DV'].to_numpy())
            result = np.fliplr(result)
            self.splineWarped.loc[start + 1: end + 1, 'AP'] = pd.Series(result[0], index=list(range(384)))
            self.splineWarped.loc[start + 1: end + 1, 'ML'] = pd.Series(result[1], index=list(range(384)))
            print('After', self.splineWarped.iloc[start + 1: end + 1])

    # helper function for refining the anchor based on arrow key
    def refineAnchorHelper(self, orientation, dist):
        view = self.imageRefine.getView()
        y = [p[1] for p in self.channels]

        if len(self.anchorPos) > 1:
            popped = self.anchorPos.pop(-1)
            new_x = popped[0]
            new_y = popped[1]

            view.removeItem(self.anchorItems.pop(-1))

            if orientation == 'v': # up or down
                new_y +=  dist
                self.refineSpline(new_y, dist)
            else:
                new_x += dist

            item = pg.ScatterPlotItem(pos=[[new_x, new_y]], pen=QtGui.QPen(QColor('orange')), brush=QtGui.QBrush(QColor('orange')), size=5)
            view.addItem(item)
            self.anchorPos.append([new_x, new_y])
            self.anchorItems.append(item)
        else:
            original = self.anchorPos[0]
            new_x = original[0]
            new_y = original[1]

            if orientation == 'v': # up or down
                new_y +=  dist
                self.refineSpline(new_y, dist)
            else:
                new_x += dist

            item = pg.ScatterPlotItem(pos=[[new_x, new_y]], pen=QtGui.QPen(QColor('orange')), brush=QtGui.QBrush(QColor('orange')), size=5)
            view.addItem(item)
            self.anchorPos.append([new_x, new_y])
            self.anchorItems.append(item)

    # helper function when using keyboard to refine anchors
    # returns the channel that the alignment corresponds to
    def removeAnchorHelper(self):
        print('Helper', self.y_coord)
        y = [p[1] for p in self.plots['unit_density'].channels]
        channel = y.index(self.y_coord)
        view = self.image.getView()
        ind = self.pointsAdded.index(self.y_coord)
        view.removeItem(self.lineItems.pop(ind))
        view.removeItem(self.textItems.pop(ind))
        self.pointsAdded.pop(ind)
        self.anchorPts.pop(ind)

        return channel
    
    # warps the selected anchor
    def warpSelectedAnchor(self):
        y = [p[1] for p in self.plots['unit_density'].channels]
        point = self.coords[self.y_coord + 85]
        self.warpAnchorPoint(point)

    # clears the ccf image
    def clearCCFImage(self, clear_anchor=False):
        view_ccf = self.imageRefine.getView()
        for item in self.anchorItems:
            view_ccf.removeItem(item)

        for item in self.anchorText:
            view_ccf.removeItem(item)

        if clear_anchor:
           self.anchorText.clear()

        self.anchorItems.clear()
        self.imageRefine.getImageItem().clear()

    # performs action based on what key is hit
    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Left: # delete last anchor done
            popped = self.lineItems.pop(-1) # last alignment done
            popped_text = self.textItems.pop(-1)
            popped_img_plot = self.plotImageItems.pop(-1)

            view = self.image.getView()
            view.removeItem(popped)
            view.removeItem(popped_text)
            view.removeItem(popped_img_plot)

            self.anchorPts.pop(-1)
            self.pointsAdded.pop(-1)
            self.plots['unit_density'].linearSpacePoints(self.pointsAdded)
            #self.plots[self.metrics.currentText().lower()].linearSpacePoints(self.pointsAdded)

            anchor = [self.plots['unit_density'].channels, None, self.anchorPts.copy(), 
                                              self.pointsAdded.copy()] 
            self.alignments[self.probeDropDown.currentText()] = anchor
            self.saveAnchor(anchor)

            y = [point[1] for point in self.plots['unit_density'].channels]
            self.removeImagePlotCCFAreas()
            self.displayImageCCFAreas(y)
        elif event.key() == Qt.Key_Delete: # delete selected anchor
            if self.y_coord in self.pointsAdded:
                ind = self.pointsAdded.index(self.y_coord)
                print('ImageItems', len(self.plotImageItems))
                popped = self.lineItems.pop(ind) # last alignment done
                popped_text = self.textItems.pop(ind)
                popped_img_plot = self.plotImageItems.pop(ind)

                view = self.image.getView()
                view.removeItem(popped)
                view.removeItem(popped_text)
                view.removeItem(popped_img_plot)
                self.anchorPts.pop(ind)
                self.pointsAdded.pop(ind)
                self.plots['unit_density'].linearSpacePoints(self.pointsAdded)
                #self.plots[self.metrics.currentText().lower()].linearSpacePoints(self.pointsAdded)

                anchor = [self.plots['unit_density'].channels, None, self.anchorPts.copy(), 
                                                  self.pointsAdded.copy()] 
                self.alignments[self.probeDropDown.currentText()] = anchor
                self.saveAnchor(anchor)

                y = [point[1] for point in self.plots['unit_density'].channels]
                self.removeImagePlotCCFAreas()
                self.displayImageCCFAreas(y)
        elif event.key() == Qt.Key_Down: # refine anchor down
            #self.refineAnchorHelper('v', 1)
            channel = self.removeAnchorHelper()
            self.onClickProbeHelper(self.plots['unit_density'], [80, self.y_coord + 1], channel=channel, scatter_point=self.plots['unit_density'].channels[channel], 
                                    color='white',  arrow_key=True)
            self.y_coord += 1
            print('Key', self.y_coord)

            y = [point[1] for point in self.plots['unit_density'].channels]
            self.removeImagePlotCCFAreas()
            self.displayImageCCFAreas(y)
        elif event.key() == Qt.Key_Up: # refine anchor up
            #self.refineAnchorHelper('v', -1)
            channel = self.removeAnchorHelper()
            self.onClickProbeHelper(self.plots['unit_density'], [80, self.y_coord - 1], channel=channel, scatter_point=self.plots['unit_density'].channels[channel], 
                                    color='white', arrow_key=True)
            self.y_coord -= 1

            y = [point[1] for point in self.plots['unit_density'].channels]
            self.removeImagePlotCCFAreas()
            self.displayImageCCFAreas(y)

    # toggles the mask 
    def toggleCCFRegions(self):
        if self.showMask:
            rot = np.rot90(self.volArray)
            flip = np.flipud(rot)
            self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            self.showMask = False
        else:
            self.showMaskHelper()
            self.showMask = True

    # generates the line based on the annotated points
    def generateLine(self, probe):
        x = self.probeAnnotations[self.probeAnnotations.probe_name == probe].ML 
        y = self.probeAnnotations[self.probeAnnotations.probe_name == probe].DV 
        
        z = self.probeAnnotations[self.probeAnnotations.probe_name == probe].AP

        # get trajectory
        if len(z) > 0: # 
            data = np.vstack((z,y,x)).T
            datamean = data.mean(axis=0)

            D = data - datamean
            m1 = np.min(D[:,1]) * 2
            m2 = np.max(D[:,1]) * 2
            uu,dd,vv = np.linalg.svd(D)

            linepts = vv[0] * np.mgrid[-530:600:1][:,np.newaxis]
            linepts += datamean
            
            if linepts[-1,1] - linepts[0,1] < 0:
                linepts = np.flipud(linepts)
            
        print(linepts.shape[0]) 
        intensity_values_red = np.zeros((linepts.shape[0],160)) # dummy array to be used later on for class fields below
        self.coords = {}
        self.linepts = linepts
        self.intensityValues = intensity_values_red

        ap_scale = 1
        dv_scale = 1/0.94
        lr_scale = 1
 
        # get vectors in xyz
        self.dz = (linepts[-1, 0] - linepts[0, 0])
        self.dy = (linepts[-1, 1] - linepts[0, 1])
        self.dx = (linepts[-1, 2] - linepts[0, 2])
        self.vector = np.array([self.dz, self.dy, self.dx])
            
        # apply scale
        self.dzScale = self.dz * ap_scale
        self.dyScale = self.dy * dv_scale
        self.dxScale = self.dx * lr_scale
        self.vectorScale = np.array([self.dzScale, self.dyScale, self.dxScale])

        # get scale factor
        self.scale = np.linalg.norm(self.vectorScale) / np.linalg.norm(self.vector)

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2]) # dictionary to store 3d coord at probe y coord

        return self.coords

    # displays the ccf areas for each point on the line trajectory
    def displayCCFAreas(self):
        color = 'cyan'
        view = self.image.getView()
        prev_area = ''
        self.allAreas = {}

        for i in range(100, self.linepts.shape[0]):
            if i + 85 in self.ccfAreas: # check if area in dictionary
                area = self.ccfAreas[i + 85][1]

                if area != prev_area: # new area seen
                    if area == 'N/A':
                        color = 'white'
                    else:
                        if color == 'cyan':
                            color = 'pink'
                        else:
                            color = 'cyan'

                    text = pg.TextItem(area, color=color)
                    text.setPos(100, i)
                    view.addItem(text)
                    self.ccfAreaItems.append(text)

                    prev_area = area

            item = pg.ScatterPlotItem(pos=[[80, i]], pen=QtGui.QPen(QColor(color)), brush=QtGui.QBrush(QColor(color)), size=5)
            item.sigClicked.connect(self.onclickProbe)
            self.probeItems.append(item)
            view.addItem(item)

    def removeImagePlotCCFAreas(self):
        view = self.image.getView()
        for item in self.plotImageCCFItems:
            view.removeItem(item)

        for item in self.plotImageTextItems:
            view.removeItem(item)

        self.plotImageCCFItems.clear()
        self.plotImageTextItems.clear()

    # displays the ccf areas for the channels with the image plot
    def displayImageCCFAreas(self, y: list):
        color = 'cyan'
        view = self.image.getView()
        prev_area = ''
        initial_pos = self.plots['unit_density'].ogChannelsShift[0]
        for i in range(len(y)):
            y_coord = int(np.round(y[384 - i - 1])) + 85
            key = self.closestValue(list(self.ccfAreas.keys()), y_coord)
            area = self.ccfAreas[key][1]

            channel_pos = self.plots['unit_density'].ogChannelsShift[384 - i - 1]

            if area != prev_area: # new area seen
                if area == 'N/A':
                    color = 'white'
                else:
                    if color == 'cyan':
                        color = 'pink'
                    else:
                        color = 'cyan'
            
                text = pg.TextItem(area, color=color)
                text.setPos(initial_pos[0] + 60, channel_pos[1])
                self.plotImageTextItems.append(text)
                view.addItem(text)
                prev_area = area

            item = pg.ScatterPlotItem(pos=[[initial_pos[0] + 50, channel_pos[1]]], pen=QtGui.QPen(QColor(color)), brush=QtGui.QBrush(QColor(color)), size=5)
            self.plotImageCCFItems.append(item)
            view.addItem(item)

    def updatePlotImage(self):
        view = self.image.getView()
        transform = [5., 0., 0., 0., 5., 0., -700.,
                        self.plots['unit_density'].channels[-1][1], 1.]
        self.plotImage.setTransform(QtGui.QTransform(*transform))
        #view.addItem(self.plotImage)

    # displays the plot image: correlation, etc.
    # data: dict with plot info
    def displayPlotImage(self):
        probe = self.probeDropDown.currentText()
        probe = probe[probe.index(' ')+1:]
        with open(pathlib.Path(self.storageDirectory, 'image_plots', '{}_corr.pickle'.format(probe)), 'rb') as f:
            data = pickle.load(f)

        view = self.image.getView()
        self.plotImage = pg.ImageItem()
        #print(data)
        rot = np.rot90(np.rot90(data['img']))
        flip = np.flipud(rot)
        #print(flip.shape)
        self.plotImage.setImage(flip)
        transform = [data['scale'][0], 0., 0., 0., data['scale'][1], 0., -700.,
                        self.plots['unit_density'].ogChannelsShift[-1][1], 1.]
        self.plotImage.setTransform(QtGui.QTransform(*transform))
        cmap = data.get('cmap', [])
        
        if cmap:
            color_bar = cb.ColorBar(data['cmap'])
            lut = color_bar.getColourMap()
            self.plotImage.setLookupTable(lut)
            self.plotImage.setLevels((data['levels'][0], data['levels'][1]))
            #cbar = color_bar.makeColourBar(20, 5, self.fig_img_cb, min=data['levels'][0],
                                               #max=data['levels'][1], label=data['title'])
        else:
            self.plotImage.setLevels((1, 0))

        view.addItem(self.plotImage)

    # displays the region along the probe track
    # probe: string, the probe to be displayed from the drop down
    def updateDisplay(self, probe, restore=False):
        self.showMask = False
        self.showProbe = True
        #print('Probe', probe)
        self.generateLine(probe)

        p = probe.replace(' ', '_')
        self.volArray = np.array(Image.open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_slice.png'.format(p)))) # read slice 
        self.mask = Image.open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_mask.png'.format(p))) # read mask
        self.blended = np.array(Image.open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_overlay.png'.format(p)))) # read overlay
        with open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_labels.pickle'.format(p)), 'rb') as handle:
            self.labels_pos = pickle.load(handle)
        with open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_areas.pickle'.format(p)), 'rb') as handle:
            self.ccfAreas = pickle.load(handle)

        self.int_arrays = {}
        self.int_arrays['red'] = self.volArray[:, :, 0].copy()
        self.int_arrays['green'] = self.volArray[:, :, 1].copy()
        # modify based on if red and green channel toggle have been checked or not
        if self.isRedChecked:
            self.redOld = self.volArray[:, :, 0].copy()
            self.volArray[:, :, 0] = 0

        if self.isGreenChecked:
            self.greenOld = self.volArray[:, :, 1].copy()
            self.volArray[:, :, 1] = 0

        
        if not self.isRedChecked:
            self.volArray[:, :, 0] = self.imageClip('red', self.redSlider.value())
        if not self.isGreenChecked:
            self.volArray[:, :, 1] = self.imageClip('green', self.greenSlider.value())

        #self.volArray = np.array(vol_mask)
        #self.volArray *= 2
        rot = np.rot90(self.volArray)
        flip = np.flipud(rot)
        rot_mask = np.rot90(self.mask)
        flip_mask = np.flipud(rot_mask)

        self.imageMask.setImage(flip_mask, levels=(0, 255), autoRange=False)
        self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
        view = self.image.getView()

        self.points = [[80, t] for t in range(self.linepts.shape[0])]
        self.plItem = pg.ScatterPlotItem(pos=self.points, pen=QtGui.QPen(QColor('red')), brush=QtGui.QBrush(QColor('red')))
        #self.plItem.setClickable(True)
        #self.plItem.sigClicked.connect(self.onclickProbe)

        self.displayCCFAreas()

        if not restore:
            if self.path != '':
                view.addItem(self.plots['unit_density'].channelsPlot)
                view.addItem(self.plots['unit_density'].ogChannelsPlot)
                self.displayPlotImage()
                #view.addItem(self.plots[self.metrics.currentText().lower()].channelsPlot)

                #view.addItem(self.plItem)
        else: # read from dictionary storing saved plots for the probe
            #view.addItem(self.plItem)
            plot_items = self.alignments[probe]
           
            self.plots['unit_density'].channels = plot_items[0].copy()
            self.plots['unit_density'].channelsPlot.setData(pos=np.array(plot_items[0], dtype=float), adj=np.array(self.plots['unit_density'].adj, dtype=int))
            #self.plots[self.metrics.currentText().lower()].channelsPlot.setData(pos=np.array(plot_items[1], dtype=float), adj=np.array(self.plots['unit_density'].adj, dtype=int))
            #self.plots[self.metrics.currentText().lower()].channels = plot_items[1].copy()


            view.addItem(self.plots['unit_density'].channelsPlot) # add unit plot
            self.plots['unit_density'].updatePlotChannelPositions(self.plots['unit_density'].channels)
            view.addItem(self.plots['unit_density'].ogChannelsPlot)
            #view.addItem(self.plots[self.metrics.currentText().lower()].channelsPlot) # add metric plot
            self.displayPlotImage()
            y = [t[1] for t in self.plots['unit_density'].channels]
            self.anchorPts = plot_items[2].copy() # add line pts

            for pts in self.anchorPts:
                if pts[0][1] in y:
                    line_item = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=3)
                    line_item.sigClicked.connect(self.clickLine)
                    self.lineItems.append(line_item)

            for pts in self.anchorPts:
                if pts[0][1] not in y:
                    self.anchorPts.remove(pts)
        
            for item in self.lineItems:
                view.addItem(item)

            self.pointsAdded = plot_items[3].copy() # restore alignment
            self.displayImageCCFAreas(y)

            for point in self.pointsAdded:
                if point not in y:
                    self.pointsAdded.remove(point)

            for point in self.pointsAdded:
                if point in y:
                    ind = y.index(point)
                    text = pg.TextItem(str(ind))
                    text.setPos(80, point)
                    self.textItems.append(text)

                    plotItem = pg.ScatterPlotItem(pos=[[i, self.plots['unit_density'].ogChannelsShift[ind][1]] for i in range(-700, 
                                                                                                    int(self.plots['unit_density'].ogChannelsShift[ind][0]))], 
                                                    pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=3)
                    view.addItem(text)
                    view.addItem(plotItem)
                    self.plotImageItems.append(plotItem)
            
            print('Length of points added, anchor points', len(self.pointsAdded), len(self.anchorPts))
        #view.addItem(self.textItem)
        
if __name__ == '__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID

    app = QApplication(sys.argv)
    if args.templeton:
        if not args.record:
            v = VolumeAlignment(mouse_id, templeton=True,  base_path=args.templeton)
        elif args.oldDirStruct:
            if args.synology:
                v = VolumeAlignment(mouse_id, templeton=True, old_struct=True, synology=True, base_path=args.templeton, record_node=args.record)
            else:
                v = VolumeAlignment(mouse_id, templeton=True, old_struct=True, base_path=args.templeton, record_node=args.record)
        else:
            v = VolumeAlignment(mouse_id, templeton=True, base_path=args.templeton, record_node=args.record)
    elif args.dr:
        v = VolumeAlignment(mouse_id, dr=True)
    else:
        v = VolumeAlignment(mouse_id)
        
    sys.exit(app.exec_())