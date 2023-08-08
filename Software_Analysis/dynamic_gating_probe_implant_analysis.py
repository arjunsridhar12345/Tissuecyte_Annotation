from importlib.resources import path
from itertools import groupby
import struct
import numpy as np
import SimpleITK as sitk
import visvis as vis

from PyQt5.QtCore import Qt, QAbstractTableModel
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, QTableView, \
    QVBoxLayout, QWidget, QPushButton, QGridLayout, QCheckBox, QFormLayout, QFileDialog, QComboBox, QMessageBox, QLineEdit
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5.QtGui import QColor

import pandas as pd
import pathlib
import sys
import itertools
from scipy.spatial.transform import Rotation as R

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from collections import Counter

import seaborn as sns
import pickle

import sqlite3
from sqlalchemy import  create_engine

backend = 'pyqt5'
app = vis.use(backend)

# ------------------------------------------------------------------------------------
# allow integers >8 bytes to be stored in sqlite3
sqlite3.register_adapter(np.int64, lambda val: int(val))
sqlite3.register_adapter(np.int32, lambda val: int(val))
# ------------------------------------------------------------------------------------
DB_PATH = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/dr_master.db')
# with contextlib.suppress(OSError):
#     DB_PATH.unlink()
sqlite3.connect(DB_PATH).close()
DB = f"sqlite:///{DB_PATH}"
ENGINE = create_engine(DB, echo=False)

class ProbeHoleImplantViewer(QWidget):
    def __init__(self):
        super().__init__()
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior')
        #self.probeHoleDataFrame = pd.read_excel(pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/dr_master_sheet.xlsx'))
        self.probeHoleDataFrame = pd.read_sql_table('channel_ccf_coords', con=ENGINE)
        self.sessionDataFrame = pd.read_sql_table('session_metadata', con=ENGINE)
        self.getProbeHoleImplantCombinations(self.probeHoleDataFrame)
        self.annotationImage = sitk.ReadImage(pathlib.Path(self.workingDirectory, 'tissuecyte', 'field_reference', 'average_template_25.nrrd'))

        with open(pathlib.Path(self.workingDirectory, 'tissuecyte', 'field_reference', 'acrnm_map.pkl'), 'rb') as f:
            self.acrnm_map = pickle.load(f)

        self.anno = sitk.GetArrayFromImage(sitk.ReadImage(pathlib.Path(self.workingDirectory, 'tissuecyte', 'field_reference', 'ccf_ano.mhd')))
        self.volume = sitk.GetArrayFromImage(self.annotationImage).T

        self.initializeUI()

    def initializeUI(self):
        self.mainLayout = QVBoxLayout()

        self.plotLayout = QHBoxLayout()
        self.mainFrame = QWidget()

        # Make figure using "self" as a parent
        self.Figure = app.GetFigureClass()
        self.fig = self.Figure(self)
        self.ax = vis.subplot(111)
        self.plotLayout.addWidget(self.fig._widget)

        self.createMainFrame(self.mainFrame)
        self.plotLayout.addWidget(self.mainFrame)

        self.mainLayout.addLayout(self.plotLayout)
        # probe drop downs
        self.probeDropDown = QComboBox()
        self.probeDropDown.addItem('Probes')
        self.probeDropDown.currentTextChanged.connect(self.probeDropDownChanged)
        for probe in self.uniqueProbes:
            self.probeDropDown.addItem(probe)

        self.holeDropDown = QComboBox()
        self.holeDropDown.addItem('Hole')
        #self.holeDropDown.currentTextChanged.connect(self.holeDropDownChanged)
        for hole in self.uniqueHoles:
            self.holeDropDown.addItem(str(hole))

        self.implantChanged = True
        self.implantDropDown = QComboBox()
        self.implantDropDown.addItem('Implant')
        self.implantDropDown.currentTextChanged.connect(self.implantDropDownChanged)
        for implant in self.uniqueImplants:
            self.implantDropDown.addItem(str(implant))

        self.probeHoleImplantLayout = QHBoxLayout()
        self.probeHoleImplantLayout.addWidget(self.probeDropDown)
        self.probeHoleImplantLayout.addWidget(self.implantDropDown)
        self.probeHoleImplantLayout.addWidget(self.holeDropDown)
        self.mainLayout.addLayout(self.probeHoleImplantLayout)

        self.viewProbeHoleImplantCCFButton = QPushButton('View probe implant hole regions')
        self.viewProbeHoleImplantCCFButton.clicked.connect(self.viewProbeHoleImplantRegions)
        self.probeHoleImplantLayout.addWidget(self.viewProbeHoleImplantCCFButton)

        self.updateDisplay()
        self.setLayout(self.mainLayout)
        self.showMaximized()

    def viewProbeHoleImplantRegions(self):
        plt.close()
        probe = self.probeDropDown.currentText()
        hole = self.holeDropDown.currentText()
        implant = self.implantDropDown.currentText()

        mouse_ids_probe_days = []

        for item in self.probeHoleImplantList:
            if hole == 'nan':
                if item[1] == probe and item[3] == implant:
                    # MID, probe, probe+day, session
                    mouse_ids_probe_days.append((item[0], item[1], item[1]+str(item[4]), item[5]))
            else:
                if item[1] == probe and item[2] == hole and str(item[3]) == implant:
                    # MID, probe, probe+day, session
                    mouse_ids_probe_days.append((item[0], item[1], item[1]+str(item[4]), item[5]))

        self.updateDisplay(mouse_ids_probe_days=mouse_ids_probe_days)

    def createMainFrame(self, frame:QWidget):
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        
        self.dpi = 100
        self.fig_2d = Figure((10, 10), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig_2d)
        self.canvas.setParent(frame)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig_2d.add_subplot(111)

        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, frame)

    # explicity close visvis widget
    def closeEvent(self, event):
        self.fig._widget.close()

    # updates the graph histogram
    def updateGraphPlot(self, regions_data:dict, colors:list):
        self.axes.clear()

        if len(regions_data['CCF_Region']) > 0:
            df = pd.DataFrame(regions_data)
            sns.barplot(x='CCF_Region', y='Channels_Per_CCF_Area', hue='Probe_Day_Mouse', orient='v', data=df, ax=self.axes, dodge=False,
                        palette=colors)

            plt.setp(self.axes.get_xticklabels(), rotation=30, horizontalalignment='center')
            self.axes.set_xticklabels(self.axes.get_xticklabels(), fontsize=6)
            self.axes.legend(title='Probe Day Mouse', loc='upper right')

        self.canvas.draw()
        plt.show()
        """
        self.axes.bar(list(region_counts.keys()), list(region_counts.values()), width=0.4)
        self.axes.set_xlabel('CCF Area Region')
        self.axes.set_ylabel('Total Count For Probe Hole Implant Combo')
        plt.setp(self.axes.get_xticklabels(), rotation=30, horizontalalignment='right')
        self.canvas.draw()
        plt.show()
        """
    # on change for probe drop down
    def probeDropDownChanged(self):
        if self.probeDropDown.currentText() != 'Probes':
            if not self.implantChanged:
                self.implantChanged = True

            self.holeDropDown.clear()
            self.implantDropDown.clear()

            self.df_probe_based = self.probeHoleDataFrame.loc[self.probeHoleDataFrame['Probe'] == self.probeDropDown.currentText()]
            unique_implants = self.df_probe_based['Implant'].unique()

            for implant in unique_implants:
                self.implantDropDown.addItem(str(implant))

            self.implantChanged = False
            self.implantDropDownChanged()

    def implantDropDownChanged(self):
        if not self.implantChanged:
            self.holeDropDown.clear()
            if self.implantDropDown != 'Implant':
                self.df_probe_implant_based = self.df_probe_based.loc[self.df_probe_based['Implant'].astype(str) == self.implantDropDown.currentText()]

                unique_holes = self.df_probe_implant_based['Hole'].unique()
                for hole in unique_holes:
                    self.holeDropDown.addItem(str(hole))

            #self.holeChanged = True
    
    # helper to get the probe locator info
    def getProbeLocatorsInfo(self, mouse_ids_probe_days:list, scatter_colors:dict) -> dict:
        probe_locators_dict = {'Probe_Day_Mouse': [], 'Probe_Loc_x': [], 'Probe_Loc_y': [], 'colors':[]}
        np_exp_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')

        for item in mouse_ids_probe_days:
            mouse_id = item[0]
            session = str(item[1])
            probe_day = item[2]
            #probe_loc_x = item[2]
            #probe_loc_y = item[3]

            df_path = pathlib.Path(pathlib.Path(self.workingDirectory, 'tissuecyte', str(mouse_id), 'Probe_{}_channels_{}_warped.csv'.format(probe_day, str(mouse_id))))
            if df_path.exists():
                probe_locator_csv_path = pathlib.Path(np_exp_path, session, 'probelocator_{}_insertions_in_rig_image_space.csv'.format(session))

                if probe_locator_csv_path.exists():
                    df_probe_locator = pd.read_csv(probe_locator_csv_path)
                    locator_values = df_probe_locator[probe_day[0]].values

                    probe_day_label = 'Probe{} Day{} Mouse{}'.format(probe_day[0], probe_day[1], mouse_id)
                    probe_locators_dict['Probe_Day_Mouse'].append(probe_day_label)
                    probe_locators_dict['Probe_Loc_x'].append(locator_values[0])
                    probe_locators_dict['Probe_Loc_y'].append(locator_values[1])
                    probe_locators_dict['colors'].append(scatter_colors[probe_day_label])

        return probe_locators_dict

    def plot2DProbeLocators(self, mouse_ids_probe_days:list, scatter_colors:dict):
        probe_locators_dict = self.getProbeLocatorsInfo(mouse_ids_probe_days, scatter_colors)

        if len(probe_locators_dict['Probe_Loc_x']) > 0:
            temp = probe_locators_dict['Probe_Loc_x'].copy()
            probe_locators_dict['Probe_Loc_x'] = probe_locators_dict['Probe_Loc_y']
            probe_locators_dict['Probe_Loc_y'] = temp

            probe_locators_df = pd.DataFrame(probe_locators_dict)

            sns.scatterplot(x='Probe_Loc_x', y='Probe_Loc_y', hue='Probe_Day_Mouse', data=probe_locators_df, palette=probe_locators_dict['colors'],
                            legend=False)
            plt.show()
    
    # gets the coordinates of the given area
    def getBrainAreaCoordinates(self, area:str) -> np.ndarray:
        structure_id = self.acrnm_map[area]
        brain_area_coordinates = np.where(self.anno == structure_id)

        return brain_area_coordinates

    # updates the 3d display
    def updateDisplay(self, mouse_ids_probe_days: list=None):
        if mouse_ids_probe_days is not None:
            self.ax.bgcolor='k'
            self.ax.axis.axisColor = 'w'
            vis.cla()

            vol = np.rot90(self.volume, k=2, axes=(0,2))
            vis.volshow3(vol)

            temp = self.volume.shape[2] 
            legend = []
            colors = []

            area = 'MD'
            brain_coordinates = self.getBrainAreaCoordinates(area)
            center_of_mass = np.average(brain_coordinates, axis=1)

            brain_coordinates_color = (1, 0, 0)

            mouse_id_probe_days_used = []

            r = R.from_rotvec([0, np.pi, 0])
            regions_data = {'CCF_Region': [], 'Channels_Per_CCF_Area': [], 'Probe_Day_Mouse': []}
            last_mouse_id = -1
            last_color = None

            cm = plt.get_cmap('jet')
            clrs = [cm(1.*i/len(mouse_ids_probe_days)) for i in range(len(mouse_ids_probe_days))]
            scatter_colors = {}
            i = 0

            coords = None

            channel_area_distance = {'Probe': [], 'AP_Distance': [], 'DV_Distance': [], 'ML_Distance': [], 'Area': []}

            for item in mouse_ids_probe_days:
                mouse_id = item[0]
                probe_day = item[2]
                session = item[3]
                session = session[0:session.index('_')]

                df_path = pathlib.Path(pathlib.Path(self.workingDirectory, 'tissuecyte', str(mouse_id), 'Probe_{}_channels_{}_warped.csv'.format(probe_day, str(mouse_id))))
                if df_path.exists():
                    df_channels = pd.read_csv(pathlib.Path(self.workingDirectory, 'tissuecyte', str(mouse_id), 'Probe_{}_channels_{}_warped.csv'.format(probe_day, str(mouse_id))))
                    df_channels = df_channels.loc[df_channels['AP'] >= 0]
                    df_channels = df_channels.loc[~pd.isna(df_channels['region'])]
                    coords = df_channels[['AP', 'DV', 'ML']].to_numpy().astype(np.float64)

                    probe_day_label = 'Probe{} Day{} Mouse{}'.format(probe_day[0], probe_day[1], mouse_id)
                    region_counts = Counter(df_channels['region'])

                    for key in region_counts.keys():
                        regions_data['CCF_Region'].append(key)
                        regions_data['Channels_Per_CCF_Area'].append(region_counts[key])
                        regions_data['Probe_Day_Mouse'].append(probe_day_label)

                    coords[:, 0] -= self.volume.shape[0] / 2
                    coords[:, 1] -= self.volume.shape[1] / 2
                    coords[:, 2] -= self.volume.shape[2] / 2
                    coords = r.apply(coords)

                    coords[:, 0] += self.volume.shape[0] / 2
                    coords[:, 1] += self.volume.shape[1] / 2
                    coords[:, 2] += self.volume.shape[2] / 2

                    min_ap_distance = np.linalg.norm(coords[:, 0] - center_of_mass[0]).min()
                    min_dv_distance = np.linalg.norm(coords[:, 1] - center_of_mass[1]).min()
                    min_ml_distance = np.linalg.norm(coords[:, 2] - center_of_mass[2]).min()
                    channel_area_distance['Probe'].append('{}{}'.format(mouse_id,probe_day))
                    channel_area_distance['AP_Distance'].append(min_ap_distance)
                    channel_area_distance['DV_Distance'].append(min_dv_distance)
                    channel_area_distance['ML_Distance'].append(min_ml_distance)
                    channel_area_distance['Area'].append(area)

                    color = clrs[i][0:3]
                    colors.append(color)
                    i += 1
                    scatter_colors[probe_day_label] = color

                    vis.plot(temp - coords[:, 2], coords[:, 1], coords[:, 0], lw=0, mw=10, ms='x', mc=color, axes=self.ax)
                    legend.append('{} Mouse{} Probe{} Day{}'.format(session, mouse_id, probe_day[0], probe_day[1]))
            
            #vis.plot(temp - brain_coordinates[2], brain_coordinates[1], brain_coordinates[0], lw=0, mc=brain_coordinates_color, mw=10, ms='*', axes=self.ax)
            #legend.append('MD')
            self.ax.legend = legend
            self.ax.axis.xLabel = 'ML'
            self.ax.axis.yLabel = 'DV'
            self.ax.axis.zLabel = 'AP'
            
            #regions_unpacked = list(itertools.chain.from_iterable(regions_hit))
            self.updateGraphPlot(regions_data, colors)
            #self.plotClosestDistance(channel_area_distance)
            #self.plot2DProbeLocators(mouse_ids_probe_days, scatter_colors)
        else:
            vol = np.rot90(self.volume, k=2, axes=(0,2))
            vis.volshow3(vol)

    # plots the closest distance in each direction to the area
    def plotClosestDistance(self, channel_area_distance: dict):
        plt.close()
        df_distance = pd.DataFrame(channel_area_distance)

        fig, ax = plt.subplots(1, 3)

        sns.barplot(x='Area', y='AP_Distance', hue='Probe', orient='v', data=df_distance, ax=ax[0], errorbar=None, dodge=False)
        sns.barplot(x='Area', y='DV_Distance', hue='Probe', orient='v', data=df_distance, ax=ax[1], errorbar=None, dodge=False)
        sns.barplot(x='Area', y='ML_Distance', hue='Probe', orient='v', data=df_distance, ax=ax[2], errorbar=None, dodge=False)

        plt.show()

    # gets the combination of probe, hole, implant and mouse ids
    def getProbeHoleImplantCombinations(self, probe_hole_implant:pd.DataFrame):
        #self.probeHoleDataFrame = self.probeHoleDataFrame.loc[self.probeHoleDataFrame['annotated'] == True]
        self.probeHoleDataFrame = self.probeHoleDataFrame.merge(self.sessionDataFrame[['session', 'MID', 'day']], 
                                                                left_on=['MID', 'Day'], right_on=['MID', 'day'])
        self.probeHoleImplantList = self.probeHoleDataFrame[['MID', 'Probe', 'Hole', 'Implant', 'Day', 'session']].to_numpy().tolist()
        df_unique_probe_hole_implant = probe_hole_implant.apply(lambda col: col.unique())
        self.uniqueMouseIDs = df_unique_probe_hole_implant['MID']
        self.uniqueProbes = df_unique_probe_hole_implant['Probe']
        self.uniqueHoles = df_unique_probe_hole_implant['Hole']
        self.uniqueImplants = df_unique_probe_hole_implant['Implant']

        """
        combination_list = [self.uniqueMouseIDs, self.uniqueProbes, self.uniqueHoles, self.uniqueImplants]
        self.probeHoleImplantCombinations = list(itertools.product(*combination_list))
        self.dfProbeHoleImpantCombinations = pd.DataFrame(self.probeHoleImplantCombinations, columns=['MID', 'probe', 'hole', 'implant'])
        """

if __name__ == '__main__':
    app.Create()
    analysis_viewer = ProbeHoleImplantViewer()
    app.Run()