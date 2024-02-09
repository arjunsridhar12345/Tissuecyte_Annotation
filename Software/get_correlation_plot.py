import glob
from re import M
import numpy as np
import pathlib
import argparse
from PyQt5 import QtGui
from matplotlib import cm
import matplotlib.pyplot as plt
from generate_metrics_paths import generate_metrics_path_days, generate_metrics_path_days_codeocean
import pickle
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')

# class for qc on kilosort output
class qcChecker():
    def __init__(self, kilo_sort_path:pathlib.Path, mouse_id:str, probe:str, scale:int=384, channel_max:int=3840):
        self.kiloPath = kilo_sort_path
        self.mouseID = mouse_id
        self.probe = probe
        self.scale = scale
        self.channel_max = channel_max

        self.spikeDepths = np.load(pathlib.Path(self.kiloPath, 'spike_depths.npy'),  mmap_mode='r') # spike depths
        self.spikeTimes = np.load(pathlib.Path(self.kiloPath, 'spike_times.npy'), mmap_mode='r')
        self.spikeIdx = np.arange(self.spikeTimes.size)
        # Filter for nans in depths and also in amps
        self.kpIdx = np.where(~np.isnan(self.spikeDepths[self.spikeIdx]))[0]
    
    # helper function for firing rate 
    def bincount2D(self, x, y, xbin=0, ybin=0, xlim=None, ylim=None, weights=None):
        """
        Computes a 2D histogram by aggregating values in a 2D array.
        :param x: values to bin along the 2nd dimension (c-contiguous)
        :param y: values to bin along the 1st dimension
        :param xbin:
            scalar: bin size along 2nd dimension
            0: aggregate according to unique values
            array: aggregate according to exact values (count reduce operation)
        :param ybin:
            scalar: bin size along 1st dimension
            0: aggregate according to unique values
            array: aggregate according to exact values (count reduce operation)
        :param xlim: (optional) 2 values (array or list) that restrict range along 2nd dimension
        :param ylim: (optional) 2 values (array or list) that restrict range along 1st dimension
        :param weights: (optional) defaults to None, weights to apply to each value for aggregation
        :return: 3 numpy arrays MAP [ny,nx] image, xscale [nx], yscale [ny]
        """
        # if no bounds provided, use min/max of vectors
        if xlim is None:
            xlim = [np.min(x), np.max(x)]
        if ylim is None:
            ylim = [np.min(y), np.max(y)]

        def _get_scale_and_indices(v, bin, lim):
            # if bin is a nonzero scalar, this is a bin size: create scale and indices
            if np.isscalar(bin) and bin != 0:
                scale = np.arange(lim[0], lim[1] + bin / 2, bin)
                ind = (np.floor((v - lim[0]) / bin)).astype(np.int64)
            # if bin == 0, aggregate over unique values
            else:
                scale, ind = np.unique(v, return_inverse=True)
            return scale, ind

        xscale, xind = _get_scale_and_indices(x, xbin, xlim)
        yscale, yind = _get_scale_and_indices(y, ybin, ylim)
        # aggregate by using bincount on absolute indices for a 2d array
        nx, ny = [xscale.size, yscale.size]
        ind2d = np.ravel_multi_index(np.c_[yind, xind].transpose(), dims=(ny, nx))
        r = np.bincount(ind2d, minlength=nx * ny, weights=weights).reshape(ny, nx)

        # if a set of specific values is requested output an array matching the scale dimensions
        if not np.isscalar(xbin) and xbin.size > 1:
            _, iout, ir = np.intersect1d(xbin, xscale, return_indices=True)
            _r = r.copy()
            r = np.zeros((ny, xbin.size))
            r[:, iout] = _r[:, ir]
            xscale = xbin

        if not np.isscalar(ybin) and ybin.size > 1:
            _, iout, ir = np.intersect1d(ybin, yscale, return_indices=True)
            _r = r.copy()
            r = np.zeros((ybin.size, r.shape[1]))
            r[iout, :] = _r[ir, :]
            yscale = ybin

        return r, xscale, yscale

    def get_correlation_data_img(self):
        T_BIN = 0.05
        D_BIN = 40

        chn_min = np.min(np.r_[0, self.spikeDepths[self.spikeIdx][self.kpIdx]])
        chn_max = np.max(np.r_[self.channel_max, self.spikeDepths[self.spikeIdx][self.kpIdx]])
        R, times, depths = self.bincount2D(self.spikeTimes[self.spikeIdx][self.kpIdx],
                                          self.spikeDepths[self.spikeIdx][self.kpIdx],
                                          T_BIN, D_BIN, ylim=[chn_min, chn_max])
        corr = np.corrcoef(R)
        corr[corr < 0] = 0
        corr[np.isnan(corr)] = 0
        scale = self.scale / corr.shape[0]
   
        data_img = {
                'img': corr,
                'scale': np.array([scale, scale]),
                'levels': np.array([np.min(corr), np.max(corr)]),
                'offset': np.array([0, 0]),
                'xrange': np.array([0, self.channel_max / 10]),
                'cmap': 'viridis',
                'title': 'Correlation',
                'xaxis': 'Distance from probe tip (um)'
        }

        path = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/image_plots'.format(self.mouseID))

        if not path.exists():
            path.mkdir()

        with open(pathlib.Path(path, '{}_corr.pickle'.format(self.probe)), 'wb') as f:
            pickle.dump(data_img, f)

def get_correlation_data(mouse_id:str):
    metrics_paths = generate_metrics_path_days_codeocean(pathlib.Path('//allen/programs/mindscope/workgroups/np-exp'), mouse_id)
    probe_letters = ['A', 'B', 'C', 'D', 'E', 'F']

    for day in metrics_paths:
        metric_paths_day = metrics_paths[day]

        for metric_path in metric_paths_day:
            peak_channels = pd.read_csv(metric_path)['peak_channel']
            kilo_sort_path = pathlib.Path(metric_path).parent
            letter = [probe_letter for probe_letter in probe_letters if 'probe{}'.format(probe_letter) in str(kilo_sort_path)][0]
            probe = letter + str(day)
                    
            if not pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/image_plots/{}_corr.pickle'.format(mouse_id, probe)).exists():
                if peak_channels.max() > 383:
                    qcChecker(kilo_sort_path, mouse_id, probe, scale=768, channel_max=7680).get_correlation_data_img()
                else:
                    qcChecker(kilo_sort_path, mouse_id, probe).get_correlation_data_img()

if __name__ == '__main__':
    args = parser.parse_args()
    #mouse = args.mouseID
    mouse = '686740'
    get_correlation_data(mouse)
