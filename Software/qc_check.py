import glob
from re import M
import numpy as np
import pathlib
import argparse
from PyQt5 import QtGui
from matplotlib import cm
import matplotlib.pyplot as plt
from generate_metrics_paths import generate_metrics_path_days
import pickle

# class for qc on kilosort output
class qcChecker():
    def __init__(self, kilo_sort_path:pathlib.Path, mouse_id:str, probe:str):
        self.kiloPath = kilo_sort_path
        self.mouseID = mouse_id
        self.probe = probe

        self.channelCoords = np.load(pathlib.Path(self.kiloPath, 'channel_positions.npy'), mmap_mode='r') # local coordinates file
        print(self.channelCoords.dtype)
        self.rawInd = np.load(pathlib.Path(self.kiloPath, 'channel_map.npy'), mmap_mode='r') # raw ind file
        print(self.rawInd.dtype)
        self.chnMin = np.min(self.channelCoords[:, 1])
        self.chnMax = np.max(self.channelCoords[:, 1])

        self.spikeAmps = np.load(pathlib.Path(self.kiloPath, 'amplitudes.npy'), mmap_mode='r')[:, 0] # spikes amps file
        print(self.spikeAmps.dtype)
        self.spikeClusters = np.load(pathlib.Path(self.kiloPath, 'spike_clusters.npy'), mmap_mode='r') # spike clusters file
        print(self.spikeClusters.dtype)
        self.spikeDepths = np.load(pathlib.Path(self.kiloPath, 'spike_depths.npy'),  mmap_mode='r') # spike depths
        print(self.spikeDepths.dtype)
        self.spikeTimes = np.load(pathlib.Path(self.kiloPath, 'spike_times.npy'), mmap_mode='r') / 30000. # spike times
        print(self.spikeTimes.dtype)
        self.spikeIdx = np.arange(self.spikeClusters.size)
        # Filter for nans in depths and also in amps
        self.kpIdx = np.where(~np.isnan(self.spikeDepths[self.spikeIdx]) &
                               ~np.isnan(self.spikeAmps[self.spikeIdx]))[0]
    

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

    # firing rate drift plot, uses bincount helper
    def get_fr(self):
        T_BIN = 0.05
        D_BIN = 5
        chn_min = np.min(np.r_[0, self.spikeDepths[self.spikeIdx][self.kpIdx]])
        chn_max = np.max(np.r_[3840, self.spikeDepths[self.spikeIdx][self.kpIdx]])
        n, times, depths = self.bincount2D(self.spikeTimes[self.spikeIdx][self.kpIdx],
                                          self.spikeDepths[self.spikeIdx][self.kpIdx],
                                          T_BIN, D_BIN, ylim=[chn_min, chn_max])

        img = n.T / T_BIN
        xscale = (times[-1] - times[0]) / img.shape[0]
        yscale = (depths[-1] - depths[0]) / img.shape[1]

        n = np.flipud(n)
        plt.imshow(n, cmap='binary', extent=[0, 6000, 0, 4000], vmin=0, vmax=1)
        plt.show()

    # gets the amplitude
    def get_amp(self):
        A_BIN = 10
        amp_range = np.quantile(self.spikeAmps[self.spikeIdx][self.kpIdx], [0, 0.9])
        amp_bins = np.linspace(amp_range[0], amp_range[1], A_BIN)
        colour_bin = np.linspace(0.0, 1.0, A_BIN + 1)
        colours = ((cm.get_cmap('BuPu')(colour_bin)[np.newaxis, :, :3][0]))
        spikes_colours = np.empty(self.spikeAmps[self.spikeIdx][self.kpIdx].size,
                                      dtype=object)
        spikes_size = np.empty(self.spikeAmps[self.spikeIdx][self.kpIdx].size)
        for iA in range(amp_bins.size):
            if iA == (amp_bins.size - 1):
                idx = np.where((self.spikeAmps[self.spikeIdx][self.kpIdx] >
                                    amp_bins[iA]))[0]
                # Make saturated spikes a very dark purple
                spikes_colours[idx] = [(64 / 255, 0 / 255, 128 / 255) for i in range(len(spikes_colours[idx]))]
            else:
                idx = np.where((self.spikeAmps[self.spikeIdx][self.kpIdx] >
                                    amp_bins[iA]) &
                                   (self.spikeAmps[self.spikeIdx][self.kpIdx] <=
                                    amp_bins[iA + 1]))[0]
                spikes_colours[idx] = [colours[iA] for i in range(len(spikes_colours[idx]))]

            spikes_size[idx] = iA / (A_BIN / 4)

        plt.scatter(self.spikeTimes[self.spikeIdx][self.kpIdx][0:-1:100], self.spikeDepths[self.spikeIdx][self.kpIdx][0:-1:100], s=spikes_size[0:-1:100],
                    c=spikes_colours[0:-1:100])
        plt.show()

    def get_correlation_data_img(self):
        T_BIN = 0.05
        D_BIN = 40

        chn_min = np.min(np.r_[0, self.spikeDepths[self.spikeIdx][self.kpIdx]])
        chn_max = np.max(np.r_[3840, self.spikeDepths[self.spikeIdx][self.kpIdx]])
        R, times, depths = self.bincount2D(self.spikeTimes[self.spikeIdx][self.kpIdx],
                                          self.spikeDepths[self.spikeIdx][self.kpIdx],
                                          T_BIN, D_BIN, ylim=[chn_min, chn_max])
        corr = np.corrcoef(R)
        corr[np.isnan(corr)] = 0
        scale = 384 / corr.shape[0]
   
        data_img = {
                'img': corr,
                'scale': np.array([scale, scale]),
                'levels': np.array([np.min(corr), np.max(corr)]),
                'offset': np.array([0, 0]),
                'xrange': np.array([0, 384]),
                'cmap': 'viridis',
                'title': 'Correlation',
                'xaxis': 'Distance from probe tip (um)'
        }

        path = pathlib.Path('/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/image_plots'.format(self.mouseID))

        if not path.exists():
            path.mkdir()

        with open(pathlib.Path(path, '{}_corr.pickle'.format(self.probe)), 'wb') as f:
            pickle.dump(data_img, f)

if __name__ == '__main__':
    probe = 'probeB'
    session = 'DRpilot_644866_20230210'
    path = r'\\allen\programs\mindscope\workgroups\dynamicrouting\datajoint\inbox\ks_paramset_idx_1\{}\{}_{}_sorted\continuous\Neuropix-PXI-100.0'.format(
                    session, session, probe)
    kilo_sort_path = pathlib.Path(path)
    qc = qcChecker(kilo_sort_path)
    qc.get_fr()