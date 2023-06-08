import numpy as np
import pathlib
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
from qc_check import qcChecker
from generate_metrics_paths import generate_metrics_path_ephys, generate_metrics_path_days
import contextlib
import numba

# gets the spike depth for a given kilo sort directory and saves numpy file in same directory
def generate_spike_depth(kilo_sort_dir:pathlib.Path) -> None:
    print(kilo_sort_dir)
    output_path = pathlib.Path(kilo_sort_dir, 'spike_depths.npy') # output path for spike depths

    if not output_path.exists():
        sparse_features = np.load(pathlib.Path(kilo_sort_dir, 'pc_features.npy'), mmap_mode='r').squeeze().transpose((0, 2, 1)) # pc features
        print(sparse_features.shape)
        sparse_features_ind = np.load(pathlib.Path(kilo_sort_dir, 'pc_feature_ind.npy'), mmap_mode='r') # pc features ind
        spike_templates = np.load(pathlib.Path(kilo_sort_dir, 'spike_templates.npy'), mmap_mode='r')[:, 0]
        print(spike_templates.shape)
        spike_times = np.load(pathlib.Path(kilo_sort_dir, 'spike_times.npy'), mmap_mode='r')[:, 0]
        spike_times = spike_times / 30000.
        print(spike_times.shape)
        channel_positions = np.load(pathlib.Path(kilo_sort_dir, 'channel_positions.npy'), mmap_mode='r')

        nbatch = 50000
        c = 0
        spikes_depths = np.zeros_like(spike_times)
        nspi = spikes_depths.shape[0]

        while True:
            ispi = np.arange(c, min(c + nbatch, nspi))
            # take only first component
            features = sparse_features[ispi, :, 0]
            features = np.maximum(features, 0) ** 2  # takes only positive values into account


            ichannels = sparse_features_ind[spike_templates[ispi]].astype(np.uint32)
            # features = np.square(self.sparse_features.data[ispi, :, 0])
            # ichannels = self.sparse_features.cols[self.spike_templates[ispi]].astype(np.int64)
            ypos = channel_positions[ichannels, 1]
            #ypos = ypos[:, 0, :]
            with np.errstate(divide='ignore'):
                print('Features', features.shape)
                print('Ypos', ypos.shape)
                spikes_depths[ispi] = (np.sum(np.transpose(ypos * features) /np.sum(features, axis=1), axis=0))
            c += nbatch
            if c >= nspi:
                break

        np.save(output_path, spikes_depths)


if __name__ == '__main__':
    mouse_ids = ['614608', '615047', '615563', '615564', '623319', '623322', '623784', '623786', '626279', '632295', '632296', 
                 '633232', '636740', '637483', '638387', '640887', '640890', '642504']
    # TODO: 615048, 636740
    for mouse in ['636740']:
        metrics_paths = generate_metrics_path_days(pathlib.Path('/allen/programs/mindscope/workgroups/np-exp'), mouse)
        probe_letters = ['A', 'B', 'C', 'D', 'E', 'F']

        for day in metrics_paths:
            metric_paths_day = metrics_paths[day]

            for metric_path in metric_paths_day:
                kilo_sort_path = pathlib.Path(metric_path).parent
                letter = [probe_letter for probe_letter in probe_letters if probe_letter in str(kilo_sort_path)][0]
                probe = letter + str(day)
                    
                if not pathlib.Path('/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/image_plots/{}_corr.pickle'.format(mouse, probe)).exists():
                    generate_spike_depth(kilo_sort_path)
                    qcChecker(kilo_sort_path, mouse, probe).get_correlation_data_img()
                    
