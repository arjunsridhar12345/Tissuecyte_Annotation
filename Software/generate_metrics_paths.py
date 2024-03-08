# generates a dictionary with the paths for the metrics file

import pandas as pd
import argparse
import os
import pathlib
import glob
from typing import Union
import npc_lims


parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)

# gets the directories for relevatn mouse
def get_metrics_directory(base_path: Union[str, pathlib.Path],  mouse_id: str):
    directories = os.listdir(base_path)
    probe_directories = []

    for d in directories:
        if mouse_id in d and 'json' not in d and 'Shortcut' not in d:
            probe_directories.append(d)
    
    return probe_directories

def generate_templeton_metric_path_days(mouse_id: str, base_path: str, record_node: str, old_struct: bool=False, synology=False):
    record_node = 'Record Node {}'.format(record_node)
    base_path = pathlib.Path(base_path)
    mouse_dirs = get_metrics_directory(base_path, mouse_id)
    
    probe_metrics_dirs = {}
    metrics_path_days = {}

    for directory in mouse_dirs:
        if not old_struct:
            date = directory[directory.rindex('_')+1:]
            probe_sorted_dir = [d for d in os.listdir(os.path.join(base_path, directory)) if os.path.isdir(os.path.join(base_path, directory, d))][0]

            probe_dirs = [d for d in os.listdir(os.path.join(base_path, directory, probe_sorted_dir, record_node, 'experiment1', 'recording1', 'continuous')) 
                      if os.path.isdir(os.path.join(base_path, directory, probe_sorted_dir, record_node, 'experiment1', 'recording1', 'continuous'))]

            probe_metrics_dirs [date] = [os.path.join(base_path, directory, probe_sorted_dir, record_node, 'experiment1', 'recording1', 'continuous', d) 
                                     for d in probe_dirs]
        else:
            if not synology:
                date = directory[0:directory.index('_')]
                probe_dirs = [d for d in os.listdir(os.path.join(base_path, directory, record_node, 'experiment1', 'recording1', 'continuous')) 
                          if os.path.isdir(os.path.join(base_path, directory, record_node, 'experiment1', 'recording1', 'continuous'))]
                probe_metrics_dirs [date] = [os.path.join(base_path, directory, record_node, 'experiment1', 'recording1', 'continuous', d) for d in probe_dirs]
            else:
                date = directory[0:directory.index('_')]
                probe_dirs = [d for d in os.listdir(os.path.join(base_path, directory, directory, record_node, 'experiment1', 'recording1', 'continuous')) 
                          if os.path.isdir(os.path.join(base_path, directory, directory, record_node, 'experiment1', 'recording1', 'continuous'))]
                probe_metrics_dirs [date] = [os.path.join(base_path, directory, directory, record_node, 'experiment1', 'recording1', 'continuous', d) for d in probe_dirs]
            
    for date in probe_metrics_dirs:
        for d in probe_metrics_dirs[date]:
            files = [os.path.join(d, f) for f in os.listdir(d)]
            for f in files:
                if 'waveform_metrics.csv' in f:
                    waveform_metrics_path = f
                    if date not in metrics_path_days:
                        metrics_path_days[date] = [waveform_metrics_path]
                    else:
                        metrics_path_days[date].append(waveform_metrics_path)

    #print(metrics_path_days)
    return metrics_path_days

# gets the path for the metrics csv
def generate_metrics_path_days_codeocean(base_path: pathlib.Path, mouse_id: str) -> dict[int, str]:
    """
    >>> np_exp_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    >>> metrics_paths = generate_metrics_path_days_codeocean(np_exp_path, '681532')
    >>> metrics_paths[1][0]
    '//allen/programs/mindscope/workgroups/np-exp/681532_2023-10-16_0/probeA/continuous/Neuropix-PXI-100.0/metrics.csv'
    """
    sessions_mouse: list[npc_lims.SessionInfo] = []
    for session in npc_lims.get_session_info():
        if session.is_ephys and session.subject == mouse_id:
            sessions_mouse.append(session)

    dates = sorted([session.date for session in sessions_mouse])
    mouse_dirs = get_metrics_directory(base_path, mouse_id)
    days = [(dates.index(date) + 1, f'{mouse_id}_{date}_0') for date in dates if f'{mouse_id}_{date}_0' in mouse_dirs]
    metrics_path_days = {}

    for pair in days:
        session = pair[1]
        day = pair[0]
        metrics_files = list(base_path.glob(f'{session}/*/*/*/metrics.csv'))
        metrics_files = [metric_file.as_posix() for metric_file in metrics_files]
        metrics_path_days[day] = metrics_files
 
    return metrics_path_days

# gets the path for the metrics csv
def generate_metrics_path_days(base_path, mouse_id):
    mouse_dirs = get_metrics_directory(base_path, mouse_id)
    probe_metrics_dirs = []

    for directory in mouse_dirs:
        probe_dirs = [d for d in os.listdir(os.path.join(base_path, directory)) if os.path.isdir(os.path.join(base_path, directory, d))]

        probe_metrics_dirs.append([os.path.join(base_path, directory, d, 'continuous') for d in probe_dirs if os.path.exists(os.path.join(base_path, directory, d, 'continuous'))])
    
    metrics_path_days = {}
    days = [1, 2, 3, 4]
    i = 0
    for directory in sorted(probe_metrics_dirs):
        if len(directory) > 0:
            for d in directory:
                #date = d[d.index('_')+1:]
                #date = date[date.index('_')+1:date.index('\\')]
                date = days[i]
                files = [os.path.join(d, f) for f in os.listdir(d)]
                for f in files:
                    metrics = os.listdir(os.path.join(d, f))

                    waveform_metrics_path = [os.path.join(d, f, m) for m in metrics if m == 'metrics.csv']
   
                    if len(waveform_metrics_path) > 0:
                        waveform_metrics_path = waveform_metrics_path[0]
                        if date not in metrics_path_days:
                            metrics_path_days[date] = [waveform_metrics_path]
                        else:
                            metrics_path_days[date].append(waveform_metrics_path)
                
        i += 1

    #print(metrics_path_days)
    return metrics_path_days

def generate_metrics_path_ephys(base_path: pathlib.Path, mouse_id: str):
    dj_path = '//allen/programs/mindscope/workgroups/dynamicrouting/datajoint/inbox/ks_paramset_idx_1'
    mouse_dirs = get_metrics_directory(base_path, mouse_id)

    metrics_path_days = {}

    for directory in mouse_dirs:
        probe_dirs = [d for d in os.listdir(os.path.join(base_path, directory)) if os.path.isdir(os.path.join(base_path, directory, d))]

        if len(probe_dirs) > 0:
            behavior_dict = pd.read_pickle(os.path.join(base_path, directory, '{}.behavior.pkl'.format(directory)))
            ephys_day = behavior_dict['items']['behavior']['params']['stage']
            day = int(ephys_day[ephys_day.index('_') + 1:])

            metrics_dj = glob.glob('{}'.format(str(pathlib.Path(str(dj_path), directory, '*', '*', '*', 'metrics.csv'))))

            if len(metrics_dj) > 0:
                metrics_path_days[day] = metrics_dj
            else:
                metrics = glob.glob('{}'.format(str(pathlib.Path(str(base_path), directory, '*', '*', '*', 'metrics.csv'))))
                metrics_path_days[day] = metrics

    #print(metrics_path_days)
    return metrics_path_days

if __name__ == '__main__':
    import doctest

    doctest.testmod(
        optionflags=(doctest.IGNORE_EXCEPTION_DETAIL | doctest.NORMALIZE_WHITESPACE)
    )
    """
    base_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    output_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior')

    #args = parser.parse_args()
    #mouse_id = args.mouseID
    mouse_id = '681532'
    #generate_metrics_path_days(base_path, output_path, mouse_id)
    generate_metrics_path_days_codeocean(base_path, mouse_id)
    """