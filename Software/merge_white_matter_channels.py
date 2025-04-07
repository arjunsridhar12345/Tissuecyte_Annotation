import pandas as pd
import numpy.typing as npt
import numpy as np
import pathlib
import pickle
import matplotlib.pyplot as plt
from typing import Any
from scipy.spatial.distance import cdist
from clean_structure_acronym import clean_channel_annotations
import argparse
import warnings
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')

TISSUECYTE_PATH = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')

def get_alignment_paths(mouse_id: str) -> tuple[pathlib.Path, ...]:
    """
    >>> anchors = get_alignment_paths('626791')
    >>> len(anchors)
    11
    """
    return tuple(TISSUECYTE_PATH.glob(f'{mouse_id}/anchors/*'))

def open_alignment(anchor_file: pathlib.Path) -> list[Any]:
    """
    >>> anchor_path = get_alignment_paths('626791')[0]
    >>> anchor_data = open_alignment(anchor_path)
    >>> len(anchor_data)
    4
    """
    with open(anchor_file, 'rb') as f:
        return pickle.load(f)

def get_channel_file_paths(mouse_id: str) -> tuple[pathlib.Path, ...]:
    """
    >>> channel_file_paths = get_channel_file_paths('626791')
    >>> len(channel_file_paths)
    12
    """
    return tuple(TISSUECYTE_PATH.glob(f'{mouse_id}/*_channels*_warped.csv'))

def get_closest_non_white_matter_area(df_channels_by_anchor: pd.DataFrame, white_matter_coordinate: npt.NDArray[np.int64]) -> str:
    distance = cdist(df_channels_by_anchor[['AP', 'DV', 'ML']].to_numpy(dtype=np.int64), white_matter_coordinate, metric='euclidean')
    return df_channels_by_anchor['region'].to_numpy()[distance.argmin(axis=0)]

def assign_area(df_channels: pd.DataFrame, df_channels_between_anchors: pd.DataFrame) -> None:
    if pd.isna(df_channels_between_anchors['region']).all():
        return
    
    df_non_white_matter_between_anchors = df_channels_between_anchors[~pd.isna(df_channels_between_anchors['region'])]
    if len(df_non_white_matter_between_anchors) == len(df_channels_between_anchors): # no white matter
        return

    df_white_matter_only = df_channels_between_anchors[pd.isna(df_channels_between_anchors['region'])]
    for index, row in df_white_matter_only.iterrows():
        area = get_closest_non_white_matter_area(df_non_white_matter_between_anchors, np.array([row[['AP', 'DV', 'ML']].to_numpy()], 
                                                                                                dtype=np.int64))[0]
        df_channels.at[index, 'postprocessed_region'] = area

def execute_white_matter_assignment(df_channels: pd.DataFrame, channel_anchors: list[int]) -> None:
    df_channels['postprocessed_region'] = df_channels['region']
    df_channels['raw_location'] = df_channels['region']

    for i in range(len(channel_anchors) - 1):
        upper_channel = channel_anchors[i]
        lower_channel = channel_anchors[i + 1]

        df_channels_between_anchors = df_channels[(df_channels['channel'] > lower_channel) & (df_channels['channel'] < upper_channel)]  
        assign_area(df_channels, df_channels_between_anchors)
    
    last_anchor = channel_anchors[-1]
    # assumes if one anchor, use below since above is out of brain
    df_channels_below_last_anchor = df_channels[df_channels['channel'] < last_anchor] 
    assign_area(df_channels, df_channels_below_last_anchor)

def plot_postprocessed_regions(df_channels: pd.DataFrame, channel_anchors: list[int]) -> None:
    fig, ax = plt.subplots()
    
    for index, row in df_channels.iterrows():
        ax.text(0, index, row['region'])
        ax.text(5, index, row['postprocessed_region'], color='red')
    
    ax.set_ylim(top=384, bottom=0)
    ax.set_xlim(left=0, right=6)

    for channel in channel_anchors:
        ax.axhline(y=channel, xmin=0, xmax=6)
    fig.show()

def process_white_matter_channels(mouse_id: str) -> None:
    alignment_paths = sorted(get_alignment_paths(mouse_id))
    channel_file_paths = sorted(get_channel_file_paths(mouse_id))

    for i in range(len(alignment_paths)):
        alignment_path_str = alignment_paths[i].stem
        probe_day = alignment_path_str[0:alignment_path_str.index('_')+3]

        alignments = open_alignment(alignment_paths[i])
        anchors = alignments[3]
        y_position = [position[1] for position in alignments[0]]

        channel_anchors = [y_position.index(anchor) for anchor in anchors if anchor in y_position]
        channel_anchors = sorted(channel_anchors, reverse=True)
        channel_file_path = tuple(channel_file_path for channel_file_path in channel_file_paths if probe_day in str(channel_file_path))

        if not channel_file_path:
            warnings.warn(f'No channels csv for subject {mouse_id} and probe day {probe_day}', stacklevel=2)
            continue

        if len(channel_anchors) == 0:
            warnings.warn(f'No anchors aligned for subject {mouse_id} and probe day {probe_day}', stacklevel=2)
            continue
        
        df_channels = pd.read_csv(channel_file_path[0])

        execute_white_matter_assignment(df_channels, channel_anchors)
        clean_channel_annotations(mouse_id, channel_file_path[0], df_channels)
        #plot_postprocessed_regions(df_channels, channel_anchors)
        
def test_assignment() -> None:
    df_channels = pd.read_csv(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\366122\Probe_A1_channels_366122_warped.csv")
    channel_anchors = [236, 232, 224]
    execute_white_matter_assignment(df_channels, channel_anchors)
    plot_postprocessed_regions(df_channels, channel_anchors)
    print()

def run_tests() -> None:
    import doctest

    doctest.testmod(
        optionflags=(doctest.IGNORE_EXCEPTION_DETAIL | doctest.NORMALIZE_WHITESPACE)
    )

if __name__ == '__main__':
    #run_tests()
    dest_path = pathlib.Path(r"\\allen\programs\mindscope\workgroups\dynamicrouting\arjun")
    source_path = pathlib.Path(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte")

    args = parser.parse_args()
    mouse_id = args.mouseID

    process_white_matter_channels(mouse_id)
    shutil.copytree(source_path / mouse_id / 'images', dest_path / 'slice_images' / mouse_id, dirs_exist_ok=True)
    shutil.copytree(source_path / mouse_id / 'anchors', dest_path / 'alignment_anchors' / mouse_id, dirs_exist_ok=True)
    shutil.copytree(source_path / mouse_id / 'image_plots', dest_path / 'correlation_plots' / mouse_id, dirs_exist_ok=True)
