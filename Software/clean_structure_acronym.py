import numpy as np
import pathlib
import pandas as pd
import pickle
import SimpleITK as sitk
import argparse

STRUCTURE_TREE = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/ccf_structure_tree_2017.csv'))
with open(pathlib.Path(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\field_reference\acrnm_map.pkl"), 'rb') as f:
    ACRONYM_MAP = pickle.load(f)

ANNOTATION_VOLUME = sitk.GetArrayFromImage(sitk.ReadImage(pathlib.Path(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\field_reference\ccf_ano.mhd")))

ANNOTATION_PATH = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')

def strip_subregions_list(areas: list):
    areas_abbreviated = []

    for area in areas:
        areas_abbreviated.append(strip_subregions_layers(area))

    return areas_abbreviated

def strip_subregions_layers(areastr):
    if not isinstance(areastr, str) or areastr == 'root' or areastr == 'No Area' or areastr == 'Track not annotated':
        return areastr
    
    areastr = areastr.split('-')[0]
    #remove layer stuff
    areaname = STRUCTURE_TREE[STRUCTURE_TREE['acronym']==areastr]['name'].values
    if len(areaname)==0:
        return areastr
    else:
        areaname = areaname[0]
        
    if 'layer' in areaname:
        layer = areaname.split('layer')[-1].split(' ')[-1]
        areastr = areastr.replace(layer, '')
        
    #hack for ACA and MOp which doesn't play nice
    if 'ACAd' in areastr:
        areastr = 'ACAd'
        
    if 'ACAv' in areastr:
        areastr = 'ACAv'

    if 'MOp' in areastr:
        areastr = 'MOp'
        
    return areastr

def get_structure_acronym(point:tuple[int, int, int]) -> str:
    if point[1] < 0:
        return 'out of brain'
    
    structure_ids = tuple(ACRONYM_MAP.values())
    labels = tuple(ACRONYM_MAP.keys())

    structure_id = ANNOTATION_VOLUME[point[0], point[1], point[2]]
    if structure_id in structure_ids:
        index = structure_ids.index(structure_id)
        label = labels[index]
    else:
        label = 'root'
    
    return label

def clean_channel_annotations(mouse_id: str) -> None:
    annotation_path = pathlib.Path(f'//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{mouse_id}')
    channels_paths = tuple(annotation_path.glob(f'*_channels_{mouse_id}_warped.csv'))

    for channel_path in channels_paths:
        df_channels = pd.read_csv(channel_path)
        for index, row in df_channels.iterrows():
            if pd.isna(row.region) or row.region == 'out of brain':
                label = get_structure_acronym((row.AP, row.DV, row.ML))
                df_channels.loc[index, 'region'] = label
        
        df_channels['region_stripped'] = strip_subregions_list(df_channels['region'].tolist())

        output_path = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/arjun')
        output_dir = output_path / mouse_id
        if not output_dir.exists():
            output_dir.mkdir()
        
        file_name = channel_path.stem + '_processed.csv'
        df_channels.to_csv(output_dir / file_name, index=False)
    
if __name__ == '__main__':
    args = parser.parse_args()
    #mouse_id = args.mouseID
    mouse_ids = ['660023']
    for mouse_id in mouse_ids:
        clean_channel_annotations(mouse_id)