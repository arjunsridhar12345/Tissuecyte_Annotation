import pathlib
import pandas as pd
import pickle
import SimpleITK as sitk
import argparse
import numpy.typing as npt
import logging

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

def get_probe_trajectory_for_areas(areas_of_interest: list[str], area_trajectories:dict[str, npt.NDArray]) -> list[npt.NDArray]:
    all_probes_areas = []

    for area in areas_of_interest:
        if area not in area_trajectories:
            logging.debug(f'{area} has not been hit yet')
        else:
            for area_trajectories_per_area in area_trajectories[area]:
                area_trajectories_per_area_filtered = area_trajectories_per_area[area_trajectories_per_area[:, 1] > 0]
                all_probes_areas.append(area_trajectories_per_area_filtered)
    
    return all_probes_areas