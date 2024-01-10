from abc import get_cache_token
import numpy as np
import pandas as pd
import nrrd
from scipy import ndimage
import pathlib
import json

structure_tree = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/ccf_structure_tree_2017.csv'))
annotationData = nrrd.read('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/annotation_25.nrrd')[0]

def get_area_volume(area_acro, areavol = np.zeros_like(annotationData)):
    #areavol = np.zeros_like(annotationData)
    area_acro = np.copy(area_acro)
    areavol = np.copy(areavol)
    areaid = structure_tree[structure_tree['acronym']==area_acro]['id'].values[0]
    print(areaid)
    if np.sum(annotationData==areaid)==0:
        #find children
        children = structure_tree[structure_tree['parent_structure_id']==areaid]['acronym'].values
        for child in children:
            areavol = get_area_volume(child, areavol=areavol)
    
    else:
        areavol[annotationData==areaid] = 1
    
    areavol[:, :, int(456/2):] = 0 #only left hemi
    return areavol


def get_area_center_of_mass(area_acronym):
    areavolume = get_area_volume(area_acronym)
    com = np.array(ndimage.center_of_mass(areavolume))
    return (com*25).tolist()

if __name__ == '__main__':
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/all_regions_check.json'), 'r') as f:
        areas_check = json.load(f)

    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'r') as f:
        area_com_dict = json.load(f)

    for area in areas_check:
        if area not in area_com_dict:
            area_com_dict[area] = get_area_center_of_mass(area)

    #area_com_dict = {a:get_area_center_of_mass(a) for a in areas_check}
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'w') as f:
        json.dump(area_com_dict, f)
