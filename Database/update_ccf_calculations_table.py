import pandas as pd
import sqlite3
from sqlalchemy import create_engine
import pathlib
import numpy as np
import logging
import json

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

structure_tree = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/ccf_structure_tree_2017.csv'))

def strip_subregions_list(areas: list):
    areas_abbreviated = []

    for area in areas:
        areas_abbreviated.append(strip_subregions_layers(area))

    return areas_abbreviated

def strip_subregions_layers(areastr):
    if not isinstance(areastr, str):
        return areastr
    
    areastr = areastr.split('-')[0]
    #remove layer stuff
    areaname = structure_tree[structure_tree['acronym']==areastr]['name'].values
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

def count_channels_in_region(group, region):
     if region in group:
         return 1
        
     return 0

def min_distance_to_ccf_point(group:np.array, point:list):
    channel_coords = group*25
    distance = [np.linalg.norm(c-point) for c in channel_coords]
    
    closest_point_ind = np.argmin(distance)
    closest_point = channel_coords[closest_point_ind]
    displacement = point - closest_point
    
    return displacement, np.min(distance)

def get_annotation_file(mid: str, probe: str, day: int) -> str:
    path = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/Probe_{}_channels_{}_warped.csv'
    full_path = pathlib.Path(path.format(mid, probe + str(day), mid)).as_posix()

    return full_path

def update_ccf_and_calculation_tables(df_insertion_channel_info:pd.DataFrame, df_insertion_hit_rate:pd.DataFrame, df_min_displacement_vector:pd.DataFrame,
                                      df_min_distance:pd.DataFrame, area_center_of_mass:dict):
    for index, row in df_insertion_channel_info.iterrows():
        #if pd.isna(row.Channel_0_AP):
        annotation_file = get_annotation_file(row.MID, row.Probe, row.Day)
        try:
            df = pd.read_csv(annotation_file)
        except FileNotFoundError:
            logging.warning('{} not found'.format(annotation_file))
        else:
            for ccf_index, ccf_row in df.iterrows():
                df_insertion_channel_info.at[index,'Channel_{}_AP'.format(ccf_index)] = ccf_row.AP
                df_insertion_channel_info.at[index,'Channel_{}_DV'.format(ccf_index)] = ccf_row.DV
                df_insertion_channel_info.at[index,'Channel_{}_ML'.format(ccf_index)] = ccf_row.ML
            
            regions = df['region'].tolist()
            abbreviated_regions = strip_subregions_list(regions)
            ccf_coords = df_insertion_channel_info.iloc[index, 8:].to_numpy().reshape((384, 4))

            for area in area_center_of_mass:
                channel_hit_in_region = count_channels_in_region(abbreviated_regions, area)
                df_insertion_hit_rate.at[index, area] = channel_hit_in_region

                displacement, min_distance = min_distance_to_ccf_point(ccf_coords[:, 0:3], area_center_of_mass[area])
                df_min_distance.at[index, area] = min_distance
                df_min_displacement_vector.at[index, '{}_AP'.format(area)] = displacement[0]
                df_min_displacement_vector.at[index, '{}_DV'.format(area)] = displacement[1]
                df_min_displacement_vector.at[index, '{}_ML'.format(area)] = displacement[2]

    df_insertion_channel_info.to_sql('channel_ccf_coords', con=ENGINE, schema=None, if_exists='replace')
    df_min_distance.to_sql('min_distance_to_region', con=ENGINE, schema=None, if_exists='replace')
    df_min_displacement_vector.to_sql('vector_to_region_com', con=ENGINE, schema=None, if_exists='replace')
    df_insertion_hit_rate.to_sql('hit_region', con=ENGINE, schema=None, if_exists='replace')

if __name__ == '__main__':
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'r') as f:
        area_center_of_mass = json.load(f)

    df_insertion_channel_info = pd.read_sql_table('channel_ccf_coords', con=ENGINE)
    df_insertion_hit_rate = pd.read_sql_table('hit_region', con=ENGINE)
    df_min_displacement_vector = pd.read_sql_table('vector_to_region_com', con=ENGINE)
    df_min_distance = pd.read_sql_table('min_distance_to_region', con=ENGINE)

    update_ccf_and_calculation_tables(df_insertion_channel_info, df_insertion_hit_rate, df_min_displacement_vector, df_min_distance, area_center_of_mass)

