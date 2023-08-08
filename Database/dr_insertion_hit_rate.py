import numpy as np
import pandas as pd
import nrrd
from scipy import ndimage
import pathlib
import json
import sqlite3
from sqlalchemy import create_engine

# ------------------------------------------------------------------------------------
# allow integers >8 bytes to be stored in sqlite3
sqlite3.register_adapter(np.int64, lambda val: int(val))
sqlite3.register_adapter(np.int32, lambda val: int(val))
# ------------------------------------------------------------------------------------

DB_PATH = pathlib.Path("//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/dr_master.db")
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
    if not isinstance(areastr, str) or areastr == 'No Area' or areastr == 'Track not annotated':
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

def insertion_hit_rate_table(df_insertion_channel_coords:pd.DataFrame):
    dict_insertion_hit_rate = {'MID': [], 'Day': [], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': []}

    for area in area_center_of_mass:
        dict_insertion_hit_rate[area] = []

    for index, row in df_insertion_channel_coords.iterrows():
        try:
            df_annotations = pd.read_csv(row.Channel_annotation_file)
        except FileNotFoundError:
            for area in area_center_of_mass:
                dict_insertion_hit_rate[area].append(np.nan)
        else:
            regions = df_annotations['region'].tolist()

            regions_abbreviated = strip_subregions_list(regions)

            for area in area_center_of_mass:
                number_of_hits = count_channels_in_region(regions_abbreviated, area)
                dict_insertion_hit_rate[area].append(number_of_hits)

        dict_insertion_hit_rate['MID'].append(row.MID)
        dict_insertion_hit_rate['Day'].append(row.Day)
        dict_insertion_hit_rate['Probe'].append(row.Probe)
        dict_insertion_hit_rate['Implant'].append(row.Implant)
        dict_insertion_hit_rate['Hole'].append(row.Hole)
        dict_insertion_hit_rate['Rig'].append(row.Rig)

    df_insertion_hit_rate = pd.DataFrame(dict_insertion_hit_rate)
    df_insertion_hit_rate.to_sql('hit_region', con=ENGINE, schema=None, if_exists='replace')

def check_insertion_table(df_insertion:pd.DataFrame):
    df_probe_hole_implant = df_insertion.drop(columns=['index', 'MID', 'Day'])
    df_probe_hole_implant_test = df_probe_hole_implant.loc[(df_probe_hole_implant.Probe == 'A') & (df_probe_hole_implant.Implant == 'TS5')
                                                           & (df_probe_hole_implant.Hole == 'C') & (df_probe_hole_implant.Rig == 'NP1')]
    print(df_probe_hole_implant_test.groupby(['Probe', 'Implant', 'Hole', 'Rig']).mean()['ACAd'])

if __name__ == '__main__':
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'r') as f:
        area_center_of_mass = json.load(f)

    df_insertion_channel_coords = pd.read_sql_table('channel_ccf_coords', con=ENGINE.connect())
    insertion_hit_rate_table(df_insertion_channel_coords)
    #df_insertion_hit_rate = pd.read_sql_table('hit_rate_by_insertion', con=ENGINE.connect(), schema=None)
    #check_insertion_table(df_insertion_hit_rate)
