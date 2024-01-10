import pandas as pd
import numpy as np
import json
import pathlib
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

def min_distance_to_ccf_point(group:np.array, point:list):
    channel_coords = group*25
    distance = [np.linalg.norm(c-point) for c in channel_coords]
    
    closest_point_ind = np.argmin(distance)
    closest_point = channel_coords[closest_point_ind]
    displacement = point - closest_point
    
    return displacement, np.min(distance)

def min_displacement_vector_table():
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'r') as f:
        area_center_of_mass = json.load(f)

    df_insertion_channel_coords = pd.read_sql_table('channel_ccf_coords', con=ENGINE.connect()).to_numpy()
    #print(df_insertion_channel_coords)
    dict_closest_point_insertion = {'session': [], 'MID': [], 'Day': [], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': []}
    directions = ['AP', 'DV', 'ML']

    for area in area_center_of_mass:
        for direction in directions:
            dict_closest_point_insertion['{}_{}'.format(area, direction)] = []

    for insertion in df_insertion_channel_coords:
        ccf_coords = insertion[9:].reshape((384, 4))
        
        for area in area_center_of_mass:
            displacement_vector = min_distance_to_ccf_point(ccf_coords[:, 0:3], area_center_of_mass[area])[0]

            for i in range(3):
                dict_closest_point_insertion['{}_{}'.format(area, directions[i])].append(displacement_vector[i])

        dict_closest_point_insertion['session'].append(insertion[1])
        dict_closest_point_insertion['MID'].append(insertion[2])
        dict_closest_point_insertion['Day'].append(insertion[3])
        dict_closest_point_insertion['Probe'].append(insertion[4])
        dict_closest_point_insertion['Implant'].append(insertion[5])
        dict_closest_point_insertion['Hole'].append(insertion[6])
        dict_closest_point_insertion['Rig'].append(insertion[7])

    df_min_distance_insertion = pd.DataFrame(dict_closest_point_insertion)
    df_min_distance_insertion.to_sql('vector_to_region_com', con=ENGINE, schema=None, if_exists='replace')

def min_distance_table():
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'r') as f:
        area_center_of_mass = json.load(f)

    df_insertion_channel_coords = pd.read_sql_table('channel_ccf_coords', con=ENGINE.connect()).to_numpy()

    dict_min_distance_insertion = {'session': [], 'MID': [], 'Day': [], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': []}

    for area in area_center_of_mass:
        dict_min_distance_insertion[area] = []

    for insertion in df_insertion_channel_coords:
        ccf_coords = insertion[9:].reshape((384, 4))
        
        for area in area_center_of_mass:
            displacement = min_distance_to_ccf_point(ccf_coords[:, 0:3], area_center_of_mass[area])[1]

            dict_min_distance_insertion[area].append(displacement)

        dict_min_distance_insertion['session'].append(insertion[1])
        dict_min_distance_insertion['MID'].append(insertion[2])
        dict_min_distance_insertion['Day'].append(insertion[3])
        dict_min_distance_insertion['Probe'].append(insertion[4])
        dict_min_distance_insertion['Implant'].append(insertion[5])
        dict_min_distance_insertion['Hole'].append(insertion[6])
        dict_min_distance_insertion['Rig'].append(insertion[7])

    df_min_distance_insertion = pd.DataFrame(dict_min_distance_insertion)
    df_min_distance_insertion.to_sql('min_distance_to_region', con=ENGINE, schema=None, if_exists='replace')

def check_min_distance_table(df_min_distance:pd.DataFrame):
    df_probe_hole_implant = df_min_distance.drop(columns=['index', 'MID', 'Day'])
    df_probe_hole_implant_test = df_probe_hole_implant.loc[(df_probe_hole_implant.Probe == 'B') & (df_probe_hole_implant.Implant == 'football')
                                                           & (df_probe_hole_implant.Hole == 'G') & (df_probe_hole_implant.Rig == 'NP1')]
    print(df_probe_hole_implant_test.groupby(['Probe', 'Implant', 'Hole', 'Rig']).median()['LGd'])

def check_min_displacement_table(df_min_displacement:pd.DataFrame):
    df_probe_hole_implant = df_min_displacement.drop(columns=['index', 'MID', 'Day'])
    df_probe_hole_implant_test = df_probe_hole_implant.loc[(df_probe_hole_implant.Probe == 'B') & (df_probe_hole_implant.Implant == 'football')
                                                           & (df_probe_hole_implant.Hole == 'G') & (df_probe_hole_implant.Rig == 'NP1')]

    area = 'LGd'
    print(df_probe_hole_implant_test.groupby(['Probe', 'Implant', 'Hole', 'Rig']).median()['{}_AP'.format(area)])
    print(df_probe_hole_implant_test.groupby(['Probe', 'Implant', 'Hole', 'Rig']).median()['{}_DV'.format(area)])
    print(df_probe_hole_implant_test.groupby(['Probe', 'Implant', 'Hole', 'Rig']).median()['{}_ML'.format(area)])

if __name__ == '__main__':
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'r') as f:
        area_center_of_mass = json.load(f)

    df_insertion_channel_coords = pd.read_sql_table('channel_ccf_coords', con=ENGINE.connect()).to_numpy()
    #min_distance_table(df_insertion_channel_coords)
    min_displacement_vector_table()
    #df_min_displacement = pd.read_sql_table('closest_point_by_insertion', con=ENGINE.connect(), schema=None)
    #check_min_displacement_table(df_min_displacement)