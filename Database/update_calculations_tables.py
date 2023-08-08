import pandas as pd
import sqlite3
import pathlib
import numpy as np

def min_distance_to_ccf_point(group:np.ndarray, point:list):
    channel_coords = group*25
    distance = [np.linalg.norm(c-point) for c in channel_coords]
    
    closest_point_ind = np.argmin(distance)
    closest_point = channel_coords[closest_point_ind]
    displacement = point - closest_point
    
    return displacement, np.min(distance)

def strip_subregions_list(areas: list, structure_tree:pd.DataFrame):
    areas_abbreviated = []

    for area in areas:
        areas_abbreviated.append(strip_subregions_layers(area, structure_tree))

    return areas_abbreviated

def strip_subregions_layers(areastr, structure_tree:pd.DataFrame):
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

def update_hit_rate_table(ccf_coordinates_to_add:pd.DataFrame, area_center_of_mass:dict, structure_tree:pd.DataFrame, engine):
    dict_insertion_hit_rate:dict = {'MID': [], 'Day': [], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': []}

    for area in area_center_of_mass:
        dict_insertion_hit_rate[area] = []

    for index, row in ccf_coordinates_to_add.iterrows():
        try:
            df_annotations = pd.read_csv(row.Channel_annotation_file)
        except FileNotFoundError:
            for area in area_center_of_mass:
                dict_insertion_hit_rate[area].append(np.nan)
        else:
            regions = df_annotations['region'].tolist()

            regions_abbreviated = strip_subregions_list(regions, structure_tree)

            for area in area_center_of_mass:
                number_of_hits = count_channels_in_region(regions_abbreviated, area)
                dict_insertion_hit_rate[area].append(number_of_hits)

        dict_insertion_hit_rate['MID'].append(row.MID)
        dict_insertion_hit_rate['Day'].append(row.Day)
        dict_insertion_hit_rate['Probe'].append(row.Probe)
        dict_insertion_hit_rate['Implant'].append(row.Implant)
        dict_insertion_hit_rate['Hole'].append(row.Hole)
        dict_insertion_hit_rate['Rig'].append(row.Rig)
    
    df_hit_region_rows_to_add = pd.DataFrame(dict_insertion_hit_rate)
    df_hit_region_rows_to_add.index = ccf_coordinates_to_add.index
    df_hit_region_rows_to_add.to_sql('hit_region', con=engine, schema=None, if_exists='append')

def update_min_distance_vector_displacement_tables(ccf_coordinates_to_add:pd.DataFrame, area_center_of_mass:dict, engine):
    dict_min_distance_insertion:dict = {'MID': [], 'Day': [], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': []}
    dict_closest_point_insertion:dict = {'MID': [], 'Day': [], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': []}
    directions = ['AP', 'DV', 'ML']

    for area in area_center_of_mass:
        for direction in directions:
            dict_closest_point_insertion['{}_{}'.format(area, direction)] = []

        dict_min_distance_insertion[area] = []

    for insertion in ccf_coordinates_to_add.to_numpy():
        ccf_coords = insertion[7:].reshape((384, 4))
        
        for area in area_center_of_mass:
            displacement_vector, displacement = min_distance_to_ccf_point(ccf_coords[:, 0:3], area_center_of_mass[area])
            dict_min_distance_insertion[area].append(displacement)

            for i in range(3):
                dict_closest_point_insertion['{}_{}'.format(area, directions[i])].append(displacement_vector[i])
        
        dict_min_distance_insertion['MID'].append(insertion[0])
        dict_min_distance_insertion['Day'].append(insertion[1])
        dict_min_distance_insertion['Probe'].append(insertion[2])
        dict_min_distance_insertion['Implant'].append(insertion[3])
        dict_min_distance_insertion['Hole'].append(insertion[4])
        dict_min_distance_insertion['Rig'].append(insertion[5])

        dict_closest_point_insertion['MID'].append(insertion[0])
        dict_closest_point_insertion['Day'].append(insertion[1])
        dict_closest_point_insertion['Probe'].append(insertion[2])
        dict_closest_point_insertion['Implant'].append(insertion[3])
        dict_closest_point_insertion['Hole'].append(insertion[4])
        dict_closest_point_insertion['Rig'].append(insertion[5])
    
    df_min_distance_rows_to_add = pd.DataFrame(dict_min_distance_insertion)
    df_min_displacement_vector_rows_to_add = pd.DataFrame(dict_closest_point_insertion)
    df_min_distance_rows_to_add.index = ccf_coordinates_to_add.index
    df_min_displacement_vector_rows_to_add.index = ccf_coordinates_to_add.index

    df_min_distance_rows_to_add.to_sql('min_distance_to_region', con=engine, schema=None, if_exists='append')
    df_min_displacement_vector_rows_to_add.to_sql('vector_to_region_com', con=engine, schema=None, if_exists='append')
