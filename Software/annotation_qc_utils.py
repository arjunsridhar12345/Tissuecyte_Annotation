import pandas as pd
import sqlite3
import pathlib
from sqlalchemy import  create_engine
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import pathlib
from dr_analysis_tools import get_exp_paths
import math

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

NUM_DAYS = 4
PROBES = ['A', 'B', 'C', 'D', 'E', 'F']

##in bregma (anything in mm will do)
HOLE_COORDS_2002 = {
    'A1': [-.75, -.3],
    'A2': [-.95, 2.05],
    'A3': [-1.25, 2.5], #guess
    'B1': [-.92, -1.19],
    'B2': [-.96, -3.1],
    'B3': [-.75, .9],
    'B4': [-1.5, .325],
    'C1': [-1.53, -3.95],
    'C2': [-1.59, -2.39],
    'C3': [-2.3, -1.9],
    'C4': [-.7, -2.05],
    'D1': [-3.85, -3.7],
    'D2': [-3.08, -3.55],
    'D3': [-2.2, -3.2],
    'E1': [-3.44, -.98],
    'E2': [-3.56, -2.83],
    'E3': [-3.25, 0.0],
    'E4': [-3.6, -1.9],
    'F1': [-2.6, 1.8],
    'F2': [-1.74, 2.1],
    'F3': [-2.8, .9],
    None : [np.nan, np.nan], 
    'default': [0, 0]
}


##only works if mouse has probe_insertions jsons
def insertion_holes_from_json(mouse, probes):
    paths = get_exp_paths(mouse)
    json_list = []
    for path in paths:
        p_list = os.listdir(path)
        for p in p_list:
            if 'probe_insertions' in p:
                json_list.append(os.path.join(path, p))
    
    insertions = {}
    for j in json_list:
        with open(j) as file:
            data = json.load(file)
        for probe in probes:
            hole = data['probe_insertions']['probe' + probe]['hole']
            insertions[probe + str(data['probe_insertions']['day_1-4'])] = hole
    return insertions

##only works if mouse is in channel_ccf_coords table
def insertion_holes_from_db(mouse, df_sessions):
    mouse_df = df_sessions[df_sessions['MID'] == mouse]
    insertions = {}
    for index, row in mouse_df.iterrows():
        insertions[row['Probe'] + str(row['Day'])] = row['Hole']
    return insertions

def insertion_holes_from_db_metadata(df_sessions_metadata: pd.DataFrame) -> dict:
    insertions:dict = {}
    
    for index, row in df_sessions_metadata.iterrows():
        if row['implant'] != '2002':
            raise ValueError(f"Implant {row['implant']} different from 2002")
        
        for probe in PROBES:
            probe_day = probe+str(row['day'])
            insertions[probe_day] = row[f"Probe{probe}"]
    print(insertions)
    return insertions

##returns 'probeDay - probeDay' : [x, y] for each probe's implant hole across days
def get_implant_vectors(insertions:dict, implant='2002') -> dict:
    for probe in PROBES:
        for i in range(4):
            if probe+str(i+1) not in insertions:
                insertions[probe+str(i)] = 'default'

    implant_coords = {}
    implant_vectors = {} 
    for k, v in insertions.items():
        #for multiple implants, 
        #hole_dict = 'hole_coords_' + implant
        #implant_coords[k] = np.array(hole_dict[insertions[k]])
        implant_coords[k] = np.array(HOLE_COORDS_2002[insertions[k]])

    for probe in PROBES:
        probedict = {}
        for k, v in implant_coords.items():
            if probe in k:
                probedict[k] = v
            for i in range(len(probedict)):
                for k,v in probedict.items():
                    implant_vectors[probe + str(i+1) + '-' + k] = v - probedict[probe + str(i+1)]
    return implant_vectors

def generate_line(probe_annotations_day:pd.DataFrame) -> np.ndarray:
    # 25 micron affine space
    x = probe_annotations_day.ML / 2.5
    y = probe_annotations_day.DV / 2.5
    z = probe_annotations_day.AP / 2.5

    # get trajectory
    if len(z) > 0:
        data = np.vstack((z,y,x)).T
        datamean = data.mean(axis=0)
        D = data - datamean
        m1 = np.min(D[:,1]) * 2
        m2 = np.max(D[:,1]) * 2
        uu,dd,vv = np.linalg.svd(D)

        linepts = vv[0] * np.mgrid[-200:200:1][:,np.newaxis]
        linepts += datamean
    
        if linepts[-1,1] - linepts[0,1] < 0:
            linepts = np.flipud(linepts)
    
    return linepts

def get_surface_coords(probe_annotations: pd.DataFrame) -> dict:
    probes_annotated = probe_annotations['probe_name'].unique()
    probes = ['A', 'B', 'C', 'D', 'E', 'F']
    surface_channel_coords: dict = {}

    for probe in probes:
        for i in range(1, NUM_DAYS+1):
            probe_day_key = f'Probe {probe+str(i)}'
            probe_day=probe+str(i)
            if probe_day_key in probes_annotated:
                probe_annotations_day = probe_annotations[probe_annotations['probe_name'] == probe_day_key]
                extrapolated_line = generate_line(probe_annotations_day)
                surface_coord = extrapolated_line[extrapolated_line[:, 1] > 0][0]
                surface_channel_coords[probe_day] = (np.array([(surface_coord[2]/40), -(surface_coord[0]/40)]), 'True')
            else:
                surface_channel_coords[probe_day] = (np.array([-1, -1]), 'No annotation file')
    
    return surface_channel_coords

##returns 'probeDay - probeDay' : [x, y] vectors for each probe's annotation across days
def get_annotation_vectors(coords:dict) -> dict:
    ann_vectors:dict = {}
    for probe in PROBES:
        probedict:dict = {}
        for k, v in coords.items():
            if probe in k:
                probedict[k] = v
        for i in range(4):
            for k,v in probedict.items():
                if probe + str(i+1) in probedict:
                    ann_vectors[probe + str(i+1) + '-' + k] = (v[0] - probedict[probe + str(i+1)][0], v[1])
                else:
                    ann_vectors[probe + str(i+1) + '-' + k] = (np.array([-1, -1]), v[1])
    return ann_vectors

def plot_vectors_arjun(annotation_vectors, implant_vectors, surface_coords):
    for probe in PROBES:
        vector_probes = np.array([vector for vector in list(annotation_vectors.keys()) if probe in vector]).reshape((4, 4))

        x_vectors = [annotation_vectors[vector][0][0] for vector in annotation_vectors.keys()]
        y_vectors = [annotation_vectors[vector][0][1] for vector in annotation_vectors.keys()]

        fig = plt.figure(figsize=(25,25))
        gs = fig.add_gridspec(ncols = 4, nrows = 4)

        for row in range(vector_probes.shape[0]):
            for column in range(vector_probes.shape[1]):
                vector_key = vector_probes[row, column]
                index = vector_key.index('-')
                first_probe = vector_key[0:index]
                second_probe = vector_key[index+1:]

                annotation_file_first_probe = surface_coords[first_probe][1]
                annotation_file_second_probe = surface_coords[second_probe][1]
                ax = fig.add_subplot(gs[row, column])
                ax.set_title(vector_key)
                ax.set_xlim(min(x_vectors)-0.5, max(x_vectors)+0.5)
                ax.set_ylim(min(x_vectors)-0.5, max(x_vectors)+0.5)
                if annotation_file_first_probe != 'No annotation file' and annotation_file_second_probe != 'No annotation file':
                     ax.quiver([0, 0], [0, 0], [annotation_vectors[vector_key][0][0], implant_vectors[vector_key][0]],
                               [annotation_vectors[vector_key][0][1], implant_vectors[vector_key][1]],
                               angles='xy', scale_units='xy', scale=1, color=['r', 'b'])
                else:
                    ax.text(-2, 0, 'Nothing', size=18)
        
        fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    mouse_id = 660023
    annotation_path = pathlib.Path(f'//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{mouse_id}/probe_annotations_{mouse_id}.csv')
    probe_annotations = pd.read_csv(annotation_path)
    df_sessions_metadata = pd.read_sql_table('session_metadata', con=ENGINE.connect())
    df_sessions_metadata = df_sessions_metadata[df_sessions_metadata['MID'] == mouse_id]

    insertions_holes_db = insertion_holes_from_db_metadata(df_sessions_metadata)
    implant_vectors = get_implant_vectors(insertions_holes_db)
    coords = get_surface_coords(probe_annotations)
    ann_vectors = get_annotation_vectors(coords)

    plot_vectors_arjun(ann_vectors, implant_vectors, coords)
    plt.show()
