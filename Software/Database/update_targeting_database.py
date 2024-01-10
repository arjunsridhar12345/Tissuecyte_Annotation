import pandas as pd
import sqlite3
import pathlib
from sqlalchemy import  create_engine
import numpy as np
import json
import argparse
from dr_master_tables import get_annotation_file
from update_calculations_tables import update_min_distance_vector_displacement_tables, update_hit_rate_table
from typing import Union
import npc_session
from clean_region import assign_label
import pickle
import SimpleITK as sitk


with open(pathlib.Path(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\field_reference\acrnm_map.pkl"), 'rb') as f:
    ACRONYM_MAP = pickle.load(f)

ANNOTATION_VOLUME = sitk.GetArrayFromImage(sitk.ReadImage(pathlib.Path(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\field_reference\ccf_ano.mhd")))

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

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)
parser.add_argument('--genotype', help='Genotype of mouse', required=True)
parser.add_argument('--project', help='Project of experiment', required=True)
parser.add_argument('--rig', help='Experiment rig', required=True)
parser.add_argument('--date', help='Date of experiment', required=True)
parser.add_argument('--day', help='Experiment Day', required=True)

def get_insertion_path(project:str, date:str) -> pathlib.Path:
    insertion_base_path = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/ben')
    if project == 'Templeton':
        insertion_path = pathlib.Path(insertion_base_path, 'templeton', 'insertion_records')
    else:
        insertion_path = pathlib.Path(insertion_base_path, 'implants', 'insertion_records')
    
    insertion_json_path = list(insertion_path.glob('{}*'.format(date.replace('_', ''))))[0]

    return insertion_json_path

def get_dye_path(project:str, date:str, mouse_id:str) -> Union[pathlib.Path, None]:
    if project == 'Templeton':
        dye_base_path = pathlib.Path('//allen/programs/mindscope/workgroups/templeton/TTOC/pilot recordings')
        dye_json_path_list = list(dye_base_path.glob('{}*_{}/*_dye.json'.format(date.replace('_', '-'), mouse_id)))
        if len(dye_json_path_list) > 0:
            dye_json_path = dye_json_path_list[0]
        else:
            return None
    else:
        dye_base_path = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot')
        dye_json_path_list = list(dye_base_path.glob('*_{}_{}/*_dye.json'.format(mouse_id, date.replace('_', ''))))
        if len(dye_json_path_list) > 0:
            dye_json_path = dye_json_path_list[0]
        else:
            return None
    
    return dye_json_path

def update_session_metadata_table(mouse_id:str, genotype:str, project:str, rig:str, date:str, day:int):
    session_metadata = pd.read_sql_query('''SELECT session FROM session_metadata''', con=ENGINE)
    meta_data_index = len(session_metadata)
    insertion_json_path = get_insertion_path(project, date)
    with open(insertion_json_path, 'r') as f:
        insertion_json = json.load(f)
    
    dye_json_path = get_dye_path(project, date, mouse_id)
    if dye_json_path is not None:
        with open(dye_json_path, 'r') as f:
            dye_json = json.load(f)
            dye = dye_json['dye']
    else:
        dye = None
    

    session_metadata_dict:dict = {}
    date = date.replace('_', '')
    session_metadata_dict['session'] = [f'{mouse_id}_{date}']
    session_metadata_dict['MID'] = [mouse_id]
    session_metadata_dict['implant'] = '2002' if '2002' in [insertion_json['probe_insertions']['implant']] else [insertion_json['probe_insertions']['implant']]
    session_metadata_dict['day'] = [day]
    session_metadata_dict['dye'] = dye
    session_metadata_dict['Rig'] = [rig]
    session_metadata_dict['genotype'] = [genotype]
    session_metadata_dict['ProbeA'] = [insertion_json['probe_insertions']['probeA']['hole']]
    session_metadata_dict['ProbeB'] = [insertion_json['probe_insertions']['probeB']['hole']]
    session_metadata_dict['ProbeC'] = [insertion_json['probe_insertions']['probeC']['hole']]
    session_metadata_dict['ProbeD'] = [insertion_json['probe_insertions']['probeD']['hole']]
    session_metadata_dict['ProbeE'] = [insertion_json['probe_insertions']['probeE']['hole']]
    session_metadata_dict['ProbeF'] = [insertion_json['probe_insertions']['probeF']['hole']]
    annotation_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}'.format(mouse_id))

    if session_metadata_dict['session'] in session_metadata['session'].values:
        raise ValueError('Session {} already in table'.format(session_metadata['session']))
    
    if not annotation_path.exists():
        raise FileNotFoundError(f'{annotation_path.as_posix()} path not found')

    session_metadata_dict['annotation_directory'] = [annotation_path.as_posix()]

    df_row_to_add = pd.DataFrame(session_metadata_dict)
    df_row_to_add.index = [meta_data_index]
    df_row_to_add.to_sql('session_metadata', con=ENGINE, schema=None, if_exists='append')

def update_ccf_channel_table(mouse_id:str, day:int) -> pd.DataFrame:
    sessions_to_update = pd.read_sql_query(f'''SELECT * FROM session_metadata sm WHERE sm.MID = {mouse_id} AND sm.day = {day}''',
                                           con=ENGINE)
    ccf_mouse_ids = pd.read_sql_query(f'''SELECT ccf.MID, ccf.Day FROM channel_ccf_coords ccf WHERE ccf.MID = {mouse_id} AND ccf.day = {day}''',
                                      con=ENGINE)
    all_ccf_mouse_ids = pd.read_sql_query(f'''SELECT ccf.MID FROM channel_ccf_coords ccf''',
                                      con=ENGINE)
    print(sessions_to_update)

    if len(ccf_mouse_ids) > 0:
        raise ValueError('Entry already in database')
    
    if len(sessions_to_update) == 0:
        raise ValueError('No record to update')
    
    ccf_channel_row: dict = {'session': [],'MID': [], 'Day':[], 'Probe': [], 'Implant': [], 'Hole': [], 
                             'Rig': [], 'Channel_annotation_file': []}
    for i in range(384):
        for position in ['AP', 'DV', 'ML', 'region']:
            ccf_channel_row['Channel_{}_{}'.format(i, position)] = []
        
        ccf_channel_row['Channel_{}_{}'.format(i, 'region')]

    probes = [probe for probe in sessions_to_update.columns if 'Probe' in probe]

    for probe in probes:
        hole = sessions_to_update[probe].values[0]
        ccf_channel_row['session'].append(sessions_to_update['session'])
        ccf_channel_row['MID'].append(mouse_id)
        ccf_channel_row['Day'].append(day)
        ccf_channel_row['Probe'].append(probe[-1])
        ccf_channel_row['Implant'].append(sessions_to_update['implant'].values[0])
        ccf_channel_row['Hole'].append(hole)
        ccf_channel_row['Rig'].append(sessions_to_update['Rig'].values[0])
        annotation_file = get_annotation_file(mouse_id, probe[-1], day)
        ccf_channel_row['Channel_annotation_file'].append(annotation_file)

        try:
            print(mouse_id)
            df = pd.read_csv(annotation_file)
        except FileNotFoundError as e:
            for i in range(384):
                ccf_channel_row['Channel_{}_AP'.format(i)].append(np.nan)
                ccf_channel_row['Channel_{}_DV'.format(i)].append(np.nan)
                ccf_channel_row['Channel_{}_ML'.format(i)].append(np.nan)
                ccf_channel_row['Channel_{}_region'.format(i)].append('Track not annotated')
                #logging.warning(e)
        else:
            for index, r in df.iterrows():
                ccf_channel_row['Channel_{}_AP'.format(index)].append(r.AP)
                ccf_channel_row['Channel_{}_DV'.format(index)].append(r.DV)
                ccf_channel_row['Channel_{}_ML'.format(index)].append(r.ML)

                if pd.isna(r.region) or r.region == 'out of brain':
                    region = assign_label(ACRONYM_MAP, ANNOTATION_VOLUME, (r.AP, r.DV, r.ML))
                    ccf_channel_row['Channel_{}_region'.format(index)].append(region)
                else:
                    ccf_channel_row['Channel_{}_region'.format(index)].append(r.region)

    df_rows_to_add = pd.DataFrame(ccf_channel_row)
    index = range(len(all_ccf_mouse_ids), len(all_ccf_mouse_ids) + len(df_rows_to_add))
    df_rows_to_add.index = index
    print(df_rows_to_add)
    df_rows_to_add.to_sql('channel_ccf_coords', con=ENGINE, schema=None, if_exists='append')

    return df_rows_to_add

if __name__ == '__main__':
    args = parser.parse_args()
    with open(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/area_center_of_mass.json'), 'rb') as f:
        area_center_of_mass = json.load(f)

    structure_tree = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/ccf_structure_tree_2017.csv'))

    update_session_metadata_table(args.mouseID, args.genotype, args.project, args.rig, args.date, int(args.day))
    #ccf_rows_added = update_ccf_channel_table(args.mouseID, int(args.day))
    #update_min_distance_vector_displacement_tables(ccf_rows_added, area_center_of_mass, ENGINE)
    #update_hit_rate_table(ccf_rows_added, area_center_of_mass, structure_tree, ENGINE)


