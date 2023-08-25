from __future__ import annotations

import pathlib
import sqlite3
from turtle import window_height
import numpy as np

from sqlalchemy import create_engine
import pandas as pd
import logging
import pickle
import SimpleITK as sitk
from clean_region import assign_label

with open(pathlib.Path(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\field_reference\acrnm_map.pkl"), 'rb') as f:
    ACRONYM_MAP = pickle.load(f)

ANNOTATION_VOLUME = sitk.GetArrayFromImage(sitk.ReadImage(pathlib.Path(r"\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\field_reference\ccf_ano.mhd")))

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


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        logging.info('connection succesful')
    except sqlite3.Error as e:
        logging.warning(e)

    return conn
"""
class Base(DeclarativeBase):
    pass
"""
"""
class MultiSessionRecording(Base):
    # week_start_date
    pass

class Mouse(Base):
    # implant
    pass
"""
def get_channel_annotations(df:pd.DataFrame) -> list:
    df_mouse_probes = df[['MID', 'probe', 'day']]
    annotation_files = []
    path = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/Probe_{}_channels_{}_warped.csv'

    for index, row in df_mouse_probes.iterrows():
        full_path = pathlib.Path(path.format(row.MID, row.probe + str(row.day), row.MID)).as_posix()
        annotation_files.append(full_path)

    return annotation_files

def get_annotation_directory(mid: str) -> str:
    path = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}'
    full_path = pathlib.Path(path.format(mid)).as_posix()

    return full_path

def get_annotation_file(mid: str, probe: str, day: int) -> str:
    path = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/Probe_{}_channels_{}_warped.csv'
    full_path = pathlib.Path(path.format(mid, probe + str(day), mid)).as_posix()

    return full_path

def populate_session_table(path:pathlib.Path):
    dict_session_table = {'session': [], 'MID': [], 'implant': [], 'day': [], 'dye': [], 'Rig': [], 'genotype': [], 
                          'ProbeA': [], 'ProbeB': [], 'ProbeC': [], 'ProbeD': [], 
                          'ProbeE': [], 'ProbeF': [], 'annotation_directory': []}

    master_excel = pd.read_excel(path)

    unique_mids = master_excel['MID'].unique()
    for mid in unique_mids:
        df_mid = master_excel.loc[master_excel['MID'] == mid]
        days = df_mid['day'].unique()

        for day in days:
            df_mid_day = df_mid.loc[df_mid['day'] == day]

            for probe in ['A', 'B', 'C', 'D', 'E', 'F']:
                print(probe)
                df_mid_day_probe = df_mid_day.loc[df_mid_day['probe'] == probe]

                if not df_mid_day_probe.empty:
                    dict_session_table['Probe{}'.format(probe)].append(df_mid_day_probe['hole'].values[0])
                else:
                    dict_session_table['Probe{}'.format(probe)].append(np.nan)

            dict_session_table['session'].append(df_mid_day['session'].values[0])
            dict_session_table['MID'].append(df_mid_day['MID'].values[0])
            dict_session_table['implant'].append(df_mid_day['implant'].values[0])
            dict_session_table['day'].append(df_mid_day['day'].values[0])
            dict_session_table['dye'].append(df_mid_day['dye'].values[0])
            dict_session_table['Rig'].append(df_mid_day['Rig'].values[0])
            dict_session_table['genotype'].append(df_mid_day['genotype'].values[0])
            dict_session_table['annotation_directory'].append(get_annotation_directory(df_mid_day['MID'].values[0]))

    print([[key, len(dict_session_table[key])] for key in dict_session_table])
    df_session = pd.DataFrame(dict_session_table)
    df_session.to_sql('session_metadata', con=ENGINE, schema=None, if_exists='replace')

def create_session_table():
    dummy_excel = pd.read_excel('dummy_session.xlsx')

    dummy_excel.to_sql('dr_session', con=ENGINE, schema=None)

def update_session_table():
    session_metadata = pd.read_sql_table('session_metadata', con=ENGINE)
    session_metadata.drop(columns=['level_0', 'index'], inplace=True)
    print(session_metadata)
    """
    
    for index, row in session_metadata.iterrows():
        session = row.session
        session_id = session[0:session.index('_')]

        if len(session_id) == 10:
            session_metadata.at[index, 'session'] = session_id
    """
    session_metadata.to_sql('session_metadata', con=ENGINE, schema=None, if_exists='replace')

def create_insertion_channel_table():
    master_table = pd.read_sql_table('session_metadata', con=ENGINE)

    insertion_channel = {'session': [], 'MID': [], 'Day':[], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': [], 'Channel_annotation_file': []}
    for i in range(384):
        for position in ['AP', 'DV', 'ML', 'region']:
            insertion_channel['Channel_{}_{}'.format(i, position)] = []
        
        insertion_channel['Channel_{}_{}'.format(i, 'region')]

    for index, row in master_table.iterrows():
        for probe in ['A', 'B', 'C', 'D', 'E', 'F']:
            insertion_channel['session'].append(row.session)
            insertion_channel['MID'].append(row.MID)
            insertion_channel['Probe'].append(probe)
            insertion_channel['Day'].append(row.day)
            insertion_channel['Implant'].append(row.implant)
            insertion_channel['Hole'].append(row[f'Probe{probe}'])
            insertion_channel['Rig'].append(row.Rig)

            annotation_file = get_annotation_file(row.MID, probe, row.day)
            
            try:
                print(row.MID)
                df = pd.read_csv(annotation_file)
            except FileNotFoundError as e:
                for i in range(384):
                    insertion_channel['Channel_{}_AP'.format(i)].append(np.nan)
                    insertion_channel['Channel_{}_DV'.format(i)].append(np.nan)
                    insertion_channel['Channel_{}_ML'.format(i)].append(np.nan)
                    insertion_channel['Channel_{}_region'.format(i)].append('Track not annotated')
                insertion_channel['Channel_annotation_file'].append('No annotation file')
                    #logging.warning(e)
            else:
                for index, r in df.iterrows():
                    insertion_channel['Channel_{}_AP'.format(index)].append(r.AP)
                    insertion_channel['Channel_{}_DV'.format(index)].append(r.DV)
                    insertion_channel['Channel_{}_ML'.format(index)].append(r.ML)

                    if pd.isna(r.region):
                        region = assign_label(ACRONYM_MAP, ANNOTATION_VOLUME, (r.AP, r.DV, r.ML))
                        insertion_channel['Channel_{}_region'.format(index)].append(region)
                    else:
                        insertion_channel['Channel_{}_region'.format(index)].append(r.region)

                insertion_channel['Channel_annotation_file'].append(annotation_file)

    df_channel_coords = pd.DataFrame(insertion_channel)
    df_channel_coords.to_sql('channel_ccf_coords', con=ENGINE, schema=None, if_exists='replace')

if __name__ == "__main__":
    #populate_session_table(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/dr_master_sheet_ccb.xlsx'))
    create_insertion_channel_table()
    #update_session_table()