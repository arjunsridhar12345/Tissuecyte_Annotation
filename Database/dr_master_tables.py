from __future__ import annotations

import contextlib
import csv
import datetime
import enum
import hashlib
#from lib2to3.pgen2 import driver
import pathlib
import sqlite3
from collections.abc import Sequence
from turtle import window_height
from typing import ClassVar, Optional
import json

from typing_extensions import Literal
import uuid
#import np_config
import numpy as np

from sqlalchemy import create_engine
import pandas as pd
import logging

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

def create_insertion_channel_table(path:pathlib.Path):
    master_excel = pd.read_excel(path)

    insertion_channel = {'MID': [], 'Day':[], 'Probe': [], 'Implant': [], 'Hole': [], 'Rig': [], 'Channel_annotation_file': []}
    for i in range(384):
        for position in ['AP', 'DV', 'ML', 'region']:
            insertion_channel['Channel_{}_{}'.format(i, position)] = []
        
        insertion_channel['Channel_{}_{}'.format(i, 'region')]

    for index, row in master_excel.iterrows():
        insertion_channel['MID'].append(row.MID)
        insertion_channel['Probe'].append(row.probe)
        insertion_channel['Day'].append(row.day)
        insertion_channel['Implant'].append(row.implant)
        insertion_channel['Hole'].append(row.hole)
        insertion_channel['Rig'].append(row.Rig)

        annotation_file = get_annotation_file(row.MID, row.probe, row.day)
        insertion_channel['Channel_annotation_file'].append(annotation_file)
        
        try:
            print(row.MID)
            df = pd.read_csv(annotation_file)
        except FileNotFoundError as e:
            for i in range(384):
                insertion_channel['Channel_{}_AP'.format(i)].append(np.nan)
                insertion_channel['Channel_{}_DV'.format(i)].append(np.nan)
                insertion_channel['Channel_{}_ML'.format(i)].append(np.nan)
                insertion_channel['Channel_{}_region'.format(i)].append('Track not annotated')
                #logging.warning(e)
        else:
            for index, r in df.iterrows():
                insertion_channel['Channel_{}_AP'.format(index)].append(r.AP)
                insertion_channel['Channel_{}_DV'.format(index)].append(r.DV)
                insertion_channel['Channel_{}_ML'.format(index)].append(r.ML)

                if pd.isna(r.region):
                    insertion_channel['Channel_{}_region'.format(index)].append('No Area')
                else:
                    insertion_channel['Channel_{}_region'.format(index)].append(r.region)

    df_channel_coords = pd.DataFrame(insertion_channel)
    df_channel_coords.to_sql('channel_ccf_coords', con=ENGINE, if_exists='replace')
"""
class DRMasterTable(Base):
    __tablename__ = 'dr_master'
    #id = mapped_column(Uuid(as_uuid=True), primary_key=True, default=uuid.uuid4)
    MID: Mapped[int]
    implant: Mapped[str] = mapped_column(primary_key=True)
    day: Mapped[int]
    probe: Mapped[str] = mapped_column(primary_key=True)
    hole: Mapped[Optional[str]] = mapped_column(primary_key=True)
    dye: Mapped[str]
    session: Mapped[int] = mapped_column(primary_key=True)
    tissuecyte: Mapped[bool]
    annotated: Mapped[bool]
    genotype: Mapped[str]
    dipped5x: Mapped[bool]
    retracted: Mapped[bool]
    probeandday: Mapped[str]
    probeloc_x: Mapped[Optional[float]]
    probeloc_y: Mapped[Optional[float]]
    channel_annotation_file: Mapped[str]
    rig_id: Mapped[Optional[str]]

    mean_channel_counts = relationship('mean_channel_counts', back_populates='dr_master')

    @classmethod
    def getChannelAnnotations(cls, df:pd.DataFrame) -> list:
        df_mouse_probes = df[['MID', 'probeandday']]
        annotation_files = []
        path = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/Probe_{}_channels_{}_warped.csv'

        for index, row in df_mouse_probes.iterrows():
            full_path = pathlib.Path(path.format(row.MID, row.probeandday, row.MID)).as_posix()
            annotation_files.append(full_path)

        return annotation_files

    @classmethod
    def getRigIDs(cls, df:pd.DataFrame) -> list:
        sessions = df['session'].tolist()
        rig_ids = []
        path = '//allen/programs/mindscope/workgroups/np-exp/{}/{}_platformD1.json'

        for session in sessions:
            with open(pathlib.Path(path.format(session, session))) as f:
                platform_json = json.load(f)
                rig_id = platform_json['rig_id']
                rig_ids.append(rig_id)

        return rig_ids

    @classmethod
    def insertDataFromExcel(cls, path:str):
        path_excel = pathlib.Path(path)
        df = pd.read_excel(path_excel)
        df.loc[df['probeloc_x'] == 'noprobeloc', 'probeloc_x'] = -1
        df.loc[df['probeloc_y'] == 'noprobeloc', 'probeloc_y'] = -1
        annotation_files = cls.getChannelAnnotations(df)
        rig_ids = cls.getRigIDs(df)
        df['channel_annotation_file'] = annotation_files
        df['rig_id'] = rig_ids

        return tuple(
            cls(
                **{
                    metric: df.loc[idx][metric]
                    for metric in df.columns
                }
            )
            for idx in df.index
        )


class MeanChannelCounts():
    __tablename__ = 'mean_channel_counts'

    @classmethod
    def populateTable(cls):
        df = pd.read_csv(pathlib.Path("//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/mean_channel_counts.csv"))
        df['Rig'] = ['NP.1' for i in range(len(df['Probe']))]
        df.to_sql('mean_channel_counts', con=ENGINE, schema=None, if_exists='replace')

class InsertionChannels():
    __table__name = 'insertion_channels'

    @classmethod
    def populateTable(cls):
        conn = create_connection(DB_PATH)
        query = 'SELECT dr.probe, dr.implant, dr.hole, dr.mid, dr.channel_annotation_file, mc.insertion_count \
                FROM DR_MASTER dr, MEAN_CHANNEL_COUNTS mc \
                WHERE dr.probe=mc.probe AND dr.implant=mc.implant AND dr.hole=mc.hole'

        cur = conn.cursor()
        cur.execute(query)

        rows = cur.fetchall()

        insertion_channel = {'Probe': [], 'Implant': [], 'Hole': [], 'Insertion_count': [], 'Channel_annotation_file': []}
        for i in range(384):
            for position in ['AP', 'DV', 'ML']:
                insertion_channel['Channel_{}_{}'.format(i, position)] = []

        for row in rows:
            insertion_channel['Probe'].append(row[0])
            insertion_channel['Implant'].append(row[1])
            insertion_channel['Hole'].append(row[2])
            insertion_channel['Insertion_count'].append(row[-1])
            insertion_channel['Channel_annotation_file'].append(row[-2])

            try:
                df = pd.read_csv(row[-2])
            except FileNotFoundError as e:
                for i in range(384):
                    insertion_channel['Channel_{}_AP'.format(i)].append(np.nan)
                    insertion_channel['Channel_{}_DV'.format(i)].append(np.nan)
                    insertion_channel['Channel_{}_ML'.format(i)].append(np.nan)
                logging.warning(e)
            else:
                for index, r in df.iterrows():
                    insertion_channel['Channel_{}_AP'.format(index)].append(r.AP)
                    insertion_channel['Channel_{}_DV'.format(index)].append(r.DV)
                    insertion_channel['Channel_{}_ML'.format(index)].append(r.ML)

        df_insertion = pd.DataFrame(insertion_channel)
        df_insertion.to_sql('insertion_channel_info', con=ENGINE, schema=None)
        conn.close()
"""
def md5(path: str | pathlib.Path) -> str:
    return hashlib.md5(pathlib.Path(path).read_bytes()).hexdigest()



if __name__ == "__main__":

    print('Yes')
    # with SESSION as session:
    #     # probe = NeuropixelsProbe(serial_number=18005117142)
    #     session.merge(Recording.dummy())
    #     session.merge(NeuropixelsProbe.dummy())
    #     session.add(SortedProbeRecording.dummy())
    #     session.add(LIMSEcephysSession.dummy())
    #     session.add_all(
    #         SortedUnit.from_csv_path(r'\\allen\programs\mindscope\workgroups\np-ultra\0_0_20230123\0_0_20230123_probeF_sorted\continuous\Neuropix-PXI-100.0\metrics.csv')
    #     )
    #     session.commit()


    # stmt = select(NeuropixelsProbe).where(NeuropixelsProbe.serial_number.in_([18005117142]))

    # for probe in SESSION.scalars(stmt):
    #     print(probe)
        
    """
    print(tuple(SESSION.scalars(select(LIMSEcephysSession).where(NeuropixelsProbe.serial_number.in_([18005117142])))))
    tuple(SESSION.scalars(select(SortedUnit).outerjoin(SortedProbeRecording)))
    
    # requires pandas ver that doesn't support 3.7
    df = pd.read_sql_table('sorted_units', con=ENGINE.connect(), schema=None)
    print(df)
    """
    
    #InsertionChannels.populateTable()
    #populate_session_table(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/dr_master_sheet_ccb.xlsx'))
    create_insertion_channel_table(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/dynamic_gating_insertions/dr_master_sheet_ccb.xlsx'))