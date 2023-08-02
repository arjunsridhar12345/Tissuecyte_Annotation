import pandas as pd
import sqlite3
import pathlib
from sqlalchemy import  create_engine
import numpy as np

if __name__ == '__main__':
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
    """
    TABLES
        - session_metadata
        - vector_to_region_com
        - hit_region
        - channel_ccf_coords
        - min_distance_to_region
    """
    # read dr session table as a dataframe
    table_name = 'channel_ccf_coords'
    df_sessions = pd.read_sql_table(table_name, con=ENGINE)
    print(df_sessions)