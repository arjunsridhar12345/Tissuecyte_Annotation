from psycopg2 import connect, extras
    
# query lims and return result from given query
def query_lims(query_string):
    con = connect(
    dbname='lims2',
    user='limsreader',
    host='limsdb2',
    password='limsro',
    port=5432,
    )
    con.set_session(
        readonly=True, 
        autocommit=True,
    )
    cursor = con.cursor(
        cursor_factory=extras.RealDictCursor,
    )
    cursor.execute(query_string)
    result = cursor.fetchall()

    return result

# gets the tissuecyte info for the mouse id
def get_tc_info(mouse_id):
    TISSUECYTE_QRY = '''
            SELECT *
            FROM image_series im
            WHERE im.specimen_id = {}
        '''

    storage_directory = ''
    tc = query_lims(TISSUECYTE_QRY.format(get_specimen_id_from_labtracks_id(int(mouse_id))))
    for row in tc:
        d = dict(row)
        if d['alignment3d_id'] != None:
            storage_directory = d['storage_directory']

    print(storage_directory)
    return storage_directory

def get_specimen_id_from_labtracks_id(labtracks_id):
    SPECIMEN_QRY = '''
            SELECT *
            FROM specimens sp
            WHERE sp.external_specimen_name=cast({} as character varying)
        '''

    mouse_info = query_lims(SPECIMEN_QRY.format(int(labtracks_id)))
    return mouse_info[0]['id']

if __name__ == '__main__':
    get_tc_info('671439')