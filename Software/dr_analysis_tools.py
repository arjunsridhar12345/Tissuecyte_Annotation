import os

dr_path = r'\\allen\programs\mindscope\workgroups\dynamicrouting\PilotEphys\Task 2 pilot'
dg_path = r'\\allen\programs\mindscope\workgroups\np-exp'

def get_exp_list_dr(mouse_id):
    os.chdir(dr_path)
    files = os.listdir(dr_path)
    exp_list = []
    for file in files:
        if file.__contains__(str(mouse_id)) and 'surface_channels' not in file:
            exp_list.append(file)
    return exp_list

def get_exp_list_dg(mouse_id):
    os.chdir(dg_path)
    files = os.listdir(dg_path)
    exp_list = []
    for file in files:
        if file.__contains__(str(mouse_id)) and len(file) == 26:
            exp_list.append(file)
    return exp_list

def get_exp_paths(mouse_id):
    os.chdir(dr_path)
    files = os.listdir(dr_path)
    path_list = []
    for file in files:
        if file.__contains__(str(mouse_id)) and 'surface_channels' not in file:
            path_list.append(os.path.join(dr_path, file))
    return path_list