import os
import pathlib
import requests
import npc_lims
import warnings
import argparse
from get_correlation_plot import get_correlation_data

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')

def download_annotation_data(download_response:requests.Response, probe: str, session:str):
    url = download_response.json()['url']
    get_response = requests.get(url,stream=True)
    file_name  = url.split("/")[-1]
    file_name = file_name.split('?')[0]
    file_name = file_name[file_name.index('_')+1:]

    base_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    session_path = base_path / session
    if not session_path.exists():
        session_path.mkdir()
    
    probe_path = session_path / f'probe{probe}'
    if not probe_path.exists():
        probe_path.mkdir()

    if not (probe_path / 'continuous').exists():
        (probe_path / 'continuous').mkdir()
    
    if not (probe_path / 'continuous' / 'Neuropix-PXI-100.0').exists():
        (probe_path / 'continuous' / 'Neuropix-PXI-100.0').mkdir()
    
    session_metrics_path = probe_path / 'continuous' / 'Neuropix-PXI-100.0'
    with open((session_metrics_path / file_name).as_posix(), 'wb') as f:
        for chunk in get_response.iter_content(chunk_size=1024):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)

def download_metrics_spike_data(result_items:list[dict[str, str]], computation_id:str, session:str) -> None:
    folder = [item for item in result_items if item['type'] == 'folder'][0]
    probes = ['A', 'B', 'C', 'D', 'E', 'F']

    for probe in probes:
        download_response_metrics = npc_lims.codeocean_client.get_result_file_download_url(computation_id, pathlib.Path(folder['path'], 
                                                                                                      f'probe{probe}_metrics.csv').as_posix())

        download_response_depths = npc_lims.codeocean_client.get_result_file_download_url(computation_id, pathlib.Path(folder['path'], 
                                                                                                      f'probe{probe}_spike_depths.npy').as_posix())

        download_response_times = npc_lims.codeocean_client.get_result_file_download_url(computation_id, pathlib.Path(folder['path'], 
                                                                                                      f'probe{probe}_spike_times.npy').as_posix())
        
        if download_response_metrics.ok and download_response_depths.ok and download_response_times.ok:
            download_annotation_data(download_response_metrics, probe, session)
            download_annotation_data(download_response_depths, probe, session)
            download_annotation_data(download_response_times, probe, session)
        else:
            warnings.warn(f'Results could not be obtained for probe{probe} - Check codeocean capsule for details', stacklevel=2)

def get_capsule_results(capsule_id: str, session_id:str) -> None:
    raw_data_asset = npc_lims.codeocean.get_session_raw_data_asset(session_id)
    sorted_data_asset = npc_lims.codeocean.get_session_sorted_data_asset(session_id)

    capusle_run = npc_lims.codeocean.codeocean_client.run_capsule(capsule_id, [{'id': raw_data_asset['id'], 'mount': raw_data_asset['name']},
                                                                               {'id': sorted_data_asset['id'], 'mount': sorted_data_asset['name']}])
    
    if not capusle_run.ok:
        raise RuntimeError('Annotation Capsule failed. Check codeocean')
    
    capusle_computations = npc_lims.codeocean_client.get_capsule_computations(capsule_id).json()
    while True:
        if 'end_status' in capusle_computations[0] and capusle_computations[0]['end_status'] == 'succeeded':
            break

        capusle_computations = npc_lims.codeocean.codeocean_client.get_capsule_computations(capsule_id).json()
    
    computation_id = capusle_computations[0]['id']
    result_items = npc_lims.codeocean_client.get_list_result_items(computation_id).json()['items']
    download_metrics_spike_data(result_items, computation_id, session_id)

def get_annotation_data_for_mouse(mouse_id:str, capsule_id:str):
    sessions = npc_lims.get_sessions_with_data_assets(mouse_id)
    for session in sessions:
        get_capsule_results(capsule_id, session.id)
    
    get_correlation_data(mouse_id)

if __name__ == '__main__':
    args = parser.parse_args()
    mouse_id = args.mouseID
    get_annotation_data_for_mouse(mouse_id, '6c4dad63-7fdf-4dfe-82f1-f9b24e924d31')