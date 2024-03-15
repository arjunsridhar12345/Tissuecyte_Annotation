from __future__ import annotations

import pathlib
import requests
import npc_lims
import npc_session
import warnings
import argparse
import numpy as np
import pandas as pd
from get_correlation_plot import get_correlation_data
import npc_ephys
import re
import xmltodict
import upath
import io

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')

SAMPLING_RATE = 30000.
BASE_PATH = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/plots')
NUM_CHANNELS = 384

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
        download_response_metrics = npc_lims.get_codeocean_client().get_result_file_download_url(computation_id, pathlib.Path(folder['path'], 
                                                                                                      f'probe{probe}_metrics.csv').as_posix())

        download_response_depths = npc_lims.get_codeocean_client().get_result_file_download_url(computation_id, pathlib.Path(folder['path'], 
                                                                                                      f'probe{probe}_spike_depths.npy').as_posix())

        download_response_times = npc_lims.get_codeocean_client().get_result_file_download_url(computation_id, pathlib.Path(folder['path'], 
                                                                                                      f'probe{probe}_spike_times.npy').as_posix())
        
        if download_response_metrics.ok and download_response_depths.ok and download_response_times.ok:
            download_annotation_data(download_response_metrics, probe, session)
            download_annotation_data(download_response_depths, probe, session)
            download_annotation_data(download_response_times, probe, session)
        else:
            warnings.warn(f'Results could not be obtained for probe{probe} - Check codeocean capsule for details', stacklevel=2)

def run_capsule(session: npc_session.SessionRecord, capsule_id: str):
    raw_data_asset = npc_lims.codeocean.get_session_raw_data_asset(session)
    sorted_data_asset = npc_lims.codeocean.get_session_sorted_data_asset(session)
    sorted_data_surface_asset = None

    try:
        sorted_data_surface_asset = npc_lims.get_session_sorted_data_asset(session.with_idx(1))
    except (ValueError, IndexError):
        pass

    if sorted_data_surface_asset is None:
        capsule_run = npc_lims.codeocean.get_codeocean_client().run_capsule(capsule_id, [{'id': raw_data_asset['id'], 'mount': raw_data_asset['name']},
                                                                                {'id': sorted_data_asset['id'], 'mount': sorted_data_asset['name']}])
    else:
        capsule_run = npc_lims.codeocean.get_codeocean_client().run_capsule(capsule_id, [{'id': raw_data_asset['id'], 'mount': raw_data_asset['name']},
                                                                                {'id': sorted_data_asset['id'], 'mount': sorted_data_asset['name']},
                                                                                {'id': sorted_data_surface_asset['id'], 'mount': sorted_data_surface_asset['name']}])

    capsule_run.raise_for_status()
   
def get_capsule_results(capsule_id: str, session: npc_session.SessionRecord) -> None:   
    capusle_computations = npc_lims.get_codeocean_client().get_capsule_computations(capsule_id)
    capusle_computations.raise_for_status()
    capusle_computations = capusle_computations.json()

    for computation in capusle_computations:
        if not computation["has_results"] and computation["state"] != "completed" and computation['end_status'] != 'stopped':
            continue

        response_result_items = npc_lims.get_codeocean_client().get_list_result_items(
            computation["id"]
        )
        response_result_items.raise_for_status()
        result_items = response_result_items.json()

        session_result_item = tuple(
            item
            for item in result_items["items"]
            if re.match(  
                f"ecephys_{session.subject}_{session.date}_{npc_session.PARSE_TIME}",
                item["name"],
            )
        )

        if session_result_item:
            computation_id = computation['id']
            break

    result_items = npc_lims.get_codeocean_client().get_list_result_items(computation_id).json()['items']
    download_metrics_spike_data(result_items, computation_id, session.id)

def save_for_annotation(session: npc_session.SessionRecord, spike_times: np.ndarray,
                        spike_depths: np.ndarray, peak_channels: list[np.intp], probe: str) -> None:
    # maintain same structure as dynamic gating
    df_metrics = pd.DataFrame({'peak_channel': peak_channels})

    session_path = BASE_PATH / session.id
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
    df_metrics.to_csv(session_metrics_path / 'metrics.csv', index=False)
    np.save(session_metrics_path / 'spike_times.npy', spike_times)
    np.save(session_metrics_path / 'spike_depths.npy', spike_depths)

def get_channel_positions(settings_xml_path:upath.UPath) -> np.ndarray:
    with io.BytesIO(settings_xml_path.read_bytes()) as f:
        settings_dict = xmltodict.parse(f.read())
    
    electrode_xpositions = list(settings_dict['SETTINGS']['SIGNALCHAIN'][0]['PROCESSOR'][0]['EDITOR']['NP_PROBE'][0]['ELECTRODE_XPOS'].values())
    electrode_ypositions = list(settings_dict['SETTINGS']['SIGNALCHAIN'][0]['PROCESSOR'][0]['EDITOR']['NP_PROBE'][0]['ELECTRODE_YPOS'].values())

    electrode_xpositions = [int(electrode_xposition) for electrode_xposition in electrode_xpositions]
    electrode_ypositions = [int(electrode_yposition) for electrode_yposition in electrode_ypositions]

    channel_positions = [[electrode_xpositions[i], electrode_ypositions[i]] for i in range(len(electrode_xpositions))]
    return np.array(channel_positions)

def save_refinement_metrics(session: npc_session.SessionRecord, session_surface: npc_session.SessionRecord | None=None) -> None:
    spike_interface_data = npc_ephys.SpikeInterfaceKS25Data(session)
    probes = spike_interface_data.probes

    for probe in probes:
        peak_channels = list(npc_ephys.get_amplitudes_waveforms_channels_ks25(spike_interface_data, electrode_group_name=probe).peak_channels)
        spike_times = spike_interface_data.spike_indexes(probe) / SAMPLING_RATE
        spike_clusters = spike_interface_data.unit_indexes(probe)
        channel_positions = get_channel_positions(npc_lims.get_settings_xml_path_from_s3(session))

        if session_surface is not None:
            spike_interface_data_surface = npc_ephys.SpikeInterfaceKS25Data(session_surface)
            peak_channels_surface = list(npc_ephys.get_amplitudes_waveforms_channels_ks25(spike_interface_data_surface, electrode_group_name=probe).peak_channels)
            spike_times_surface = spike_interface_data_surface.spike_indexes(probe) / SAMPLING_RATE
            spike_clusters_surface = spike_interface_data_surface.unit_indexes(probe)
            channel_positions_surface = get_channel_positions(npc_lims.get_settings_xml_path_from_s3(session_surface))

            peak_channels_surface = [peak_channel + 384 for peak_channel in peak_channels_surface]
            spike_times_surface = spike_times_surface + (spike_times[-1] + 1)
            spike_clusters_surface = spike_clusters_surface + (np.max(spike_clusters) + 1)

            channel_positions = np.concatenate((channel_positions, channel_positions_surface))
            peak_channels = np.concatenate((np.array(peak_channels), np.array(peak_channels_surface)))
            spike_clusters = np.concatenate((spike_clusters, spike_clusters_surface))
            spike_times = np.concatenate((spike_times, spike_times_surface))

    
        cluster_difference = np.max(spike_clusters) - len(spike_interface_data.unit_locations(probe)) + 1
        if cluster_difference > 0:
            spike_clusters = spike_clusters - cluster_difference

        clusters_depths = channel_positions[peak_channels, 1]
        spike_depths = clusters_depths[spike_clusters]

        save_for_annotation(session, spike_times, spike_depths, peak_channels, probe)

def get_annotation_data_for_mouse(mouse_id:str) -> None:
    
    sessions = npc_lims.get_sessions_with_data_assets(mouse_id)
    for session in sessions:
        save_refinement_metrics(session)
    
    get_correlation_data(mouse_id)
    
if __name__ == '__main__':
    #args = parser.parse_args()
    #mouse_id = args.mouseID
    mouse_id = '667252'
    get_annotation_data_for_mouse(mouse_id)