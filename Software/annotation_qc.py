from PIL import Image, ImageEnhance
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import generate_metrics_paths
import cv2
import json
import npc_session

def get_surface_image_and_insertion_json_paths(session_id: str, templeton=False) -> tuple[pathlib.Path, ...]:
    date = session_id[session_id.rindex('_')+1:]
    session = npc_session.SessionRecord(session_id)
    base_path = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot')
    templeton_base_path = pathlib.Path('//allen/programs/mindscope/workgroups/templeton/TTOC/pilot recordings')
    insertion_records_directory = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/ben/implants/insertion_records')
    insertion_records_templeton = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/ben/templeton/insertion_records')

    
    if not templeton:
        mouse_dirs = sorted(generate_metrics_paths.get_metrics_directory(base_path, str(session.subject.id)))
        session_directory = base_path / [mouse_dir for mouse_dir in mouse_dirs if mouse_dir == f'DRpilot_{session_id}'][0]
    else:
        mouse_dirs = sorted(generate_metrics_paths.get_metrics_directory(templeton_base_path, str(session.subject.id)))
        session_directory = templeton_base_path / [mouse_dir for mouse_dir in mouse_dirs if session.date.id in mouse_dir][0]
    
    surface_image_path = next(session_directory.glob('*pre_insertion_surface_image.png'))

    if not templeton:
        insertions_json_path = next(insertion_records_directory.glob(f'{date}*'))
    else:
        insertions_json_path = next(insertion_records_templeton.glob(f'{date}*'))

    return surface_image_path, insertions_json_path


def draw_points_on_holes(image: np.ndarray, holes_letters: dict[str, str]) -> np.ndarray:
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    blur = cv2.GaussianBlur(gray, (9,9), 0)
    thresh = cv2.adaptiveThreshold(blur,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV,11,30)

    # Dilate to combine adjacent text contours
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (6,6))
    dilate = cv2.dilate(thresh, kernel, iterations=4)

    # Find contours, highlight text areas, and extract ROIs
    cnts = cv2.findContours(dilate, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cnts = cnts[0] if len(cnts) == 2 else cnts[1]
    index = 0
    print(len(cnts))
    ROI_number = 0
    holes = ['C1', 'D1', 'D2', 'D3', 'E2', 'B2', 'C2', 'E4', 'C4', 'C3', 'B1', 'E1', 'A1', 'E3', 'B4', 'F3', 'B3', 'F1', 
             'A2', 'F2', 'A3']
    
    hole_coordinates: dict = {}
    for c in cnts:
        area = cv2.contourArea(c)
        #if area > 10000:
        x,y,w,h = cv2.boundingRect(c)
        
        if index < len(holes):
            hole_coordinates[holes[index]] = (int((x + x + w) / 2), int((y + y + h) / 2))
        
        #cv2.rectangle(image, (x, y), (x + w, y + h), (36,255,12), 3)
            
        index += 1
    
    for hole in holes_letters:
        image = cv2.putText(image, holes_letters[hole], hole_coordinates[hole], cv2.FONT_HERSHEY_PLAIN, 
                   2, (255, 165, 0), 3, cv2.LINE_AA)
    
    return image

if __name__ == '__main__':
    surface_image_path, insertion_json_path = get_surface_image_and_insertion_json_paths('670181_20230718', templeton=True)
    implant_image = cv2.imread(pathlib.Path(r"\\allen\programs\mindscope\workgroups\dynamicrouting\arjun\2002_implant.png").as_posix())
    surface_image = Image.open(surface_image_path)
    enhanced_object = ImageEnhance.Brightness(surface_image)
    surface_image_enhanced = np.array(enhanced_object.enhance(2.0))
    holes_letters: dict = {}

    with open(insertion_json_path) as f:
        insertion_json = json.load(f)['probe_insertions']

    for probe in insertion_json:
        if 'probe' in probe:
            probe_dict = insertion_json[probe]
            holes_letters[probe_dict['hole']] = probe_dict['letter']

    drawing_implant = draw_points_on_holes(implant_image, holes_letters)

    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(surface_image_enhanced)
    ax[1].imshow(drawing_implant)
    plt.show()
    