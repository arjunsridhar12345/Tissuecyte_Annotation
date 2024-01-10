import numpy as np

def assign_label(acronym_map:dict[str, int], annotation_volume:np.ndarray, point:tuple[int, int, int]) -> str:
    structure_ids = tuple(acronym_map.values())
    labels = tuple(acronym_map.keys())

    structure_id = annotation_volume[point[0], point[1], point[2]]

    if point[1] < 0:
        label = 'out of brain'
    elif structure_id in structure_ids:
        index = structure_ids.index(structure_id)
        label = labels[index]
    else:
        label = 'root'
    
    return label


