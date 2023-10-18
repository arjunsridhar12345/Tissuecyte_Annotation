import pathlib
import matplotlib.pyplot as plt
import SimpleITK as sitk
import numpy as np

DEFAULT_COLOR_VALUES = [[0, 3000], [0, 3000], [0, 1000]]

def getColorVolume(int_arrays, rgb_levels=DEFAULT_COLOR_VALUES):
    level_adjusted_arrays = []
    for colori, int_level in zip(['red', 'green', 'blue'], rgb_levels):
        colarray = np.clip(int_arrays[colori], a_min=int_level[0], a_max=int_level[1]) - int_level[0]
        colarray = (colarray * 255. / (int_level[1] - int_level[0])).astype('uint8')
        level_adjusted_arrays.append(colarray)
    return np.stack(level_adjusted_arrays, axis=-1)

def view_image_slice(volume_directory: pathlib.Path, slider_value: int):
    intensity_arrays = {}

    for imcolor in ['red', 'green', 'blue']:
        resamp_image = sitk.ReadImage(volume_directory / f'resampled_{imcolor}.mhd')
        arr = sitk.GetArrayFromImage(resamp_image).T
        intensity_arrays[imcolor] = arr
        print(intensity_arrays[imcolor].shape)
    
    # save red, green, and blue arrays
    red_arr = intensity_arrays['red']
    green_arr = intensity_arrays['green']
    blue_arr = intensity_arrays['blue']

    # color arrays only done for each slice
    int_arrays = {}
    int_arrays['red'] = red_arr[slider_value, :, :]
    int_arrays['green'] = green_arr[slider_value, :, :]
    int_arrays['blue'] = blue_arr[slider_value, :, :]
    volume = getColorVolume(int_arrays)

    plt.imshow(volume)
    plt.show()

if __name__ == '__main__':
    view_image_slice(pathlib.Path(r'\\allen\programs\mindscope\workgroups\np-behavior\tissuecyte\667252'), 596)