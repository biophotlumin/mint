"""Tests `locate` parameters.
"""

#Imports
import os
import nd2
import imageio
import warnings
import trackpy as tp
import matplotlib.pyplot as plt

from output import *
from denoising import *

plt.switch_backend('TkAgg')

def test_locate(input_folder,parameters,settings):
    """Test filters and locate parameters.

    :param input_folder: Path to input folder containing raw videos.
    :type input_folder: Path or string
    :param parameters: Dictionary containing calculation parameters.
    :type parameters: dict
    :param settings: Dictionary containing calculation settings.
    :type settings: dict
    """    

    path_list = []
    name_list = []

    for path, subfolder, files in os.walk(input_folder): # Scan entire folder structure for files
            for name in files:
                if name.endswith(f'.{parameters["extension_in"]}') == False:  # Check for correct file extension
                    continue # Skip to next file if not correct extension        

                file_path = os.path.join(path, name) # Get file path of current file
                path_list.append(file_path) # Append to file path list
                name_list.append(name)

    for (path, name, j) in zip(path_list,name_list,[j for j in range(len(path_list))]): # Looping over file path list

        print(f'\nProcessing {name}')

        #Opening video file         
                                                    # match parameters['extension_in']:
        if parameters['extension_in'] == 'tif':     #     case 'tif':
            frames = imageio.volread(path)     #         frames = imageio.volread(file_path)
        elif parameters['extension_in'] == 'nd2':   #     case 'nd2':
            frames = nd2.imread(path)          #         frames = nd2.imread(file_path)
        else:
            warnings.warn(f'Extension {parameters["extension_in"]} is not supported')

        processed_frames = frames[0]

        processed_frames = frames.astype('float64') # Prevent conflict in case of wavelet filtering

        for i in range(len(frames)):
            if settings['tophat']:
                processed_frames[i] = tophat(parameters['separation'],processed_frames[i])
                
            if settings['wavelet']:
                processed_frames[i] = wavelet(processed_frames[i])

        # plt.imshow(processed_frames[0])
        # plt.title("Filter test on "+name,fontsize=10)
        # plt.xlabel("x (pixels)",fontsize=10)
        # plt.ylabel("y (pixels)",fontsize=10)
        # plt.xticks(fontsize=10)
        # plt.yticks(fontsize=10)
        # plt.show()
        # plt.close()
        
        # Locating particles and finding trajectories

        tp.quiet([True]) # Silencing TrackPy messages

        print(f'\tLocating')
        raw_coordinates = tp.batch(processed_frames, minmass=parameters['minmass'], diameter=parameters['diameter'], \
            separation=parameters['separation'],preprocess=False,engine='numba',processes='auto')
        
        plt.imshow(processed_frames[0])
        plt.title("Locate parameters test on "+name,fontsize=10)
        plt.xlabel("x (pixels)",fontsize=10)
        plt.ylabel("y (pixels)",fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        tp.annotate(raw_coordinates,processed_frames[0],plot_style={'markersize':10},color='red')
        plt.show()
        plt.close()

if __name__ == '__main__':

    parameters = {
        #trackpy.batch
        'diameter':9,
        'minmass':150,
        'separation':10,
        'extension_in':'tif',
    }

    #Optional image processing

    settings = {
        'tophat':True,
        'wavelet':True,
    }

    input_folder = r'/media/baptiste/Windows/Users/LUMIN10/Documents/TDmin'
    test_locate(input_folder,parameters,settings=settings)