"""Tests `locate` parameters.
"""

#Imports
import os
import nd2
import imageio
import warnings
import numpy as np
import trackpy as tp
import matplotlib
import matplotlib.pyplot as plt

from pathlib import Path
from utils import get_file_list
from denoising import tophat, wavelet

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

    path_list, name_list = get_file_list(input_folder, parameters['extension_in'])

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

        processed_frame = frames[0]

        processed_frame = processed_frame.astype('float64') # Prevent conflict in case of wavelet filtering

        
        if settings['tophat']:
            processed_frame = tophat(parameters['filter_separation'],processed_frame)
            
        if settings['wavelet']:
            processed_frame = wavelet(processed_frame)

        # Show filter results     

        if settings['show_filter']:
            plt.imshow(processed_frame)
            plt.title("Filter test on "+name,fontsize=10)
            plt.xlabel("x (pixels)",fontsize=10)
            plt.ylabel("y (pixels)",fontsize=10)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.show()
            plt.close()
        
        # Locating particles and finding trajectories

        tp.quiet([True]) # Silencing TrackPy messages

        print(f'\tLocating')

        if list not in {type(parameters['diameter']),type(parameters['minmass']),type(parameters['separation'])}:

            raw_coordinates = tp.locate(processed_frame, minmass=parameters['minmass'], diameter=parameters['diameter'], \
                separation=parameters['separation'],preprocess=False, engine='numba')
            
            plt.imshow(processed_frame)
            plt.title("Locate parameters test on "+name,fontsize=10)
            plt.xlabel("x (pixels)",fontsize=10)
            plt.ylabel("y (pixels)",fontsize=10)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            tp.annotate(raw_coordinates,processed_frame,plot_style={'markersize':10},color='red')
            plt.show()
            plt.close()

        else:
            if type(parameters['diameter']) == list:
                diameter_range = np.arange(parameters['diameter'][0],parameters['diameter'][1],parameters['diameter_scaling'])
            else:
                diameter_range = [parameters['diameter']]

            if type(parameters['minmass']) == list:
                minmass_range = np.arange(parameters['minmass'][0],parameters['minmass'][1],parameters['minmass_scaling'])
            else:
                minmass_range = [parameters['minmass']]

            if type(parameters['separation']) == list:
                separation_range = np.arange(parameters['separation'][0],parameters['separation'][1],parameters['separation_scaling'])
            else:
                separation_range = [parameters['separation']]

            output_folder = Path(input_folder).parent.joinpath(f'{name} - Locate parameters test')
            if not Path.is_dir(output_folder):
                os.makedirs(output_folder)

            matplotlib.use('Agg')
            for diameter in diameter_range:
                for minmass in minmass_range:
                    for separation in separation_range:
                        raw_coordinates = tp.locate(processed_frame, minmass=minmass, diameter=diameter, \
                            separation=separation,preprocess=False,engine='numba')
        
                        plt.imshow(processed_frame)
                        plt.title("Locate parameters test on "+name,fontsize=10)
                        plt.xlabel("x (pixels)",fontsize=10)
                        plt.ylabel("y (pixels)",fontsize=10)
                        plt.xticks(fontsize=10)
                        plt.yticks(fontsize=10)
                        tp.annotate(raw_coordinates,processed_frame,plot_style={'markersize':10},color='red')
                        plt.figtext(0.5,0.01,f'Number of particles found :  {len(raw_coordinates)}', ha='center')
                        plt.savefig(Path(output_folder).joinpath(f'{name}_diameter_{diameter}_minmass_{minmass}_separation_{separation}.png'))
                        plt.close()

if __name__ == '__main__':

    parameters = {
        #trackpy.batch
        'diameter':[7,12],
        'minmass':[250,350],
        'separation':[10,14],
        'filter_separation':12,
        'diameter_scaling':2,
        'minmass_scaling':10,
        'separation_scaling':1,
        'extension_in':'tif',
    }

    #Optional image processing

    settings = {
        'tophat':True,
        'wavelet':False,
        'show_filter':False,
    }

    input_folder = r'/media/lumin/SHG_tracking_data/Zebrafish data/video_benchmark_min/exp1/HET/larve4/oeil_droit/190205_nanoKTP_kif5a.lif - Series018.tif'
    test_locate(input_folder,parameters,settings=settings)