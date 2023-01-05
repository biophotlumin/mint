"""Tests `locate` parameters.
"""

#Imports
import numpy as np
import os
import trackpy as tp
from image_denoising import *
from output_files_creation import *
import matplotlib.pyplot as plt
import imageio

def test_locate(input_folder,parameters,settings):
    """Tests `locate` parameters.

    Displays results of trackpy.locate.

    :param input_folder: Path to input folder containing raw videos.
    :type input_folder: Path or string
    :param parameters: Dictionary containing calculation parameters.
    :type parameters: dict
    :param settings: Dictionary containing calculation settings.
    :type settings: dict
    """    

    for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
        for name in files:
            if name.endswith('.tif') == False:  #Check for correct file extension
                continue #Skips to next file if not .tif           

            #Building file path
            file_path = os.path.join(path, name)
            print(path)
            print(file_path)

            #Opening video file
            frames = imageio.volread(file_path)

            length = 1

            #Initializing frames array
            # frames_init = np.zeros((length,frames.shape[1],frames.shape[2]))
            frames_init = np.zeros(frames.shape)


            filtering_dict = settings['filtering']
            processed_frames = frames

            # for i in range(len(frames)):
            #     for k,v in filtering_dict.items():
            #         processed_frames = eval(filtering_dict[k](parameters,processed_frames,i))


            #Per frame denoising process
            for i in range(length):
                if settings['tophat']:
                    processed_frames = tophat(parameters,frames_init,i,frames) #Tophat denoising
                else:
                    processed_frames = frames

                if settings['wavelet']:
                    processed_frames = wavelet_denoising(processed_frames,i)
            
            plt.imshow(processed_frames[0])
            plt.show()
                 

            #Localizing particles and finding trajectories

            tp.quiet([True]) #Silencing TrackPy messages
            
            raw_coordinates = tp.batch(processed_frames, minmass=parameters['minmass'], diameter=parameters['diameter'], \
                separation=parameters['separation'],preprocess=False,engine='numba',processes='auto')
            
            plt.imshow(processed_frames[i])
            plt.title("Locate parameters test on "+name,fontsize=10)
            plt.xlabel("x (pixels)",fontsize=10)
            plt.ylabel("y (pixels)",fontsize=10)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            tp.annotate(raw_coordinates,frames[i],plot_style={'markersize':10},color='red')
            plt.show()
            plt.close()


parameters = {
    #trackpy.batch
    'diameter':9,
    'minmass':300,
    'separation':12,
}

#Optional image processing

settings = {
    'tophat':True,
    'wavelet':True,
    'filtering':{'tophat':True,'wavelet_denoising':True}
}
if __name__ == '__main__':
    input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/test'
    test_locate(input_folder,parameters,settings=settings)