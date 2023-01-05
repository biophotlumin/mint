#Imports
import numpy as np
import os
import trackpy as tp
from image_denoising import *
from output_files_creation import *
import matplotlib.pyplot as plt
import imageio


def test_locate(input_folder,parameters,whole_file,settings):
    """Displays the results of trackpy.batch on the first file found in input_folder.
        Used to test trackpy.batch parameters before analysis. 
        input_folder is the root folder containing all files to be analyzed, wether they are further divided into subfolder or not.
        parameters is a dictionary of calculation parameters, as defined in script.py.
        whole_file is boolean, where True will find features on the entire file, but False will only process the first frame.
        settings is a dictionary containing boolean values for optional data processing, as defined in script.py.
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

            #Initializing frames array
            n_frames = frames.shape[0]
            n_rows = frames.shape[1]
            n_columns = frames.shape[2]
            frames_array_init = np.zeros((n_frames,n_rows,n_columns))

            #Per frame denoising process
            if whole_file == True:
                length = n_frames
            else:
                length = 1
            for i in range(length):
                if settings['tophat']==True:
                    processed_frames = tophat(parameters,frames_array_init,i,frames) #Tophat denoising
                else:
                    processed_frames = frames_array_init
                    processed_frames[i] = frames[i]
                if settings['wavelet']==True:
                    processed_frames = wavelet_denoising(processed_frames,i) #Wavelet denoising
            
            #Localizing particles and finding trajectories
            raw_coordinates = tp.batch(processed_frames, minmass=parameters['minmass'], diameter=parameters['diameter'], \
                separation=parameters['separation'],preprocess=False,engine='numba',processes=1)
            
            plt.imshow(processed_frames[i])

            tp.annotate(raw_coordinates,frames[i],plot_style={'markersize':10},color='red')

            plt.title("Locate parameters test on "+name,fontsize=35)
            plt.xlabel("x (pixels)",fontsize=30)
            plt.ylabel("y (pixels)",fontsize=30)
            plt.xticks(fontsize=30)
            plt.yticks(fontsize=30)
            plt.close()
            break

parameters = {
    #trackpy.batch
    'diameter':5,
    'minmass':300,
    'separation':10,
}

#Optional image processing

settings = {
    'tophat':True,
    'wavelet':False,
}

#input_folder = "C:\\Users\\LUMIN10\\Desktop\\Programme Baptiste\\INSERM"
#test_locate(input_folder,parameters,whole_file=False,settings=settings)