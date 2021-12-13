#Imports
import numpy as np
import os
import imageio
import trackpy as tp
import trajectory_calculations as ft
from image_denoising import *
from output_files_creation import *
from pathlib import Path
from utils import *


def tracking(input_folder,parameters,settings,log):
    """File per file localization of particles and trajectory reconstitution. 
        input_folder is the root folder containing all files to be analyzed.
        parameters is a dictionary of calculation parameters, as defined in script.py.
        settings is a dictionary containing boolean values for optional data processing, as defined in script.py.
        log is a dictionary used to log certain values.
        Outputs .csv files, and optionally .txt and image files, as described in script.py
    """
    for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
        for name in files:
            if name.endswith('.tif') == False:  #Check for correct file extension
                continue #Skips to next file if not .tif           

            #Building output folder
            file_path = os.path.join(path, name) #Get file path of current file
            print(file_path)
            output_subfolder = file_path.replace(str(log['root_input_folder']),'')#Separate subfolder structure from root input folder
            log['output_file_path'] = str(log['output_folder']) + output_subfolder #Merge subfolder structure with root output folder

            os.makedirs(log['output_file_path']) #Create output folder for current file
            print(log['output_file_path'])

            #Opening video file
            frames = imageio.volread(file_path)

            #Per frame denoising process
            frames_init = np.zeros(frames.shape)
            for i in range(len(frames)):
                if settings['tophat']:
                    processed_frames = tophat(parameters,frames_init,i,frames) #Tophat denoising
                else:
                    processed_frames = frames

                if settings['wavelet']:
                    processed_frames = wavelet_denoising(processed_frames,i) #Wavelet denoising
            
            #Localizing particles and finding trajectories

            tp.quiet([True]) #Silencing TrackPy messages

            print('Locating '+name)
            raw_coordinates = tp.batch(processed_frames, minmass=parameters['minmass'], diameter=parameters['diameter'], \
                separation=parameters['separation'],preprocess=False,engine='numba',processes=1)

            print('Linking '+name)
            raw_trajectory = tp.link(raw_coordinates, search_range=parameters['search_range'], adaptive_step= \
                parameters['adaptive_step'], adaptive_stop=parameters['adaptive_stop'],memory=parameters['memory'])

            if settings['stub_filtering']:
                print('Stub filtering '+name)
                raw_trajectory = tp.filter_stubs(raw_trajectory,parameters['stub_filtering'])

            #Unlabeling DataFrame index to prevent future conflict with 'frame' column
            raw_trajectory.index.name = None

            #Dumping raw trajectories into csv file
            trajectory_output(log,name,"",raw_trajectory)
            
            #Optional trajectory processing
            if settings['MSD']: 
                print('MSD '+name)
                processed_trajectory = ft.MSD_filtering(raw_trajectory,parameters['threshold'])
                if len(processed_trajectory) == 0: #Check if any trajectories were found. If not, the threshold might be too high.
                    continue
            else:
                processed_trajectory = raw_trajectory

            if settings['rejoining']:
                print('Rejoining '+name)
                processed_trajectory, n_rejoined  = ft.rejoining(processed_trajectory,parameters['threshold_t'],parameters['threshold_r'])
                log['number_rejoined'] += n_rejoined
            else: 
                processed_trajectory['rejoined_particle'] = processed_trajectory['particle']

            if settings['SNR_estimation']:
                processed_trajectory = ft.SNR_spot_estimation(frames,processed_trajectory,parameters['base_level'])

            #Estimating ratio of moving particles
            first_frame=raw_coordinates[raw_coordinates.frame==0]
            n_particles=len(first_frame)
            n_particles = [n_particles]*len(raw_coordinates)    
            n_particles = pd.DataFrame(n_particles,columns=['n_particles'])
            processed_trajectory = pd.concat([processed_trajectory, n_particles],axis=1,join='inner')
            
            #Dumping rejoined trajectories into csv file
            trajectory_output(log,name,"_rejoined",processed_trajectory)

            #Per trajectory data extraction
            if settings['individual_images'] or settings['individual_txt'] or settings['group_image']:
                print("Saving plots and trajectories")
                for item in set(processed_trajectory.particle):
                    sub_trajectory = processed_trajectory[processed_trajectory.particle==item]
                    if settings['individual_images']:
                        image_output(log,name,frames,processed_trajectory,item) #Plots individual trajectory onto the first frame of the video
                    if settings['individual_txt']:
                        trajectory_separation(log,name,item,settings,sub_trajectory) #Dumps individual trajectory into txt file

                #Plots all trajectories onto the first frame of the video
                if settings['group_image']:
                    final_image_ouput(log,name,frames,processed_trajectory)

    #Writes down parameters and settings used for that run into txt files at the root of the output folder
    dict_dump(log,parameters,'parameters')
    dict_dump(log,settings,'settings')
    dict_dump(log,log,'log')