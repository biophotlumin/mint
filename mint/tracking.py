"""Main tracking function.
"""

#Imports
import os
import nd2
import imageio
import warnings
import trackpy as tp
import traj_calc as ft

from utils import print_pb
from denoising import *
from output import *

log_tracking = {
        'n_rejoined':0,
        'n_traj':0,
        'n_files':0,
    }

def tracking(input_folder,parameters,settings,log):
    """Locates features and link them into trajectories.

    :param input_folder: Path to input folder containing raw videos.
    :type input_folder: Path or string
    :param parameters: Dictionary containing calculation parameters.
    :type parameters: dict
    :param settings: Dictionary containing calculation settings.
    :type settings: dict
    :param log: Dictionary used for logging purposes.
    :type log: dict
    """    

    path_list = []
    name_list = []

    for path, subfolder, files in os.walk(input_folder): # Scan entire folder structure for files
            for name in files:
                if name.endswith(f'.{parameters["extension_in"]}') == False:  # Check for correct file extension
                    continue # Skip to next file if not correct extension  
                log_tracking['n_files'] += 1        

                file_path = os.path.join(path, name) # Get file path of current file
                path_list.append(file_path) # Append to file path list
                name_list.append(name)

    for (path, name, j) in zip(path_list,name_list,[j for j in range(len(path_list))]): # Looping over file path list

        output_subfolder = path.replace(str(log['root_input_folder']),'') # Separate subfolder structure from root input folder
        out_path = str(log['output_folder']) + output_subfolder # Merge subfolder structure with root output folder

        print_pb(f'\nProcessing {name}',j,len(path_list))
        os.makedirs(out_path) # Create output folder for current file

        #Opening video file         
                                                    # match parameters['extension_in']:
        if parameters['extension_in'] == 'tif':     #     case 'tif':
            frames = imageio.volread(path)     #         frames = imageio.volread(file_path)
        elif parameters['extension_in'] == 'nd2':   #     case 'nd2':
            frames = nd2.imread(path)          #         frames = nd2.imread(file_path)
        else:
            warnings.warn(f'Extension {parameters["extension_in"]} is not supported')

        processed_frames = frames.astype('float64') # Prevent conflict in case of wavelet filtering

        for i in range(len(frames)):
            if settings['tophat']:
                processed_frames[i] = tophat(parameters['separation'],processed_frames[i])
                
            if settings['wavelet']:
                processed_frames[i] = wavelet(processed_frames[i])
        
        # Localizing particles and finding trajectories

        tp.quiet([True]) # Silencing TrackPy messages

        print_pb(f'\tLocating',j,len(path_list))
        raw_coordinates = tp.batch(processed_frames, minmass=parameters['minmass'], diameter=parameters['diameter'], \
            separation=parameters['separation'],preprocess=False,engine='numba',processes='auto')

        print_pb(f'\tLinking',j,len(path_list))
        raw_trajectory = tp.link(raw_coordinates, search_range=parameters['search_range'], adaptive_step= \
            parameters['adaptive_step'], adaptive_stop=parameters['adaptive_stop'],memory=parameters['memory'])

        if settings['stub_filtering']:
            print_pb(f'\tStub filtering',j,len(path_list))
            raw_trajectory = tp.filter_stubs(raw_trajectory,parameters['stub_filtering'])

        # Unlabeling DataFrame index to prevent future conflict with 'frame' column
        raw_trajectory.index.name = None

        # Dumping raw trajectories into csv file
        trajectory_output(out_path,name,"",raw_trajectory)
        
        # Optional trajectory processing
        if settings['MSD']: 
            print_pb(f'\tMSD filtering',j,len(path_list))
            processed_trajectory = ft.MSD_filtering(raw_trajectory,parameters['px'],parameters['dt'],parameters['msd'])
            if len(processed_trajectory) == 0: # Check if any trajectories were found. If not, the threshold might be too high.
                warnings.warn('No trajectories retained, MSD threshold might be too high')
                continue
        else:
            processed_trajectory = raw_trajectory

        if settings['rejoining']:
            print_pb(f'\tRejoining',j,len(path_list))
            processed_trajectory, n_rejoined  = ft.rejoining(processed_trajectory,parameters['threshold_t'],parameters['threshold_r'])
            log_tracking['n_rejoined'] += n_rejoined
        else: 
            processed_trajectory['rejoined_particle'] = processed_trajectory['particle']

        if settings['SNR_estimation']:
            processed_trajectory = ft.SNR_spot_estimation(frames,processed_trajectory,parameters['base_level'])

        # Estimating ratio of moving particles
        first_frame = raw_coordinates[raw_coordinates.frame==0]
        n_particles = len(first_frame)
        n_particles = [n_particles]*len(raw_coordinates)    
        n_particles = pd.DataFrame(n_particles,columns=['n_particles'])
        processed_trajectory = pd.concat([processed_trajectory, n_particles],axis=1,join='inner')
        
        # Dumping rejoined trajectories into csv file
        trajectory_output(out_path,name,"_rejoined",processed_trajectory)

        # Per trajectory data extraction
        if settings['individual_images'] or settings['individual_txt'] or settings['group_image']:
            print_pb(f'\tSaving plots and trajectories',j,len(path_list))
            for item in set(processed_trajectory.particle):
                sub_trajectory = processed_trajectory[processed_trajectory.particle==item]
                if settings['individual_images']:
                    image_output(out_path,name,frames,processed_trajectory,item) # Plot individual trajectory onto the first frame of the video
                if settings['individual_txt']:
                    trajectory_separation(out_path,name,item,settings,sub_trajectory) # Dump individual trajectory into .txt file

            # Plot all trajectories onto the first frame of the video
            if settings['group_image']:
                image_output(out_path,name,frames,processed_trajectory,False)

        log_tracking['n_traj'] += processed_trajectory.rejoined_particle.nunique()
        
    dict_dump(log['output_folder'],log_tracking,'log')                 
    print_pb('\n',j+1,len(path_list))
    print('\n')