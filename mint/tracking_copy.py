"""Main tracking function.
"""

#Imports
import os
import nd2
import time
import imageio
import warnings
import trackpy as tp
import traj_calc as ft

from joblib import Parallel, delayed
from utils import print_pb, get_file_list
from denoising import *
from output import *

log_tracking = {
        'n_rejoined':0,
        'n_traj':0,
        'n_files':0,
        't_filt':[],
        't_msd':[],
        't_rej':[],
    }

def tracking_gen(settings,parameters,log,path,name):

    log_traj = {
        'n_rejoined':0,
        'n_traj':0,
        'n_files':0,
        't_filt':[],
        't_msd':[],
        't_rej':[],
    }

    output_subfolder = path.replace(str(log['root_input_folder']),'') # Separate subfolder structure from root input folder
    out_path = str(log['output_folder']) + output_subfolder # Merge subfolder structure with root output folder

    # print_pb(f'\nProcessing {name}',j,len(path_list))
    os.makedirs(out_path) # Create output folder for current file

    #Opening video file         
                                                # match parameters['extension_in']:
    if parameters['extension_in'] == 'tif':     #     case 'tif':
        frames = imageio.volread(path)     #         frames = imageio.volread(file_path)
    elif parameters['extension_in'] == 'nd2':   #     case 'nd2':
        frames = nd2.imread(path)          #         frames = nd2.imread(file_path)
    else:
        warnings.warn(f'Extension {parameters["extension_in"]} is not supported')

    filt_start = time.time()
    processed_frames = frames.astype('float64') # Prevent conflict in case of wavelet filtering

    if settings['parallel']:

        processed_frames_gen = Parallel(n_jobs=os.cpu_count(),return_as='generator')(delayed(tophat)(parameters['separation'],frame) for frame in processed_frames)
        
        for i, frame in zip(range(len(processed_frames)), processed_frames_gen):
            processed_frames[i] = frame
    
    else:
        for i in range(len(frames)):
            if settings['tophat']:
                processed_frames[i] = tophat(parameters['separation'],processed_frames[i])
                
            if settings['wavelet']:
                processed_frames[i] = wavelet(processed_frames[i])

    filt_time = time.time() - filt_start
    log_traj['t_filt'].append(filt_time)

    # Localizing particles and finding trajectories

    tp.quiet([True]) # Silencing TrackPy messages

    # print_pb(f'\tLocating',j,len(path_list))
    raw_coordinates = tp.batch(processed_frames, minmass=parameters['minmass'], diameter=parameters['diameter'], \
        separation=parameters['separation'],preprocess=False,engine='numba',processes='auto')
    
    del processed_frames

    # print_pb(f'\tLinking',j,len(path_list))
    raw_trajectory = tp.link(raw_coordinates, search_range=parameters['search_range'], adaptive_step= \
        parameters['adaptive_step'], adaptive_stop=parameters['adaptive_stop'],memory=parameters['memory'])

    if settings['stub_filtering']:
        # print_pb(f'\tStub filtering',j,len(path_list))
        raw_trajectory = tp.filter_stubs(raw_trajectory,parameters['stub_filtering'])

    # Unlabeling DataFrame index to prevent future conflict with 'frame' column
    raw_trajectory.index.name = None

    # Dumping raw trajectories into csv file
    msd_start = time.time()
    if settings['parallel']:
        raw_trajectory = ft.MSD_calculation_joblib(raw_trajectory,parameters['px'],parameters['dt'])
    else:
        raw_trajectory = ft.MSD_calculation(raw_trajectory,parameters['px'],parameters['dt'])
    trajectory_output(out_path,name,"",raw_trajectory)
    msd_time = time.time() - msd_start
    log_traj['t_msd'].append(msd_time)
    
    # Optional trajectory processing
    if settings['MSD']: 
        # print_pb(f'\tMSD filtering',j,len(path_list))
        processed_trajectory = ft.MSD_filtering(raw_trajectory,parameters['msd'])
        if len(processed_trajectory) == 0: # Check if any trajectories were found. If not, the threshold might be too high.
            warnings.warn('No trajectories retained, MSD threshold might be too high')
            return
    else:
        processed_trajectory = raw_trajectory

    if settings['rejoining']:
        # print_pb(f'\tRejoining',j,len(path_list))
        rej_start = time.time()
        if settings['parallel']:
            processed_trajectory, n_rejoined  = ft.rejoining(processed_trajectory,parameters['threshold_t'],parameters['threshold_r'])
        else:
            processed_trajectory, n_rejoined  = ft.rejoining(processed_trajectory,parameters['threshold_t'],parameters['threshold_r'])
        rej_time = time.time() - rej_start
        log_traj['t_rej'].append(rej_time)
        log_traj['n_rejoined'] += n_rejoined
    else: 
        processed_trajectory['rejoined_particle'] = processed_trajectory['particle']

    if settings['SNR_estimation']:
        processed_trajectory = ft.SNR_spot_estimation(frames,processed_trajectory,parameters['base_level'])

    # Estimating ratio of moving particles
    first_frame = raw_coordinates[raw_coordinates.frame==0]
    n_particles = len(first_frame)
    n_particles = [n_particles]*len(raw_trajectory)    
    n_particles = pd.DataFrame(n_particles,columns=['n_particles'])
    processed_trajectory = processed_trajectory.reset_index(drop=True)
    processed_trajectory = pd.concat([processed_trajectory, n_particles],axis=1,join='inner')

    # Number of static particles
    static = ft.lowpass_MSD_filtering(raw_trajectory,parameters['px'],parameters['dt'],9)
    static = tp.filter_stubs(static,len(frames)//10)
    n_static = static.particle.nunique()
    n_static = [n_static]*len(raw_trajectory)    
    n_static = pd.DataFrame(n_static,columns=['n_static'])
    processed_trajectory = processed_trajectory.reset_index(drop=True)
    processed_trajectory = pd.concat([processed_trajectory, n_static],axis=1,join='inner')
    
    # Dumping rejoined trajectories into csv file
    trajectory_output(out_path,name,"_rejoined",processed_trajectory)

    # Per trajectory data extraction
    if settings['individual_images'] or settings['individual_txt'] or settings['group_image']:
        # print_pb(f'\tSaving plots and trajectories',j,len(path_list))
        for item in set(processed_trajectory.particle):
            sub_trajectory = processed_trajectory[processed_trajectory.particle==item]
            if settings['individual_images']:
                image_output(out_path,name,frames,processed_trajectory,item) # Plot individual trajectory onto the first frame of the video
            if settings['individual_txt']:
                trajectory_separation(out_path,name,item,settings,sub_trajectory) # Dump individual trajectory into .txt file

        # Plot all trajectories onto the first frame of the video
        if settings['group_image']:
            image_output(out_path,name,frames,processed_trajectory,False)

    log_traj['n_traj'] += processed_trajectory.rejoined_particle.nunique()

    del frames

    # return name, log_traj

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

    path_list, name_list = get_file_list(str(input_folder), parameters['extension_in'])

    track_gen = Parallel(n_jobs=int(20),return_as='generator', backend='loky', batch_size='auto')(delayed(tracking_gen)(settings,parameters,log,path,name,) \
        for (path, name) in zip(path_list,name_list))
    
    j = 1

    for _ in track_gen:

        # log_tracking['n_rejoined'] += log_traj['n_rejoined']
        # log_tracking['n_traj'] += log_traj['n_traj']
        # log_tracking['n_files'] += log_traj['n_files']
        # log_tracking['t_filt'].append(log_traj['t_filt'])
        # log_tracking['t_msd'].append(log_traj['t_msd'])
        # log_tracking['t_rej'].append(log_traj['t_rej'])

        print_pb(f'\tProcessed track',j,len(path_list))

        j += 1

    log_tracking['msd_mean'] = np.mean(log_tracking['t_msd'])
    log_tracking['msd_std'] = np.std(log_tracking['t_msd'])
    log_tracking['filt_mean'] = np.mean(log_tracking['t_filt'])
    log_tracking['filt_std'] = np.std(log_tracking['t_filt'])
    log_tracking['rej_mean'] = np.mean(log_tracking['t_rej'])
    log_tracking['rej_std'] = np.std(log_tracking['t_rej'])

    dict_dump(log['output_folder'],log_tracking,'log')                 
    print('\n')