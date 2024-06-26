"""
Main tracking function.
"""

# Imports
import os
import time
import warnings
import trackpy as tp
import pandas as pd
import numpy as np

from pathlib import Path
from joblib import Parallel, delayed

from .input import get_frames
from .denoising import filtering, filtering_p
from .utils import print_pb, get_file_list, Path_type, logger

from .output import (trajectory_output,
                     trajectory_separation,
                     image_output,)

from .traj_calc import (MSD_calculation,
                        MSD_calculation_p,
                        MSD_filtering,
                        rejoining,
                        SNR_spot_estimation,
                        GFP_mask)

def tracking(
        input_folder: Path_type,
        parameters: dict,
        settings: dict,
        log: dict
        ) -> None:

    path_list, name_list = get_file_list(str(input_folder),
                                         f'.{parameters["extension_in"]}')

    for (path, name, j) in zip(path_list, name_list,  # Looping over file path list
                               [j for j in range(len(path_list))]):
        per_file(path, name, j, len(path_list), parameters, settings, log)

    if len(path_list) == 0:
        raise Exception(f"No files found matching extension "
                        f"{parameters['extension_in']}")

    logger.log('msd_mean', np.mean(logger.get('t_msd')), 'add')
    logger.log('msd_std', np.std(logger.get('t_msd')), 'add')
    logger.log('filt_mean', np.mean(logger.get('t_filt')), 'add')
    logger.log('filt_std', np.std(logger.get('t_filt')), 'add')
    logger.log('rej_mean', np.mean(logger.get('t_rej')), 'add')
    logger.log('rej_std', np.std(logger.get('t_rej')), 'add')
    logger.delete('t_msd')
    logger.delete('t_filt')
    logger.delete('t_rej')

    print_pb('\n', j+1, len(path_list))
    print('\n')

def p_tracking(
        input_folder: Path_type,
        parameters: dict,
        settings: dict,
        log: dict
        ) -> None:

    path_list, name_list = get_file_list(str(input_folder),
                                         f'.{parameters["extension_in"]}')

    track_gen = Parallel(
        n_jobs=os.cpu_count(),
        return_as='generator',
        batch_size='auto')(
            delayed(
                per_file)(
                    path,
                    name,
                    j,
                    len(path_list),
                    parameters,
                    settings,
                    log) for (
                        path,
                        name,
                        j) in
                            zip(path_list,
                            name_list,
                            range(len(path_list))))

    for j, _ in enumerate(track_gen):

        print_pb(f'\tProcessed track', j, len(path_list))

    logger.log('msd_mean', np.mean(logger.get('t_msd')), 'add')
    logger.log('msd_std', np.std(logger.get('t_msd')), 'add')
    logger.log('filt_mean', np.mean(logger.get('t_filt')), 'add')
    logger.log('filt_std', np.std(logger.get('t_filt')), 'add')
    logger.log('rej_mean', np.mean(logger.get('t_rej')), 'add')
    logger.log('rej_std', np.std(logger.get('t_rej')), 'add')
    logger.delete('t_msd')
    logger.delete('t_filt')
    logger.delete('t_rej')

    print('\n')

def per_file(
        path: Path_type,
        name: str,
        j: int,
        j_max: int,
        parameters: dict,
        settings: dict,
        log: dict
        ) -> None:

    # Separate subfolder structure from root input folder
    output_subfolder = str(path).replace(str(log['root_input_folder']), '')
    # Merge subfolder structure with root output folder
    out_path = str(log['output_folder']) + output_subfolder

    print_pb(f'\nProcessing {name}', j, j_max)
    os.makedirs(out_path)  # Create output folder for current file

    # Opening video file

    frames = get_frames(path)

    if len(frames) == 0:
        return

    filt_start = time.time()

    # Prevent conflict in case of wavelet filtering
    processed_frames = frames.astype('float64')

    print_pb('\tFiltering', j, j_max)
    if settings['parallel']:
        filt_func = filtering_p
    else:
        filt_func = filtering

    processed_frames = filt_func(processed_frames, settings, parameters)

    filt_time = time.time() - filt_start
    logger.log('t_filt', filt_time, 'append')

    # Localizing particles and finding trajectories

    tp.quiet(suppress=True)  # Silencing TrackPy messages

    print_pb('\tLocating', j, j_max)
    raw_coordinates = tp.batch(processed_frames,
                                minmass=parameters['minmass'],
                                diameter=parameters['diameter'],
                                separation=parameters['separation'],
                                preprocess=False,
                                engine='numba',
                                processes='auto')

    del processed_frames

    print_pb('\tLinking', j, j_max)
    raw_trajectory = tp.link(raw_coordinates,
                                search_range=parameters['search_range'],
                                adaptive_step=parameters['adaptive_step'],
                                adaptive_stop=parameters['adaptive_stop'],
                                memory=parameters['memory'],
                                link_strategy='numba')

    if settings['stub_filtering']:
        print_pb('\tStub filtering', j, j_max)
        raw_trajectory = tp.filter_stubs(raw_trajectory,
                                            parameters['stub_filtering'])

    # Unlabeling DataFrame index to prevent future conflict with 'frame' column
    raw_trajectory.index.name = None

    # Dumping raw trajectories into csv file
    msd_start = time.time()
    if settings['parallel']:
        raw_trajectory = MSD_calculation_p(raw_trajectory,
                                            parameters['px'],
                                            parameters['dt'])
    else:
        raw_trajectory = MSD_calculation(raw_trajectory,
                                            parameters['px'],
                                            parameters['dt'])
    trajectory_output(out_path, name, "", raw_trajectory)
    msd_time = time.time() - msd_start
    logger.log('t_msd', msd_time, 'append')

    # Optional trajectory processing
    if settings['MSD']:
        print_pb('\tMSD filtering', j, j_max)
        processed_trajectory = MSD_filtering(raw_trajectory, parameters['msd'])
        if len(processed_trajectory) == 0: # Check if any trajectories were found.
            # If not, the threshold might be too high.
            warnings.warn('No trajectories retained, '
                            'MSD threshold might be too high')
            return
    else:
        processed_trajectory = raw_trajectory

    if settings['rejoining']:
        print_pb('\tRejoining', j, j_max)
        rej_start = time.time()
        processed_trajectory, n_rejoined = rejoining(processed_trajectory,
                                                        parameters['threshold_t'],
                                                        parameters['threshold_r'])
        rej_time = time.time() - rej_start
        logger.log('t_rej', rej_time, 'append')
        logger.log('n_rejoined', n_rejoined, 'add')

    else:
        processed_trajectory['rejoined_particle'] = processed_trajectory['particle']

    if settings['SNR_estimation']:
        processed_trajectory = SNR_spot_estimation(frames,
                                                    processed_trajectory,
                                                    parameters['base_level'])

    # Estimating ratio of moving particles
    first_frame = raw_coordinates[raw_coordinates.frame == 0]
    n_particles = len(first_frame)
    n_particles = [n_particles]*len(raw_trajectory)
    n_particles = pd.DataFrame(n_particles, columns=['n_particles'])
    processed_trajectory = processed_trajectory.reset_index(drop=True)
    processed_trajectory = pd.concat([processed_trajectory, n_particles],
                                        axis=1, join='inner')

    # Number of static particles
    static = MSD_filtering(raw_trajectory,
                            parameters['msd'],
                            highpass=False)
    if len(static) == 0:
        n_static = 0
    else:
        static = tp.filter_stubs(static, len(frames)//10)
        n_static = static.particle.nunique()

    n_static = [n_static]*len(raw_trajectory)
    n_static = pd.DataFrame(n_static, columns=['n_static'])
    processed_trajectory = processed_trajectory.reset_index(drop=True)
    processed_trajectory = pd.concat([processed_trajectory, n_static],
                                        axis=1, join='inner')

    if settings['gfp']:
        processed_trajectory = GFP_mask(path, name, processed_trajectory)

    # Dumping rejoined trajectories into csv file
    trajectory_output(out_path, name, "_rejoined", processed_trajectory)

    # Per trajectory data extraction
    if (settings['individual_images'] or
        settings['individual_txt'] or
        settings['group_image']):
        print_pb('\tSaving plots and trajectories', j, j_max)
        for item in set(processed_trajectory.particle):
            sub_trajectory = processed_trajectory[processed_trajectory
                                                    .particle == item]
            if settings['individual_images']:
                # Plot individual trajectory onto the first frame of the video
                image_output(out_path,
                                name,
                                frames,
                                processed_trajectory,
                                item)
            if settings['individual_txt']:
                # Dump individual trajectory into .txt file
                trajectory_separation(out_path,
                                        name,
                                        item,
                                        settings,
                                        sub_trajectory)

        # Plot all trajectories onto the first frame of the video
        if settings['group_image']:
            image_output(out_path,
                            name,
                            frames,
                            processed_trajectory,
                            False)

    logger.log('n_traj', processed_trajectory.rejoined_particle.nunique(), 'add')

    del frames

if __name__ == '__main__':

    from pathlib import Path
    import time

    parameters = {
        'extension_in': 'nd2',
        # trackpy.batch
        'diameter': 7,
        'minmass': 40,
        'separation': 12,
        # trackpy.link
        'search_range': 6,
        'memory': 5,
        'adaptive_stop': 5,
        'adaptive_step': 0.9,
        'stub_filtering': 3, # Minimum length of trajectories, in points
        # trackpy.motion.msd
        'msd': 3, # Threshold for MSD filtering
        # SNR estimation
        'base_level': 0, # Base level for SNR estimation
        # Rejoining
        'threshold_t': 10, # Temporal threshold for trajectory rejoining, in frames
        'threshold_r': 40, # Spatial threshold for trajectory rejoining, in pixels
        'px': 0.11, # Pixel size, in µm
        'dt': 0.05, # Time interaval between frames, in s
    }

    settings = {
        'parallel': False,
        # Denoising
        'tophat': True,
        'wavelet': True,
        # Processing
        'stub_filtering': False,
        'MSD': True,
        'rejoining': True,
        'SNR_estimation': False,
        # Outputs
        'individual_images': False,
        'individual_txt': False,
        'group_image': True,
        }

    log = {}

    start = time.time()
    input_folder = Path(r"")
    tracking(input_folder, parameters, settings, log)
    end = time.time()
    duration = end - start
    print('%dh%s' % (int(duration//3600), f'{int((duration%3600)/60):02d}'))