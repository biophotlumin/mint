import os

from joblib import Parallel, delayed
from denoising import tophat, wavelet

def p_filtering(frames,settings,parameters):

    
    if settings['tophat']:

        tophat_gen = Parallel(n_jobs=os.cpu_count(),return_as='generator')(delayed(tophat)(parameters['separation'],frame) for frame in frames)

        for i, frame in zip(range(len(frames)), tophat_gen):
            frames[i] = frame

    if settings['wavelet']:

        wavelet_gen = Parallel(n_jobs=os.cpu_count(),return_as='generator')(delayed(wavelet)(frame) for frame in frames)

        for i, frame in zip(range(len(frames)), wavelet_gen):
            frames[i] = frame

    return frames

def phase_calculations_joblib(parameters,data,settings,condition,slide,name,animal):

    f_phase_parameters = pd.DataFrame()

    phase_calc_gen = Parallel(n_jobs=os.cpu_count(),return_as='generator')(delayed(per_phase)\
        (data, trajectory, settings, parameters, condition, slide, name, animal) for trajectory in set(data.particle))
    
    for traj in phase_calc_gen:
        if isinstance(traj,pd.DataFrame):
            f_phase_parameters.reset_index(inplace=True, drop=True)
            f_phase_parameters = pd.concat((f_phase_parameters,traj))
        elif isinstance(traj,str):
            log_analysis['r_len'] += 1
            continue
        else:
            log_analysis['r_poly'] += 1
            continue

    return f_phase_parameters