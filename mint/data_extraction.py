"""Functions used to extract transport parameters from trajectories."""

#Imports
import os
import time
import abc
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.signal import savgol_filter

from .output import dict_dump
from .traj_calc import minimization, polynomial_fit
from .utils import (csv_sniffer,
                    folder_structure_creation,
                    print_pb,
                    Path_type)

log_analysis = {
    'r_poly': 0,
    'r_len': 0,
    'r_speed': 0,
    't_filt': [],
}

list_r_conf = []

class BaseParameter(abc.ABC):
    """
    Abstract base class for trajectory parameters.

    Parameters
    ----------
    antero_retro : str, optional
        If specified, indicates whether to calculate the parameter in anterograde
        or retrograde phases
        default empty
    go_stop : str, optional
        If specified, indicates whether to calculate the parameter in GO or STOP phases
        default empty

    Attributes
    ----------
    results : list
        List of calculated parameter values for each trajectory
    name : str
        Name of the parameter
    column_name : str
        Name of the column in the output DataFrame
    antero_retro : str
        Indicates whether to calculate the parameter in anterograde or retrograde phases
    go_stop : str
        Indicates whether to calculate the parameter in GO or STOP phases

    :meta private:

    """
    def __init__(self, antero_retro: str='', go_stop: str=''):
        self.results = []
        self.name = 'Base parameter'
        self.column_name = 'base'

        self.antero_retro = antero_retro
        if self.antero_retro not in ['antero', 'retro']:
            self.antero_retro = ''

        self.go_stop = go_stop
        if self.go_stop not in ['GO', 'STOP']:
            self.go_stop = ''

    def append_to_results(self, result): # float ?
        self.results.append(result)

    def get_last_result(self):
        return self.results[-1]

    @abc.abstractmethod
    def calculate_results(self, data: pd.DataFrame):
        pass

    def get_len(self) -> int:
        return len(self.results)

    def get_df(self, results=None):
        if results is not None:
            return pd.DataFrame(results, columns=[self.column_name])
        else:
            return pd.DataFrame(self.results, columns=[self.column_name])

    def update_results(self, data, **kwargs):
        self.append_to_results(self.calculate_results(data, **kwargs))

class Intensity(BaseParameter):
    """Mean intensity
    """

    def __init__(self, go_stop: str='', **kwargs):
        super().__init__(go_stop=go_stop, **kwargs)
        self.name = (f'Mean intensity in {self.go_stop} phases'
                    if self.go_stop else 'Mean intensity')
        self.column_name = (f'intensity_{self.go_stop}'
                            if self.go_stop else 'intensity')

    def calculate_results(self, data):
        return np.mean(data.intensity)

class Variance(BaseParameter):
    """Mean variance
    """

    def __init__(self, go_stop: str='', **kwargs):
        super().__init__(go_stop=go_stop, **kwargs)
        self.name = (f'Mean variance in {self.go_stop} phases'
                    if self.go_stop else 'Mean variance')
        self.column_name = (f'variance_{self.go_stop}'
                            if self.go_stop else 'variance')

    def calculate_results(self, data):
        return np.mean(data.variance)

class DiagonalLength(BaseParameter):
    """Diagonal length, calculated as the distance between the first and last points.
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.name = 'Diagonal length'
        self.column_name = 'diag_length'

    def calculate_results(self, data):
        min_x = np.min(data.min_x)
        max_x = np.max(data.max_x)

        min_y = np.min(data.min_y)
        max_y = np.max(data.max_y)

        delta_x = max_x-min_x
        delta_y = max_y-min_y

        return np.sqrt(delta_x**2+delta_y**2)

class CurvilignLength(BaseParameter):
    """Curvilign length, calculated as the sum of run lengths.
    If `absolute` is True, retrograde run lengths are added as absolute values,
    reflecting the total travelled length of the trajectory.
    If `absolute` is False, retrograde run lengths are subtracted from the sum,
    reflecting the distance travelled between the first and last point.
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.name = 'Curvilign length'
        self.column_name = 'curv_length'

    def calculate_results(self, data, absolute=True):
        if absolute is False:
            rl = data.run_length
        else:
            rl = np.abs(data.run_length)

        return np.abs(np.sum(rl))

class PausingFrequency(BaseParameter):
    """Pausing frequency, expressed as the number of stops per minute.
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.name = 'Pausing frequency'
        self.column_name = 'pausing_frequency'
        self.n_stop = []

    def calculate_results(self, data):
        trajectory_time = np.sum(data.phase_duration)
        n_stop = len(data[data.phase == 0])
        self.n_stop.append(n_stop)

        return 60*n_stop/trajectory_time

class Duration(BaseParameter):
    """Duration, in seconds.
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.name = 'Duration'
        self.column_name = 'duration'

    def calculate_results(self, data):
        return np.sum(data.phase_duration)

class PausingTime(BaseParameter):
    """Time spent in STOP phases.
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.name = "Pausing time"
        self.column_name = 'pausing_time'

    def calculate_results(self, data):
        if len(data) == 0:
            return 0
        else:
            return np.mean(data.phase_duration)

class FractionMoving(BaseParameter):
    """Fraction of moving particles.
    If `msd` is False, the  total amount of particles is the number of features
    found in the first frame of the video.
    If `msd` is True, the total amount of particles is the number of trajectories
    above and below the MSD threshold.
    This calculation does not take into account files with no retained trajectories.
    This might artificially increase the fraction of moving particles if the amount
    of moving particles is very low.
    """

    def __init__(self, msd: bool=False):
        super().__init__()
        self.name = ('Fraction of moving particles'
                     if msd is False else
                     'Fraction of moving particles (MSD)')
        self.column_name = ('fraction_moving'
                            if msd is False else
                            'fraction_moving_msd')
        self.msd = msd

    def calculate_results(self, data, n_trajectories):
        if self.msd is False:
            n_ref = data['n_particles'].unique()
            n_ref = n_ref[0] + n_trajectories
        else:
            n_ref = data['n_static'].unique()[0]

        fraction = n_trajectories/n_ref
        return fraction

class RunLength(BaseParameter):
    """Mean distance traveled in GO phases.
    """

    def __init__(self, antero_retro='', **kwargs):
        super().__init__(antero_retro=antero_retro, **kwargs)
        self.name = (f'Mean run length in {self.antero_retro} phases'
                    if self.antero_retro else 'Mean run length')
        self.column_name = (f'run_length_{self.antero_retro}'
                            if self.antero_retro else 'run_length')

    def calculate_results(self, data):
        return np.mean(data.run_length)

class Processivity(BaseParameter):
    """Mean duration of GO phases.
    """

    def __init__(self, antero_retro='', **kwargs):
        super().__init__(antero_retro=antero_retro, **kwargs)
        self.name = (f'Mean processivity in {self.antero_retro} phases'
                    if self.antero_retro else 'Mean processivity')
        self.column_name = (f'processivity_{self.antero_retro}'
                            if self.antero_retro else 'processivity')

    def calculate_results(self, data):
        return np.mean(data.phase_duration)

class FractionPaused(BaseParameter):
    """Fraction of time in STOP phases.
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.name = 'Fraction of time paused'
        self.column_name = 'fraction_paused'

    def calculate_results(self, data, data_STOP):
        if len(data_STOP) == 0:
            return 0
        else:
            time_paused = np.sum(data_STOP.phase_duration)
            total_time = np.sum(data.phase_duration)
            return time_paused/total_time

class Directionality(BaseParameter):
    """Fraction of anterograde or retrograde transport.
    """

    def __init__(self, antero_retro='', **kwargs):
        super().__init__(antero_retro=antero_retro, **kwargs)
        self.name = f'Directionality ({self.antero_retro})'
        self.column_name = 'directionality'
        self.phase_directionality = []

    def calculate_results(self, data, data_antero, data_retro):
        if not self.antero_retro:
            warnings.warn('Directionality cannot be calculated on non-directional data')
            return None
        else:
            if self.antero_retro == 'antero':
                direction = data_antero
            elif self.antero_retro == 'retro':
                direction = data_retro

            distance = np.sum(np.abs(data.run_length))
            distance_direction = np.sum(direction.run_length)

            self.phase_directionality.append(len(direction)/len(data))

            return np.abs(distance_direction)/np.abs(distance)

class CurvilignVelocity(BaseParameter):
    """Average curvilign velocity.
    """

    def __init__(self, antero_retro='', **kwargs):
        super().__init__(antero_retro=antero_retro, **kwargs)
        self.name = (f'Curvilign velocity in {self.antero_retro} phases'
                    if self.antero_retro else 'Curvilign velocity')
        self.column_name = (f'curv_velocity_{self.antero_retro}'
                            if self.antero_retro else 'curv_velocity')

    def calculate_results(self, data):
        return np.mean(data.curvilign_velocity)

class Switch(BaseParameter):
    """Reversals of directionality and related measurements.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = 'Switch'
        self.column_name = 'switch'
        self.switch = [] # Overall number of reversals
        self.switch_a_to_r = [] # Anterograde to retrograde reversals
        self.switch_r_to_a = [] # Retrograde to anterograde reversals
        self.switch_normal = [] # Normalized number of reversals
        self.switch_var_STOP = [] # Variance of signal intensity in STOP phases
        self.pausing_time_antero = [] # Pausing time between anterograde phases
        self.pausing_time_retro = [] # Pausing time between retrograde phases
        self.pausing_time_switch = [] # Pausing time between phases of inverse direction

    def check_index(self, p_data, index):
        if p_data.iloc[0].phase == 0:
            index = index[1:]
        if p_data.iloc[-1].phase == 0:
            index = index[:-1]
        return index

    def calculate_results(self, data, duration):
        switch = 0
        switch_a_to_r = 0
        switch_r_to_a = 0
        switch_var_STOP = []
        pausing_time_switch = []
        pausing_time_antero = []
        pausing_time_retro = []

        for p in set(data.trajectory): # Switch back to non-rejoined trajectories
            p_data = data.loc[data.trajectory == p]
            p_data.sort_values('phase_number', inplace=True)
            p_data.reset_index(inplace=True, drop=True)
            zero_phases = p_data.index[p_data['phase'] == 0]
            zero_phases = self.check_index(p_data, zero_phases)

            for phase in zero_phases:
                first_phase = p_data.iloc[phase-1]
                stop_phase = p_data.iloc[phase]
                second_phase = p_data.iloc[phase+1]
                vel_first_phase = first_phase['curvilign_velocity']
                vel_sec_phase = second_phase['curvilign_velocity']

                if vel_first_phase > 0 and vel_sec_phase < 0:
                    switch += 1
                    switch_a_to_r += 1
                    switch_var_STOP.append(stop_phase.variance)
                    pausing_time_switch.append(stop_phase.phase_duration)
                elif vel_first_phase < 0 and vel_sec_phase > 0:
                    switch += 1
                    switch_r_to_a += 1
                    switch_var_STOP.append(stop_phase.variance)
                    pausing_time_switch.append(stop_phase.phase_duration)
                elif vel_first_phase > 0 and vel_sec_phase > 0:
                    pausing_time_antero.append(stop_phase.phase_duration)
                elif vel_first_phase < 0 and vel_sec_phase < 0:
                    pausing_time_retro.append(stop_phase.phase_duration)

        self.switch.append(switch)
        self.switch_a_to_r.append(switch_a_to_r)
        self.switch_r_to_a.append(switch_r_to_a)
        self.switch_normal.append(switch/duration)

        self.switch_var_STOP.append(np.mean(switch_var_STOP)
                                    if len(switch_var_STOP) > 0 else 0)
        self.pausing_time_switch.append(np.mean(pausing_time_switch)
                                        if len(pausing_time_switch) > 0 else 0)
        self.pausing_time_antero.append(np.mean(pausing_time_antero)
                                        if len(pausing_time_antero) > 0 else 0)
        self.pausing_time_retro.append(np.mean(pausing_time_retro)
                                       if len(pausing_time_retro) > 0 else 0)

        return self.switch

class Theta(BaseParameter):
    """EXPERIMENTAL : Standard deviation of the theta angle.
    """

    def __init__(self, go_stop='', **kwargs):
        super().__init__(go_stop=go_stop, **kwargs)
        self.name = (f'Standard deviation in {self.go_stop} phases'
                    if self.go_stop else 'Standard deviation')
        self.column_name = (f'theta_std_{self.go_stop}'
                            if self.go_stop else 'theta_std')
    def calculate_results(self, data):
        return np.mean(data.theta_std)

class GFPMask(BaseParameter):
    """EXPERIMENTAL : Only retains trajectories located within a GFP mask.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = 'GFP mask'
        self.column_name = 'gfp_mask'

    def calculate_results(self, data):
        return data.gfp.unique()

def inst_velocity(x: np.ndarray, y: np.ndarray, dt: float) -> np.ndarray:
    """
    Calculates instantaneous velocity.

    Parameters
    ----------
    x : numpy array
        Array of x coordinates.
    y : numpy array
        Array of y coordinates.
    dt : float
        Sampling period in seconds.

    Returns
    -------
    numpy array
        Array of point by point instantaneous velocity.
    """


    size = len(x)
    v_inst = np.zeros(size)
    # Instantaneous speed of the first point calculated from the first segment
    v_inst[0] = (np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2))/dt
    # Instantaneous speed of the last point calculated from the last segment
    v_inst[size-1] = (np.sqrt((x[size-1]-x[size-2])**2+(y[size-1]-y[size-2])**2))/dt

    for u in range(size-2):
        u = u+1
        v_inst[u] = (np.sqrt((x[u+1]-x[u-1])**2+(y[u+1]-y[u-1])**2))/(2*dt)

    return v_inst

def confinement(x: np.ndarray | pd.Series,
                y: np.ndarray | pd.Series,
                sw: int) -> list:
    """
    Calculate point-by-point confinement ratio.

    Parameters
    ----------
    x : array_like
        Array of x coordinates.
    y : array_like
        Array of y coordinates.
    sw : int
        Confinement ratio sliding window in number of points.

    Returns
    -------
    list
        List of point-by-point confinement ratio.

    """

    r_conf = []
    size = len(x)

    if (size >= 2*sw+1):
        # Extrapolate the confinement ratio at the beginning of the trajectory
        for i in range(sw):
            d_net = np.sqrt((x[sw*2]-x[0])**2+(y[sw*2]-y[0])**2)
            d_total = 0
            for u in range(2*sw):
                d_total = d_total + np.sqrt((x[u+1]-x[u])**2+(y[u+1]-y[u])**2)
            r_conf.append((d_net/d_total))

        for i in range(sw, size-sw):
            # Extrapolate a number of points equal to the sliding window at the end
            if (i == size-sw-1):
                d_total = 0
                d_net = np.sqrt((x[i+sw]-x[i-sw])**2+(y[i+sw]-y[i-sw])**2)
                for j in range(2*sw):
                    d_total = d_total + np.sqrt((x[(i-sw)+j+1]-x[(i-sw)+j])**2+(y[(i-sw)
                                                +j+1]-y[(i-sw)+j])**2)
                if (d_total != 0):
                    for p in range(sw+1):
                        r_conf.append((d_net/d_total))
                else:
                    for p in range(sw+1):
                        r_conf.append(0)
            else:

                d_total = 0
                d_net = np.sqrt((x[i+sw]-x[i-sw])**2+(y[i+sw]-y[i-sw])**2)
                for j in range(2*sw):
                    d_total = d_total + np.sqrt((x[(i-sw)+j+1]-x[(i-sw)+j])**2+(y[(i-sw)
                                                +j+1]-y[(i-sw)+j])**2)
                if (d_total != 0):
                        r_conf.append((d_net/d_total))
                else:
                        r_conf.append(0)
    else:
        pass
    return r_conf

def per_phase(data, trajectory, settings, parameters, condition, slide, name, animal):

    dt = parameters['dt']
    sw = parameters['sliding_window']

    # in µm/s, threshold defined by a minimal velocity of 10nm/dt (or 0.2µm/s)
    conf_threshold = parameters['min_thr_prec']*1e-3/(dt*(sw-1))

    subdata = data[data.particle == trajectory]
    subdata = subdata.reset_index(drop=True)

    traj_data = pd.DataFrame() ## ? Done outside of the loop for single thread

    if settings['polynomial_fit']:
            if polynomial_fit(subdata, parameters['len_cutoff'],
                              parameters['threshold_poly3']) is True:
                pass
            elif polynomial_fit(subdata, parameters['len_cutoff'],
                                parameters['threshold_poly3']) == 'len':
                return ''
            else:
                return
    if 'ratio_df' in subdata.columns:
        if float(subdata.ratio_df.unique()[0]) > 0.95:
            gfp = True
        else:
            gfp = False
    else:
        gfp = False

    if settings['minimization']: # Experimental trajectory denoising
        filt_start = time.time()
        subdata = minimization(subdata, parameters['px'], parameters['sigma'])
        filt_time = time.time() - filt_start
        log_analysis['t_filt'].append(filt_time)

    x = subdata.x
    x = x.dropna()

    y = subdata.y
    y = y.dropna()

    x = x.reset_index(drop=True)
    y = y.reset_index(drop=True)

    size = len(x)

    r_conf = confinement(x, y, sw) # Separate trajectory into phases

    if settings['conf_list']:
        list_r_conf.append(r_conf)

    # Switch from pixels to µm
    x = x*parameters['px']
    x = x.dropna()
    y = y*parameters['px']
    y = y.dropna()

    v_inst = inst_velocity(x, y, dt) # Get instantaneous velocity for each point

    phase = np.zeros(size)

    # Categorize each phase as either a GO or STOP phase
    if (size >= 2*sw+1):
        for i in range(len(r_conf)):
            if (r_conf[i] > parameters['r_conf_cut']):
                    phase[i] = 2 # GO phase

    else:
        for i in range(sw, (size-sw)): # STOP phase refinment
            vel_list = []
            for j in range((i+(-1*(sw//2))), (i+(sw//2))):
                vel = (np.sqrt((x[j+1]-x[j])**2+(y[j+1]-y[j])**2))/dt
                vel_list.append(vel)

            if (phase[i] == 2) & (np.mean(vel_list) < conf_threshold):
                phase[i] = 0

    subdata['Vinst'] = v_inst
    subdata = subdata.reset_index(drop=True)

    diff = []
    cut = []

    for i in range(size-1):
        diff.append(phase[i+1]-phase[i])

    for j in range(size-1):
        if (diff[j] == 1 or diff[j] == -1 or diff[j] == 2 or diff[j] == -2):
            cut.append(j+1)

    cut.append(len(subdata))

    min_x = x[cut[0]-1]
    max_x = x[(cut[len(cut)-1])-1]

    min_y = y[cut[0]-1]
    max_y = y[(cut[len(cut)-1])-1]

    # Per phase processing. '-1' : First and last phases are deleted
    for phase_number in range(len(cut)-1):

        start = cut[phase_number]
        stop = cut[phase_number+1]
        sub_phase = subdata.loc[start:stop-1]
        intensity = sub_phase.mass

        sub_v_inst = sub_phase.Vinst
        phase_length = len(sub_phase)
        phase_duration = phase_length*dt

        variance = (np.std(intensity))**2
        mean_intensity = np.mean(intensity)

        if (phase_length == 1):
            variance = 0

        if settings['theta']:
            # Calculate theta angle of each particle based on variation of intensity
            # Specific to nanoKTP or similarly behaving nanoparticles

            thetalist = []
            savgol = savgol_filter(intensity, window_length=9, polyorder=3,
                                   mode="nearest").astype('float64') # type: ignore

            for n in savgol:
                thetalist.append(np.arcsin(np.sqrt((n-np.min(savgol))/(np.max(savgol)-np.min(savgol))))
                        # Prevent division by zero
                        if np.max(savgol) != np.min(savgol) else np.nan)

            theta = np.array(thetalist)*180/np.pi
            theta_std = np.std(theta)

        curvilign_velocity = np.abs(np.mean(sub_v_inst))
        vectorial_velocity = np.abs((np.sqrt((x[stop-1]-x[start])**2
                                    +(y[stop-1]-y[start])**2))/(dt*phase_length))

        if settings['antero_retro']:
            # # Check wether trajectory belongs to the right or left eye
            # if slide == "oeil_droit":
            #     sign = 1
            # else:
            #     sign = -1

            # # Change the sign of the velocity accordingly
            # if ((x[stop-1]-x[start]) > 0):
            #     curvilign_velocity = -sign * curvilign_velocity
            #     vectorial_velocity = -sign * vectorial_velocity
            # else:
            #     curvilign_velocity = sign * curvilign_velocity
            #     vectorial_velocity = sign * vectorial_velocity
            if ((x[stop-1] - x[start]) > 0):
                curvilign_velocity = 1 * curvilign_velocity
                vectorial_velocity = 1 * vectorial_velocity
            else:
                curvilign_velocity = -1 * curvilign_velocity
                vectorial_velocity = -1 * vectorial_velocity

        if (phase[start] == 0):
            phase_sign = 0
        if (phase[start] == 2):
            phase_sign = 2

        run_length = curvilign_velocity*phase_duration
        curv_length = curvilign_velocity*dt

        data_dict = {'trajectory': trajectory,
                    'phase': phase_sign,
                    'phase_number': phase_number,
                    'phase_length': phase_length,
                    'vectorial_velocity': vectorial_velocity,
                    'curvilign_velocity': curvilign_velocity,
                    'phase_duration': phase_duration,
                    'run_length': run_length,
                    'intensity': mean_intensity,
                    'variance': variance,
                    'condition': condition,
                    'slide': slide,
                    'curv_length': curv_length,
                    'rejoined_trajectory': subdata.rejoined_particle.unique()[0],
                    'animal': animal,
                    'file': name,
                    'min_x': min_x,
                    'max_x': max_x,
                    'min_y': min_y,
                    'max_y': max_y,
                    'n_particles': data.n_particles.unique()[0],
                    'n_static': data.n_static.unique()[0],
                    'gfp': gfp,
                    }

        if settings['theta']:
            temp_dict = {'theta_std': theta_std}
            data_dict = {**data_dict, **temp_dict}

        traj_data = pd.concat((traj_data, pd.DataFrame([data_dict])))

    return traj_data

def phase_calculations(parameters, data, settings, condition, slide, name, animal):

    f_phase_parameters = pd.DataFrame()

    for trajectory in set(data.particle):
        subdata = data[data.particle == trajectory]
        subdata = subdata.reset_index(drop=True)
        traj_data = per_phase(data, trajectory, settings, parameters, condition,
                              slide, name, animal)
        if isinstance(traj_data, pd.DataFrame):
            f_phase_parameters.reset_index(inplace=True, drop=True)
            f_phase_parameters = pd.concat((f_phase_parameters, traj_data))
        elif isinstance(traj_data, str):
            log_analysis['r_len'] += 1
            continue
        else:
            log_analysis['r_poly'] += 1
            continue
    return f_phase_parameters

def p_phase_calculations(parameters, data, settings, condition, slide, name, animal):

    f_phase_parameters = pd.DataFrame()

    phase_calc_gen = Parallel(n_jobs=os.cpu_count(), return_as='generator')(delayed
                            (per_phase)(data, trajectory, settings, parameters,
                                        condition, slide, name, animal)
                                        for trajectory in set(data.particle))

    for traj in phase_calc_gen:
        if isinstance(traj, pd.DataFrame):
            f_phase_parameters.reset_index(inplace=True, drop=True)
            f_phase_parameters = pd.concat((f_phase_parameters, traj))
        elif isinstance(traj, str):
            log_analysis['r_len'] += 1
            continue
        else:
            log_analysis['r_poly'] += 1
            continue

    return f_phase_parameters

def trajectory_calculations(phase_parameters, parameters, settings):
    """Calculates per trajectory variables of interest.

    :param phase_parameters: DataFrame containing per phase variables of interest.
    :type phase_parameters: DataFrame
    :param settings: Dictionary containing calculation settings
    :type parameters: dict
    :return: DataFrame containing per trajectory variables of interest
    :rtype: DataFrame
    """

    trajectory_list = []
    file_list = []

    phase_parameters.sort_values(by=['file'], inplace=True)
    phase_parameters.sort_values(by=['rejoined_trajectory'], inplace=True)

    # Initialize parameters

    intensity_GO = Intensity(go_stop='GO')
    intensity_STOP = Intensity(go_stop='STOP')
    variance_GO = Variance(go_stop='GO')
    variance_STOP = Variance(go_stop='STOP')
    curvilign_velocity = CurvilignVelocity()
    processivity = Processivity()
    run_length = RunLength()
    pausing_frequency = PausingFrequency()
    pausing_time = PausingTime()
    diag_length = DiagonalLength()
    fraction_paused = FractionPaused()
    moving_particles = FractionMoving(msd=False)
    moving_particles_msd = FractionMoving(msd=True)
    duration = Duration()
    curv_length = CurvilignLength()
    pausing_time = PausingTime()
    gfp = GFPMask()

    if settings['antero_retro']:

        curvilign_velocity_antero = CurvilignVelocity(antero_retro='antero')
        curvilign_velocity_retro = CurvilignVelocity(antero_retro='retro')
        processivity_antero = Processivity(antero_retro='antero')
        processivity_retro = Processivity(antero_retro='retro')
        run_length_antero = RunLength(antero_retro='antero')
        run_length_retro = RunLength(antero_retro='retro')
        directionality = Directionality(antero_retro='retro')
        switch = Switch()

    if settings['theta']:

        theta_std_GO = Theta(go_stop='GO')
        theta_std_STOP = Theta(go_stop='STOP')

    condition = []
    animal = []
    slide = []

    print("Per trajectory calculations of :")

    for file in set(phase_parameters.file.unique()):
        print(f'\t{file}')
        file_data = phase_parameters[(phase_parameters.file == file)]
        n_trajectories = file_data.rejoined_trajectory.nunique()

        for trajectory in set(phase_parameters[(phase_parameters.file == file)]
                        .rejoined_trajectory):

            data = phase_parameters[(phase_parameters.file == file) &
                                (phase_parameters.rejoined_trajectory == trajectory)]
            data = data.reset_index(drop=True)

            # STOP phases
            data_STOP = data.loc[data['phase'] == 0, :]
            data_STOP = data_STOP.reset_index(drop=True)

            # GO phases
            data_GO = data.loc[data['phase'] == 2, :]
            data_GO = data_GO.reset_index(drop=True)

            if len(data_GO) == 0: # Check if trajectory contains at least one GO phase
                continue

            if settings['antero_retro']:
                # Antero
                data_GO_antero = data_GO.loc[data_GO['curvilign_velocity'] > 0, :]
                data_GO_antero = data_GO_antero.reset_index(drop=True)
                # Retro
                data_GO_retro = data_GO.loc[data_GO['curvilign_velocity'] < 0, :]
                data_GO_retro = data_GO_retro.reset_index(drop=True)

            trajectory_list.append(trajectory)
            file_list.append(file)

            condition.append(data.loc[data.rejoined_trajectory == trajectory,
                                       'condition'].unique()[0])
            animal.append(data.loc[data.rejoined_trajectory == trajectory,
                                   'animal'].unique()[0])
            slide.append(data.loc[data.rejoined_trajectory == trajectory,
                                  'slide'].unique()[0])

            # Intensity
            intensity_GO.update_results(data=data_GO)
            intensity_STOP.update_results(data=data_STOP)

            # Variance
            variance_GO.update_results(data=data_GO)
            variance_STOP.update_results(data=data_STOP)

            # Diagonal size
            diag_length.update_results(data=data)

            # Curvilign size
            curv_length.update_results(data=data, absolute=True)

            # Pausing frequency
            pausing_frequency.update_results(data=data)

            # Trajectory duration
            duration.update_results(data=data)

            # Pausing time
            pausing_time.update_results(data=data_STOP)

            # Ratio of moving particles
            moving_particles.update_results(data=data,
                                            n_trajectories=n_trajectories)
            moving_particles_msd.update_results(data=data,
                                                n_trajectories=n_trajectories)

            # Fraction of time paused
            fraction_paused.update_results(data=data,
                                           data_STOP=data_STOP)

            if settings['antero_retro']:

                # Curvilign velocity
                curvilign_velocity_antero.update_results(data=data_GO_antero)
                curvilign_velocity_retro.update_results(data=data_GO_retro)

                # Run length
                run_length_antero.update_results(data=data_GO_antero)
                run_length_retro.update_results(data=data_GO_retro)

                # Processivity
                processivity_antero.update_results(data=data_GO_antero)
                processivity_retro.update_results(data=data_GO_retro)

                # Directionality
                directionality.update_results(data=data_GO,
                                              data_antero=data_GO_antero,
                                              data_retro=data_GO_retro)

                # Directionality reversal
                switch.update_results(data=data,
                                      duration=duration.get_last_result())

            else:
                # Curvilign velocity
                curvilign_velocity.update_results(data)

                # Run length
                run_length.update_results(data)

                # Processivity
                processivity.update_results(data)

            if settings['theta']:

                # Theta standard deviation
                theta_std_GO.update_results(data_GO)
                theta_std_STOP.update_results(data_STOP)

            gfp.update_results(data=data)

    data_dict = {'condition': condition,
                'animal': animal,
                'slide': slide,
                'file': file_list,
                'trajectory': trajectory_list,
                'diag_length': diag_length.results,
                'curvilign_length': curv_length.results,
                'duration': duration.results,
                'pausing_time': pausing_time.results,
                'pausing_frequency': pausing_frequency.results,
                'fraction_paused': fraction_paused.results,
                'n_stop': pausing_frequency.n_stop,
                'fraction_moving': moving_particles.results,
                'fraction_moving_msd': moving_particles_msd.results,
                'intensity_go': intensity_GO.results,
                'intensity_stop': intensity_STOP.results,
                'variance_go': variance_GO.results,
                'variance_stop': variance_STOP.results,
                'gfp': gfp.results,
                }

    if settings['antero_retro']:

        temp_dict = {'curv_velocity_antero': curvilign_velocity_antero.results,
                    'curv_velocity_retro': curvilign_velocity_retro.results,
                    'processivity_antero': processivity_antero.results,
                    'processivity_retro': processivity_retro.results,
                    'run_length_antero': run_length_antero.results,
                    'run_length_retro': run_length_retro.results,
                    'directionality': directionality.results,
                    'phase_dir': directionality.phase_directionality,
                    'pausing_time_antero': switch.pausing_time_antero,
                    'pausing_time_retro': switch.pausing_time_retro,
                    'pausing_time_switch': switch.pausing_time_switch,
                    'switch': switch.switch,
                    'switch_a_to_r': switch.switch_a_to_r,
                    'switch_r_to_a': switch.switch_r_to_a,
                    'switch_normal': switch.switch_normal,
                    'switch_var_stop': switch.switch_var_STOP,
                    }

        data_dict = {**data_dict, **temp_dict}

    else:
        temp_dict = {'curvilign_velocity': curvilign_velocity.results,
                     'run_length': run_length.results,
                    'processivity': processivity.results,
                    }
        data_dict = {**data_dict, **temp_dict}

    if settings['theta']:
        temp_dict = {'theta_std_go': theta_std_GO.results,
                    'theta_std_stop': theta_std_STOP.results,
                    }
        data_dict = {**data_dict, **temp_dict}

    trajectory_parameters = pd.DataFrame(data_dict)

    return trajectory_parameters

def data_extraction(input_folder: Path_type, parameters: dict, settings: dict):
    """Extracts transport parameters from trajectories contained in a folder.

    :param input_folder: Path to a folder of .csv files containing trajectories.
    :type input_folder: str or Path
    :param parameters: Dictionary containing calculation parameters
    :type parameters: dict
    :param settings: Dictionary containing calculation settings
    :type settings: dict
    """

    output_folder, identifier = folder_structure_creation(input_folder)[0:2]

    # Guard against os.walk running on an empty folder,
    # if the input folder is placed at the root of a drive
    input_folder = Path(input_folder)
    if len(input_folder.parents) < 3:
        input_folder = input_folder.parent

    if output_folder.exists() is False:
        os.makedirs(output_folder)

    phase_parameters = pd.DataFrame()
    path_list = []

    # Define output file names
    phase_parameters_output = output_folder.joinpath((f"{identifier[11:26]}"
                                                      "_phase_parameters.csv"))
    traj_parameters_output = output_folder.joinpath((f"{identifier[11:26]}"
                                                     "_trajectory_parameters.csv"))

    for path, subfolder, files in os.walk(input_folder):
        for name in files:
            if name.endswith('_rejoined.csv') is False:
                continue

            # Build output file path
            file_path = os.path.join(path, name)
            path_list.append(file_path)

    print("Per phase calculations of :")

    for (path, j) in zip(path_list, [j for j in range(len(path_list))]):

        file_folder_path = os.path.split(Path(path).parent)[0]
        slide_path, slide = os.path.split(file_folder_path)
        animal_path, animal = os.path.split(slide_path)
        condition = os.path.split(animal_path)[1]

        data = pd.read_csv(path, sep=csv_sniffer(path))

        print_pb(f'\t{str(Path(path).name)}', j, len(path_list))

        if settings['parallel']:
            calc_func = phase_calculations
        else:
            calc_func = phase_calculations
        phase_parameters = pd.concat((phase_parameters, calc_func(parameters, data,
                                    settings, condition, slide, str(Path(path).name),
                                    animal)))

    if phase_parameters.empty:
        raise RuntimeError('No trajectories retained during analysis')
    else:
        phase_parameters.to_csv(phase_parameters_output, sep='\t')

    print('\n')

    trajectory_parameters = trajectory_calculations(phase_parameters,
                                                    parameters, settings)

    # Write results to .csv files
    trajectory_parameters.sort_values(by=['file', 'trajectory'], inplace=True)
    trajectory_parameters.to_csv(traj_parameters_output, sep='\t')
    phase_parameters.drop(['min_x', 'max_x', 'min_y', 'max_y'], axis='columns',
                          inplace=True)
    phase_parameters.to_csv(phase_parameters_output, sep='\t')

    log_analysis['n_t_file'] = trajectory_parameters.file.nunique()
    log_analysis['filt_mean'] = (np.mean(log_analysis['t_filt'])
                                if len(log_analysis['t_filt']) > 0 else 0)
    log_analysis['filt_std'] = (np.std(log_analysis['t_filt'])
                                if len(log_analysis['t_filt']) > 0 else 0)

    if settings['conf_list']:
        df_r_conf = pd.DataFrame({'r_conf': list_r_conf})
        df_r_conf.to_csv(Path(output_folder)
                         .joinpath(f'{identifier[11:26]}_Confinement ratio.csv'))

    dict_dump(Path(output_folder).parent, log_analysis, 'log')
    print('\n')


if __name__ == '__main__':
    import time
    from pathlib import Path

    parameters = {
    # Data Extraction
    'r_conf_cut': 0.817,
    'px': 0.11, # in µm
    'dt': 0.05, # in s
    'min_thr_prec': 50, # in nm
    'sliding_window': 3,
    'sigma': 129,
    'len_cutoff': 10, # Number of points
    'threshold_poly3': 1, # Deviation from third-degree polynom
    }

    settings = {
    'parallel': True,
    # Data Extraction
    'polynomial_fit': True,
    'minimization': False,
    'antero_retro': True,
    'theta': False,
    'conf_list': False,
    }

    start = time.time()
    input_folder = Path(r"/media/lumin/DATA/"
                        r"DATA_DEVRIM Results - 20240104_165320/DATA_DEVRIM")
    data_extraction(input_folder, parameters, settings)
    end = time.time()
    duration = end - start
    print('%dh%s' % (int(duration//3600), f'{int((duration%3600)/60):02d}'))
    print(f'{int(duration//3600)}h{int((duration%3600)/60):02d}')