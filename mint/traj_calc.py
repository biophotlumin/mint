"""
Trajectory processing and filtering functions.
"""
#Imports

import os
import math
import imageio
import numpy as np
import pandas as pd

from trackpy.motion import msd
from pathlib import Path
from typing import cast

from scipy import optimize
from joblib import Parallel, delayed

from .utils import Path_type

try:
    import cvxpy as cp
except ImportError:
    CVXPY_INSTALLED = False
else:
    CVXPY_INSTALLED = True
    from cvxpy.atoms.norm import norm
    from cvxpy.atoms.norm1 import norm1

pd.options.mode.chained_assignment = None

def rejoining(
        tracks: pd.DataFrame,
        threshold_t: int,
        threshold_r: int,
        ) -> tuple[pd.DataFrame, int]:
    """
    Rejoins split trajectories.

    Joins trajectories whose start and end points are within set
    spatial and temporal threshold of each other.

    This function does not generate additional data points,
    rejoined trajectories will be considered as one for statistical purposes
    thus reducing oversampling.

    Parameters
    ----------
    tracks : DataFrame
        DataFrame containing trajectories.
    threshold_t : int
        Temporal threshold, in number of frames.
    threshold_r : int
        Spatial threshold, in pixels.

    Returns
    -------
    DataFrame
        DataFrame containing rejoined trajectories.
    """

    df_start = tracks.groupby(by="particle").agg(
        {'frame': 'min'}).sort_values(by='frame', ascending=False)
    df_end = tracks.groupby(by="particle").agg(
        {'frame': 'max'}).sort_values(by='frame', ascending=False)

    tstart = tracks.merge(df_start, on=['frame', 'particle'], how='right')
    tend = tracks.merge(df_end, on=['frame', 'particle'], how='right')

    joined = tend[['x', 'y', 'frame', 'particle']
                  ].merge(tstart[['x', 'y', 'frame', 'particle']], how='cross')

    joined['r_dist'] = np.sqrt((joined.x_x - joined.x_y)**2 +
                               (joined.y_x - joined.y_y)**2)
    joined['t_dist'] = joined.frame_y - joined.frame_x

    filtered = joined.loc[(joined.t_dist > 0) & (joined.t_dist < threshold_t)]
    filtered = filtered.loc[joined.r_dist < threshold_r]
    filtered = filtered.sort_values(by=['t_dist', 'r_dist'])
    filtered = filtered.drop_duplicates('particle_y')

    n_rejoined = len(filtered)

    tracks['rejoined_particle'] = tracks.particle

    merged_df = tracks.merge(filtered[['particle_y', 'particle_x']],
                             left_on='rejoined_particle', right_on='particle_y',
                             how='left')
    merged_df.reset_index(drop=True, inplace=True)
    tracks.reset_index(drop=True, inplace=True)

    tracks['temp'] = merged_df['particle_x'].fillna(tracks['rejoined_particle'])
    tracks = tracks.drop('rejoined_particle', axis=1)
    tracks.rename(columns={'temp': 'rejoined_particle'}, inplace=True)

    tracks.rejoined_particle = tracks.rejoined_particle.astype(int)

    return tracks, n_rejoined

def SNR_spot_estimation(
        frames: np.ndarray,
        tracks: pd.DataFrame,
        base_level: int,
        ) -> pd.DataFrame:
    """
    Estimates SNR for each feature.

    Parameters
    ----------
    frames : array_like
        2D array of a given frame.
    tracks : DataFrame
        DataFrame containing trajectories.
    base_level : int
        Base level of the detector (419 for FND experiments).

    Returns
    -------
    DataFrame
        DataFrame with added SNR column.

    SNR is defined as the height of the gaussian fit divided by the noise from
    the background + shot noise.

    The signal is defined as : N = 2pi(F-Fo)*sigmaX*sigmaY (volume of 2D gaussian
    of standard deviation sigmaX and sigmaY in both directions).
    """

    _nb_frames, nb_rows, nb_columns = frames.shape[0], frames.shape[1], frames.shape[2]
    df = pd.DataFrame()

    for line in tracks.itertuples():

        Pixel_x = cast(int, line.x)
        Pixel_y = cast(int, line.y)
        N = 7 # Half the size in pixel of the square in which the gaussian is calculated

        if (Pixel_x > 7
            and Pixel_y > 7
            and Pixel_x < nb_columns-7
            and Pixel_y < nb_rows-7):
            data = frames[cast(int, line.frame),
                          int(Pixel_y-N):int(Pixel_y+N),
                          int(Pixel_x-N):int(Pixel_x+N - base_level)]
            params = _fit_spot_by_gaussian(data)
            feet, height, width_x, width_y = (params[0], params[1][0],
                                              params[1][3], params[1][4])
            # Volume of 2D Gaussian = N photons collected
            N = 2*math.pi*height* width_x *width_y
            if N < 0:
                N = 0

            nbx = data.shape[0]
            nby = data.shape[1]
            temp = np.copy(data)
            temp[2:nbx-2, 2:nby-2] = 0
            squaretemp = 0
            for i in range(nbx):
                for j in range(nby):
                    squaretemp += temp[i, j]**2

            if squaretemp/(4*(nbx+nby-4))- (temp.sum()/(4*(nbx+nby-4)))**2 < 0:
                total_noise = 0
            else:
                BN = np.sqrt(squaretemp/(4*(nbx+nby-4))
                             - (temp.sum()/(4*(nbx+nby-4)))**2)
                total_noise = np.sqrt(height+BN**2)

            if total_noise != 0:
                SNR = height/total_noise
            else:
                SNR = 0

            df = pd.concat((df, pd.DataFrame([{'N': N, 'SNR': SNR, 'feet': feet,
                             'particle': line.particle,
                             'frame': line.frame}])))

        else:
            df = pd.concat((df, pd.DataFrame([{'N': 0, 'SNR': 0, 'feet': 0,
                             'particle': line.particle,
                             'frame': line.frame}])))

    tracks = tracks.merge(df, on=['particle', 'frame'], how='left')
    return tracks

def acceleration_minimization_norm1(
        measure: np.ndarray,
        sigma0: int,
        px: float,
        nn: int = 0,
        ) -> np.ndarray:
    """
    Experimental noise reduction algorithm.

    Parameters
    ----------
    measure : array (n, 2)
        Measured data as x and y coordinates.
    sigma0 : int
        Standard deviation of localisation in nm.
    px : float
        Pixel size in µm.
    nn : int, optional
        Number of data points not taken into account
        at the extremities of the solution. For some methods,
        the extreme values are less reliable.

    Returns
    -------
    solution : array (n-2*nn, 2)
        Filtered solution with minimization of the norm 1 of the acceleration with
        difference between measured data and solution inferior or equal to the
        theoretical noise.
    """

    measure = px*measure
    n = len(measure)
    if not CVXPY_INSTALLED:
        print("CVXPY is not installed. Skipping minimization.")
        return measure

    variable = cp.Variable((n, 2))
    objective = cp.Minimize(norm1(variable[2:, 0]+variable[:-2, 0]
                                           - 2*variable[1:-1, 0])
                            +norm1(variable[2:, 1]+variable[:-2, 1]
                                            - 2*variable[1:-1, 1]))
    constraints = [norm(variable - measure, 'fro')**2 <= n*sigma0**2*10**-6]
    prob = cp.Problem(objective, constraints) #type: ignore

    prob.solve(solver='MOSEK', verbose=False) # Alternatively, 'GUROBI', 'MOSEK', 'SCS'

    solution = variable.value

    if nn == 0:
        return solution
    else:
        return solution[nn:n-nn]

def minimization(
        subdata: pd.DataFrame,
        px: float,
        sigma: int,
        ) -> pd.DataFrame:
    """
    Handles data for minimization.

    Parameters
    ----------
    subdata : DataFrame
        DataFrame containing x and y coordinates.
    px : float
        Pixel size, in microns.
    sigma : int
        Estimated precision of localisation, in nm.

    Returns
    -------
    subdata : DataFrame
        DataFrame containing denoised x and y coordinates.
    """


    #Convert coordinates to µm

    array_x = subdata['x'].to_numpy()
    array_x = array_x[:, np.newaxis]

    array_y = subdata['y'].to_numpy()
    array_y = array_y[:, np.newaxis]

    array = np.concatenate((array_x, array_y), axis=1)

    processed_array = acceleration_minimization_norm1(array, sigma, px, nn=0)

    subdata['x'] = processed_array[:, 0]
    subdata['y'] = processed_array[:, 1]

    #Convert coordinates back to pixels
    subdata['x'] = subdata['x']/px
    subdata['y'] = subdata['y']/px

    subdata.drop(subdata.head(1).index, inplace=True)
    subdata.drop(subdata.tail(1).index, inplace=True)
    subdata = subdata.reset_index(drop=True)

    return subdata

#Filtering functions

def SNR_threshold(
        tracks: pd.DataFrame,
        ) -> pd.DataFrame:
    """
    Filters trajectories based on their signal to noise ratio.

    Parameters
    ----------
    tracks : DataFrame
        DataFrame containing trajectories.

    Returns
    -------
    DataFrame
        DataFrame containing filtered trajectories.
    """

    threshold = 2 # TODO Implement as an optional filter

    df = pd.DataFrame()
    for item in set(tracks.particle):
        subtracks = tracks[tracks.particle == item]
        average_SNR = np.average(subtracks.SNR)
        if average_SNR >= threshold:
            df = pd.concat([df, subtracks], ignore_index=True)

    return df

def _gaussian(
        feet: float,
        height: float,
        center_x: float,
        center_y: float,
        width_x: float,
        width_y: float,
        ):
    """
    Returns a 2D Gaussian function with the given parameters.

    Parameters
    ----------
    feet, height, center_x, center_y, width_x, width_y : float
        Gaussian parameters.

    Returns
    -------
    gaussian : func
        A 2D Gaussian function.
    """

    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x, y: feet + height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def _spot_moments(
        data: np.ndarray,
        ):
    """
    Calculates the moments of a 2D gaussian.

    Parameters
    ----------
    data : array_like
        2D array of a given frame.

    Returns
    -------
    feet, height, center_x, center_y, width_x, width_y : float
        Gaussian parameters. `feet` is calculated by
        averaging pixels on the edge of `data`.

    """

    nbx = data.shape[0]
    nby = data.shape[1]
    temp = np.copy(data)
    temp[2:nbx-2, 2:nby-2] = 0
    feet = temp.sum()/(4*(nbx+nby-4))

    data = data - feet

    x = int(nbx/2)
    y = int(nby/2)
    width_x = 2.
    width_y = 2.
    height = data.max()
    return feet, height, x, y, width_x, width_y

def _fit_spot_by_gaussian(
        data: np.ndarray,
        ) -> tuple[float, np.ndarray]:
    """
    Fits a spot by a 2D gaussian.

    Parameters
    ----------
    data : array_like
        2D array of a given frame.

    Returns
    -------
    tuple of float
        Gaussian fit parameters (feet, height, center_x, center_y, width_x, width_y).
    """

    parameters = _spot_moments(data)
    feet = parameters[0]
    params = parameters[1:6]

    def errorfunction(p):
        return np.ravel(_gaussian(feet, *p)(*np.indices(data.shape)) - data)

    p, _success = optimize.leastsq(func=errorfunction, x0=params, maxfev=120000) #

    return feet, p

def MSD_per_traj(
        tracks: pd.DataFrame,
        traj: int,
        px: float,
        dt: float,
        ) -> pd.DataFrame:
    """
    Calculate the Mean Square Displacement (MSD) per trajectory.

    Parameters
    ----------
    tracks : pd.DataFrame
        DataFrame containing trajectories.
    traj : int
        The trajectory ID to calculate MSD for.
    px : float
        Pixel size, in microns.
    dt : float
        Sampling period, in seconds.

    Returns
    -------
    pd.DataFrame
        DataFrame with an extra column for the maximum MSD.
    """

    subtracks = tracks[tracks.particle == traj].copy()
    df2 = msd(subtracks, px, (1/dt), max_lagtime=len(subtracks))
    max_msd = [df2.msd.max()]*len(subtracks)
    subtracks['max_msd'] = max_msd

    return subtracks

def MSD_calculation(
        tracks: pd.DataFrame,
        px: float,
        dt: float,
        ) -> pd.DataFrame:
    """
    Calculates the maximum Mean Square Displacement of trajectories.

    Parameters
    ----------
    tracks : DataFrame
        DataFrame containing trajectories.
    px : float
        Pixel size, in microns.
    dt : float
        Sampling period, in seconds.

    Returns
    -------
    DataFrame
        DataFrame with an extra column for maximum MSD.
    """

    df = pd.DataFrame()
    for traj in set(tracks.particle):
        subtracks = MSD_per_traj(tracks, traj, px, dt)
        if len(subtracks) > 1:
            df = pd.concat((df, subtracks))
    return df

def MSD_calculation_p(
        tracks: pd.DataFrame,
        px: float,
        dt: float,
        ) -> pd.DataFrame:
    """
    Calculate the maximum Mean Square Displacement of trajectories.
    Parallelized with `joblib`.

    Parameters
    ----------
    tracks : DataFrame
        DataFrame containing trajectories.
    px : float
        Pixel size in microns.
    dt : float
        Sampling period in seconds.

    Returns
    -------
    DataFrame
        DataFrame with an extra column for maximum MSD.
    """

    df = pd.DataFrame()
    MSD_generator = Parallel(n_jobs=os.cpu_count(),
                             return_as='generator')(delayed(MSD_per_traj)
                                                    (tracks, traj, px, dt)
                                                      for traj in set(tracks.particle))
    for subtracks in MSD_generator:
        subtracks = cast(pd.DataFrame, subtracks) # For typing purposes
        if len(subtracks) > 1:
            df = pd.concat((df, subtracks))

    return df

def MSD_filtering(
        tracks: pd.DataFrame,
        threshold: int,
        highpass: bool = True,
        ) -> pd.DataFrame:
    """
    Filter trajectories based on Mean Square Displacement (MSD).

    Parameters
    ----------
    tracks : DataFrame
        DataFrame containing unfiltered trajectories.
    threshold : int
        MSD threshold.

    Returns
    -------
    DataFrame
        DataFrame of filtered trajectories.
    """

    df = pd.DataFrame()
    for item in set(tracks.particle):
        subtracks = tracks[tracks.particle == item]
        if len(subtracks) < 3: # MSD filtering has a built-in stub filter of 3
            continue
        if highpass is False:
            if subtracks.max_msd.unique()[0] <= threshold:
                df = pd.concat((df, subtracks))
        else:
            if subtracks.max_msd.unique()[0] > threshold:
                df = pd.concat((df, subtracks))

    return df

def f(x, a, b, c, d):
    """
    Calculates a third-degree polynomial.

    Parameters
    ----------
    x : float
        Polynomial variable.
    a, b, c, d : float
        Polynomial coefficients.

    Returns
    -------
    float
        Polynomial value.
    """
    return a*x**3+b*x**2+c*x+d

def polynomial_fit(
        data: pd.DataFrame,
        len_cutoff: int,
        threshold: float,
        ) -> bool | str:
    """
    Checks whether a trajectory fits a third-degree polynom.

    Calculates the deviation of a trajectory from a third-degree polynom,
    and returns a boolean array where True indicates that the deviation is
    below the given threshold.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing x and y coordinates.
    len_cutoff : int
        Minimum amount of points in the trajectory.
        Effectively equivalent to stub filtering.
    threshold : float
        Deviation from the third-degree polynom above which the trajectory
        is rejected.

    Returns
    -------
    bool or str
        Returns `True` if the trajectory fits a third-degree polynom,
        `False` it it doesn't, and 'len' if the trajectory rejected because of length.
    """

    # x = data.x
    # y = data.y

    # Fewer trajectories are rejected if they are rotated first
    rot = rotate_single_track(data)
    x = rot.x_rotated
    y = rot.y_rotated
    n = len(x)
    nn = 10 # Supresses the first and last 10 points
    x = x.iloc[nn:n-nn]
    x = x.reset_index(drop=True)

    y = y.iloc[nn:n-nn]
    y = y.reset_index(drop=True)

    if (len(x) >= len_cutoff):
        x = np.array(x)
        y = np.array(y)
        val, _cov = optimize.curve_fit(f, x, y)
        deviation = np.sqrt(np.mean((y-f(x, val[0], val[1], val[2], val[3]))**2))
        if (deviation < threshold):
            return True
        else:
            return False
    else:
        return 'len'

def rotate_single_track(
        data: pd.DataFrame,
        ) -> pd.DataFrame:
    """
    Rotate trajectories horizontally.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing x and y coordinates.

    Returns
    -------
    DataFrame
        DataFrame containing rotated x and y coordinates.
    """

    coords = data.loc[:, ['x', 'y']].values
    coords = coords - coords[0, :]

    distances = (coords ** 2).sum(axis=1)
    furthest = np.argmax(distances)

    theta = np.arctan2(coords[furthest, 1], coords[furthest, 0])
    rotation = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    coords = np.matmul(coords, rotation)

    return pd.DataFrame({
        'x_rotated': coords[:, 0],
        'y_rotated': coords[:, 1],
    })

## Experimental

def prec(
        N: float,
        ) -> float:
    """
    Experimental fit of the precision of localisation.

    Parameters
    ----------
    N : float
        Integrated photon count.

    Returns
    -------
    float
        Estimated precision of localization.
    """

    return(30+(600/(np.sqrt(N))))

def acc_min_norm1_pointwise_adaptative_error(measure,
                                             Signal,
                                             Noise_function,
                                             nn=0,
                                             Solver='SCS'):
    """
    EXPERIMENTAL

    Parameters
    ----------
    measure : array (n, 2)
        measured data (probably noisy) : x and y coordinates
    Signal : array(n)
        measured photon count. Higher signal means more photons
        thus better precision of localisation.
    Noise_function : function
        empirical noise estimation Signal : Noise_function(Signal) =
        array of the standard deviation of noise.
    nn : int, optional
        amount of data points not taken into account at the extremities of the solution.
        For some methods, the extreme values are less reliable.
    Solver : string, default is 'SCS'
        default solver in cvxpy is ECOS which sometimes fail on complicated problems.
        SCS is more reliable but slower.

    Returns
    -------
    solution : array (n-2*nn, 2)
        filtered solution with
        - minimization of the norm 1 of the acceleration
        - constraint : difference between measured data and solution,
        weighted by Signal, inferior or equal to the esperance of difference
        between measured data and truth, also weighted by Signal
    """

    n = len(measure)
    if not CVXPY_INSTALLED:
        print("CVXPY is not installed. Skipping minimization.")
        return measure

    variable = cp.Variable((n, 2))
    objective = cp.Minimize(norm1(variable[2:, 0] +
                                           variable[:-2, 0] - 2*variable[1:-1, 0]) +
                            norm1(variable[2:, 1] +
                                           variable[:-2, 1] - 2*variable[1:-1, 1]))
    Weights = np.zeros((n, 2))
    Estimated_Noise = Noise_function(Signal)
    Weights[:, 0] = 1/Estimated_Noise
    Weights[:, 1] = 1/Estimated_Noise
    Constrained = cp.multiply(Weights, variable - measure)

    constraints = [norm(Constrained, 'fro')**2 <= 2*n]
    prob = cp.Problem(objective, constraints) # type: ignore

    prob.solve(solver=Solver, max_iters=100000)
    solution = variable.value
    if nn == 0:
        return solution, Estimated_Noise
    else:
        return solution[nn:n-nn]

def point_minimization(
        subdata: pd.DataFrame,
        px: float,
        ) -> pd.DataFrame:
    """
    EXPERIMENTAL

    Prepare data for pointwise minimization.

    Parameters
    ----------
    subdata : DataFrame
        DataFrame containing x and y coordinates.
    px : float
        Pixel size, in microns.

    Returns
    -------
    subdata : DataFrame
        DataFrame containing denoised x and y coordinates.
    """

    #Convert coordinates to µm
    array_x = subdata['x'].to_numpy()
    array_x = array_x * px
    array_x = array_x[:, np.newaxis]

    array_y = subdata['y'].to_numpy()
    array_y = array_y * px
    array_y = array_y[:, np.newaxis]

    #Convert mass to photons
    array_mass = subdata['mass'].to_numpy()
    array_mass = array_mass/(11.4)

    array = np.concatenate((array_x, array_y), axis=1)

    (processed_array,
     estim_noise) = acc_min_norm1_pointwise_adaptative_error(array,
                                                            array_mass,
                                                            prec, nn=0,
                                                            Solver='SCS')

    subdata['x'] = processed_array[:, 0]
    subdata['y'] = processed_array[:, 1]

    #Convert coordinates back to pixels
    subdata['x'] = subdata['x']/px
    subdata['y'] = subdata['y']/px

    #Include estimated precision of localisation
    subdata['estim_noise'] = estim_noise

    subdata.drop(subdata.head(1).index, inplace=True)
    subdata.drop(subdata.tail(1).index, inplace=True)
    subdata = subdata.reset_index(drop=True)

    return subdata

def GFP_mask(
        path: Path_type,
        name: str,
        trajectories: pd.DataFrame,
        ) -> pd.DataFrame:
    """
    EXPERIMENTAL

    Loads a matching GFP image for each movie, and calculates which points of each
    trajectories are inside GFP-positive areas.

    Parameters
    ----------
    path : Path_type
        The path to the GFP image.
    name : str
        The name of the file.
    trajectories : pd.DataFrame
        The DataFrame containing trajectories.

    Returns
    -------
    pd.DataFrame
        The DataFrame containing trajectories with the GFP column added.
    """

    folder = Path(path).parent
    img_name = name[:-7]+'GFP.tif'

    try:
        img = imageio.imread(Path(folder).joinpath(img_name))
        # print('File found !')

    except FileNotFoundError:
        print(f'File not found : {img_name}')
        return trajectories

    threshold = 110
    mask = img
    mask[mask >= threshold] = 4096
    mask[mask < threshold] = 0

    mask_list = []

    for x, y in zip(trajectories.x, trajectories.y):
        if mask[int(y), int(x)] == 4096:
            mask_list.append('in')
        else:
            mask_list.append('out')

    trajectories['GFP_mask'] = mask_list

    ratio_list = []
    streak_list = []

    df = pd.DataFrame()

    for traj in trajectories.rejoined_particle.unique():
        traj = trajectories[trajectories.rejoined_particle == traj]
        ratio = len(traj[traj.GFP_mask == 'in'])/len(traj)
        ratio = [ratio]*len(traj)
        ratio_list.extend(ratio)
        ratio = pd.DataFrame(ratio, columns=['ratio_df'])
        traj['ratio_df'] = ratio.values
        streak = longest_streak(list(traj.GFP_mask))
        streak = [streak]*len(traj)
        streak = pd.DataFrame(streak, columns=['streak_df'])
        streak_list.extend(streak)
        traj['streak_df'] = streak.values
        df = pd.concat((df, traj))

    return df

def longest_streak(
        lst: list,
        ) -> int:
    """
    Calculates the longest streak of consecutive elements in the input list.

    Parameters
    ----------
    lst : list
        A list of elements to find the longest streak in.

    Returns
    -------
    int
        The length of the longest streak of consecutive elements in the input list.
    """

    max_streak = 0
    current_streak = 0
    current_element = None

    for element in lst:
        if element == current_element:
            current_streak += 1
        else:
            current_streak = 1
            current_element = element

        max_streak = max(max_streak, current_streak)

    return max_streak