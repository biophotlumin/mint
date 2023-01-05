"""Module containing trajectory processing and filtering functions.
"""
#Imports
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
from scipy import optimize
from scipy.fftpack import *
from scipy.signal import *
import pandas as pd
import trackpy as tp
import math
import cvxpy as cp
from scipy.stats import *

def rejoining(tracks,threshold_t,threshold_r):
    """Rejoins split trajectories.

    Joins trajectories whose start and end points are within set spatial and temporal threshold of each other.  
    
    This function does not generate additional data points, rejoined trajectories will be considered as one for statistical purposes thus reducing oversampling. 

    :param tracks: DataFrame containing trajectories.
    :type tracks: DataFrame
    :param threshold_t: Temporal threshold, in number of frames.
    :type threshold_t: int
    :param threshold_r: Spatial threshold, in pixels.
    :type threshold_r: int
    :return: DataFrame containing rejoined trajectories.
    :rtype: DataFrame
    """    

    df_start = pd.DataFrame()
    df_end = pd.DataFrame()
    n_rejoined = 0
    temp_tracks = tracks.copy()
    #Get first and last point of each trajectory
    for item in set(tracks.particle):

        subtrack = tracks[tracks.particle==item]
        df_temp = subtrack[subtrack.frame==np.min(subtrack.frame)]
        df_start = df_start.append(df_temp,ignore_index=True)

        df_temp = subtrack[subtrack.frame==np.max(subtrack.frame)]
        df_end = df_end.append(df_temp,ignore_index=True)

        df_start = df_start.sort_values(by = 'frame', ascending=False)
        df_end = df_end.sort_values(by = 'frame', ascending=False)

    temp_tracks = tracks.copy()
    

    for linef in df_end.itertuples():
        for lined in df_start.itertuples():
            
            timed = lined.frame
            timef = linef.frame
            particle1 = linef.particle
            particle2 = lined.particle
            
            #Rejoins trajectories if they are within spatial and temporal range
            if (timed > timef) and (timed - timef < threshold_t):
                xd,yd = lined.x,lined.y
                xf,yf = linef.x,linef.y
                r = np.sqrt((xf-xd)**2 + (yf-yd)**2) 
                if r < threshold_r and particle1 != particle2:
                    
                    df_start.loc[df_start['particle']==particle2,'particle']=particle1
                    df_end.loc[df_end['particle']==particle2,'particle']=particle1
                    temp_tracks.loc[temp_tracks['particle']==particle2,'particle']=particle1
                    df_start = df_start.loc[(df_start['frame']!=timed)&(df_start['particle']!=particle1)]
                    df_end = df_end.loc[(df_end['frame']!=timef)&(df_end['particle']!=particle1)]
                    n_rejoined += 1
                    break

    temp_tracks.rename(columns = {'particle':'rejoined_particle'}, inplace = True)     
    tracks = pd.concat([tracks,temp_tracks[['rejoined_particle']]],axis=1,join='inner')

    return tracks, n_rejoined

def SNR_spot_estimation(frames,tracks,base_level):
    """Estimates SNR for each feature.

    Returns a DataFrame with a column containing the SNR and one column containing the signal (as the integral of the 2D gaussian).

    SNR is defined as the height of the gaussian fit divided by the noise from the background + shot noise.

    The signal is defined as : N = 2pi(F-Fo)*sigmaX*sigmaY (volume of 2D gaussian of standard deviation sigmaX and sigmaY in both directions).

    :param frames: 2D array of a given frame.
    :type frames: NumPy array
    :param tracks: DataFrame containing trajectories.
    :type tracks: DataFrame
    :param base_level: Base level of the detector (419 for FNDs experiment).
    :type base_level: int
    :return: DataFrame with added SNR column.
    :rtype: DataFrame
    """    

    nb_frames,nb_rows,nb_columns = frames.shape[0],frames.shape[1],frames.shape[2]
    df = pd.DataFrame()

    for line in tracks.itertuples():
        Pixel_x = int(line.x)
        Pixel_y = int(line.y)
        N = 7 #Half the size in pixel of  the square in which the gaussian is calculated 

        if Pixel_x>7 and Pixel_y>7 and Pixel_x<nb_columns-7 and Pixel_y<nb_rows-7:
            data = frames[line.frame,Pixel_y-N:Pixel_y+N,Pixel_x-N:Pixel_x+N]-base_level
            params = _fit_spot_by_gaussian(data)
            feet,height,width_x,width_y = params[0],params[1][0],params[1][3],params[1][4]
            N = 2*math.pi*height* width_x *width_y #Volume of 2D Gaussian = N photons collected
            if N<0:
                N=0

            nbx = data.shape[0]
            nby = data.shape[1]
            temp = np.copy(data)
            temp[2:nbx-2,2:nby-2] = 0
            squaretemp=0
            for i in range(nbx):
                for j in range(nby):
                    squaretemp += temp[i,j]**2

            if squaretemp/(4*(nbx+nby-4))- (temp.sum()/(4*(nbx+nby-4)))**2<0:
                total_noise = 0
            else:
                BN = np.sqrt(squaretemp/(4*(nbx+nby-4)) - (temp.sum()/(4*(nbx+nby-4)))**2)
                total_noise = np.sqrt(height+BN**2)

            if total_noise != 0:
                SNR = height/total_noise
            else:
                SNR = 0

            df = df.append([{'N':N, 'SNR':SNR,'feet':feet,
                             'particle': line.particle,
                             'frame': line.frame}])
        else:
            df = df.append([{'N':0, 'SNR':0,'feet':0,
                             'particle':line.particle,
                             'frame':line.frame}])

    tracks = tracks.merge(df, on=['particle', 'frame'], how='left')
    return tracks

def acceleration_minimization_norm1(measure, sigma0, px, nn = 0):
    """
    Parameters
    ----------
    measure : array (n, 2)
        measured data (probably noisy) : x and y coordinates
    sigma0 : int
        standard deviation of localisation, in nm
    px : float
        pixel size in µm
    nn : int, optional
        amount of data points not taken into account at the extremities of the solution.
        For some methods, the extreme values are less reliable.

    Returns
    -------
    solution : array (n-2*nn, 2)
        filtered solution with minimization of the norm 1 of the acceleration with difference between measured data and solution inferior or equal to the theoretical noise.
    """
    measure = px*measure
    n = len(measure)       
    variable = cp.Variable((n, 2))
    objective = cp.Minimize(cp.atoms.norm1(variable[2:,0]+variable[:-2,0] - 2*variable[1:-1,0])+cp.atoms.norm1(variable[2:,1]+variable[:-2,1] - 2*variable[1:-1,1]))
    constraints = [ cp.atoms.norm(variable - measure, 'fro')**2 <= n*sigma0**2*10**-6]
    prob = cp.Problem(objective, constraints)
    
    prob.solve(solver='SCS',verbose=False,max_iters=100000) #alternatively, 'GUROBI' or 'MOSEK' #
    solution = variable.value
    if nn == 0:
        return solution
    else:
        return solution[nn:n-nn]

def prec(N):
    """Experimental fit of the precision of localisation.

    :param N: Integrated photon count.
    :type N: float
    """    

    return(30+(600/(np.sqrt(N))))

def acceleration_minimization_norm1_pointwise_adaptative_error(measure, Signal, Noise_function, nn = 0, Solver = 'SCS'):
    """
    Parameters
    ----------
    measure : array (n, 2)
        measured data (probably noisy) : x and y coordinates
    Signal : array(n)
        measured photon count. Higher signal means more photons thus better precision of localisation.
    Noise_function : function
        empirical noise estimation Signal : Noise_function(Signal) = array of the standard deviation of noise.
    nn : int, optional
        amount of data points not taken into account at the extremities of the solution.
        For some methods, the extreme values are less reliable.
    Solver : string, default is 'SCS'
        default solver in cvxpy is ECOS which sometimes fail on complicated problems. SCS is more reliable but slower.
    
    Returns
    -------
    solution : array (n-2*nn, 2)
        filtered solution with 
        - minimization of the norm 1 of the acceleration 
        - constraint : difference between measured data and solution, weighted by Signal, inferior or equal to the esperance of difference between measured data and truth, also weighted by Signal
    """

    n = len(measure)       
    variable = cp.Variable((n, 2))
    objective = cp.Minimize(cp.atoms.norm1(variable[2:,0]+variable[:-2,0] - 2*variable[1:-1,0])+cp.atoms.norm1(variable[2:,1]+variable[:-2,1] - 2*variable[1:-1,1]))
    Weights = np.zeros((n,2))
    Estimated_Noise = Noise_function(Signal)
    Weights[:,0] = 1/Estimated_Noise
    Weights[:,1] = 1/Estimated_Noise
    Constrained = cp.multiply(Weights, variable - measure)

    constraints = [ cp.atoms.norm(Constrained, 'fro')**2 <= 2*n ]
    prob = cp.Problem(objective, constraints)
    
    prob.solve(solver = Solver, max_iters=100000)
    solution = variable.value
    if nn == 0:
        return solution, Estimated_Noise
    else:
        return solution[nn:n-nn]

def minimization(subdata,parameters):
    """Prepares data for minimization.

    :param subdata: DataFrame containing x and y coordinates.
    :type subdata: DataFrame
    :param parameters: Dictionary containing the pixel size under the `'px'` key and the precision of localisation under the `'sigma'` key.
    :type parameters: dict
    :return: DataFrame containing denoised x and y coordinates.
    :rtype: DataFrame
    """    

    #Convert coordinates to µm
    px = parameters['px']
    sigma = parameters['sigma']

    array_x = subdata['x'].to_numpy()
    array_x = array_x[:, np.newaxis]
    
    array_y = subdata['y'].to_numpy()
    array_y = array_y[:, np.newaxis]

    array = np.concatenate((array_x,array_y),axis=1)

    processed_array = acceleration_minimization_norm1(array,sigma,px,nn=0)

    subdata['x'] = processed_array[:,0]
    subdata['y'] = processed_array[:,1]

    #Convert coordinates back to pixels
    subdata['x'] = subdata['x']/px
    subdata['y'] = subdata['y']/px

    subdata.drop(subdata.head(1).index,inplace=True)
    subdata.drop(subdata.tail(1).index,inplace=True)
    subdata = subdata.reset_index(drop = True)

    return subdata

def point_minimization(subdata,parameters):
    """Prepares data for pointwise minimization.

    :param subdata: DataFrame containing x and y coordinates.
    :type subdata: DataFrame
    :param parameters: Dictionary containing the pixel size under the `'px'` key.
    :type parameters: dict
    :return: DataFrame containing denoised x and y coordinates.
    :rtype: DataFrame
    """    

    px = parameters['px']

    #Convert coordinates to µm
    array_x = subdata['x'].to_numpy()
    array_x * px
    array_x = array_x[:, np.newaxis]
    
    array_y = subdata['y'].to_numpy()
    array_y * px
    array_y = array_y[:, np.newaxis]

    #Convert mass to photons
    array_mass = subdata['mass'].to_numpy()
    array_mass = array_mass/(11.4)

    array = np.concatenate((array_x,array_y),axis=1)

    processed_array, estim_noise = acceleration_minimization_norm1_pointwise_adaptative_error(array, array_mass, prec, nn = 0, Solver = 'SCS')

    subdata['x'] = processed_array[:,0]
    subdata['y'] = processed_array[:,1]

    #Convert coordinates back to pixels
    subdata['x'] = subdata['x']/px
    subdata['y'] = subdata['y']/px

    #Include estimated precision of localisation
    subdata['estim_noise'] = estim_noise

    subdata.drop(subdata.head(1).index,inplace=True)
    subdata.drop(subdata.tail(1).index,inplace=True)
    subdata = subdata.reset_index(drop = True)

    return subdata

#Filtering functions

def SNR_threshold(tracks):
    """Filters trajectories based on their Signal to Noise Ratio.

    Returns a DataFrame without trajectories for which the average SNR is below threshold.

    :param tracks: DataFrame containing trajectories.
    :type tracks: DataFrame
    :return: DataFrame containing filtered trajectories. 
    :rtype: DataFrame
    """    

    threshold = 2

    df = pd.DataFrame()
    for item in set(tracks.particle):
        subtracks = tracks[tracks.particle==item]
        average_SNR = np.average(subtracks.SNR)
        if average_SNR >= threshold:
            df = pd.concat([df,subtracks],ignore_index = True)

    return df

def _gaussian(feet,height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters.

    :param feet, height, center_x, center_y, width_x, width_y: Gaussian parameters.
    :type feet, height, center_x, center_y, width_x, width_y: float
    :return: Gaussian function.
    :rtype: func
    """    

    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: feet + height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def _spot_moments(data):
    """Calculates the moments of a 2D gaussian.

    :param data: 2D array of a given frame.
    :type data: NumPy array
    :return: Gaussian parameters. `feet` is calculated by averaging pixels on the edge of 'data'.
    :rtype: float
    """    

    nbx = data.shape[0]
    nby = data.shape[1]
    temp = np.copy(data)
    temp[2:nbx-2,2:nby-2] = 0
    feet = temp.sum()/(4*(nbx+nby-4))

    data = data - feet

    x = int(nbx/2)
    y = int(nby/2)
    width_x = 2.
    width_y = 2.
    height = data.max()
    return feet, height, x, y, width_x, width_y

def _fit_spot_by_gaussian(data):
    """Fits a spot by a 2D gaussian.

    Returns the gaussian fit parameters of a 2D distribution.

    :param data: 2D array of a given frame. 
    :type data: NumPy array
    :return: Gaussian fit parameters.
    :rtype: int
    """    

    parameters = _spot_moments(data)
    feet = parameters[0]
    params = parameters[1:6]
    errorfunction = lambda p: np.ravel(_gaussian(feet,*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(func=errorfunction, x0=params, maxfev=120000) #

    return feet, p

def MSD_filtering(tracks,threshold):
    """Filters trajectories based on their Mean Square Displacement.

    Returns a DataFrame containing trajectories whose calculated MSD is above a set threshold.

    :param tracks: DataFrame containing unfiltered trajectories.
    :type tracks: DataFrame
    :param threshold: MSD threshold.
    :type threshold: int
    :return: DataFrame of filtered trajectories.
    :rtype: DataFrame
    """    

    df = pd.DataFrame()
    for item in set(tracks.particle):
        subtracks = tracks[tracks.particle==item]
        if len(subtracks)<3:
            continue
        df2 = tp.motion.msd(subtracks,1,1,max_lagtime=len(subtracks))
        if max(df2.msd)>threshold:
            df = df.append(subtracks)
    return df

def f(x,a,b,c,d):
    """Third-degree polynom.

    :param x,a,b,c,d: Polynom parameters.
    :type x,a,b,c,d: float
    :return: Polynom.
    :rtype: float
    """
    return a*x**3+b*x**2+c*x+d

def polynomial_fit(data,parameters):
    """Checks wether a trajectory fits a third-degree polynom.

    Calculates how much a trajectory deviates from a third-degree polynom,
    returns `True` if that deviation is below a given threshold.

    :param data: DataFrame containing x and y coordinates.
    :type data: DataFrame
    :param parameters: Dictionary with the threshold under the `'threshold_poly3'` key.
    :type parameters: dict
    :return: `True` if the deviation is below the threshold, `False` if it is greater or equal.
    :rtype: Boolean.
    """    

    #x = data.x
    #y = data.y
    #Fewer trajectories are rejected if they are rotated first    
    rot = rotate_single_track(data)
    x = rot.x_rotated
    y = rot.y_rotated
    n = len(x)
    nn = 10 #Supresses the first and last 10 points
    x = x.iloc[nn:n-nn]
    x = x.reset_index(drop = True)

    y = y.iloc[nn:n-nn]
    y = y.reset_index(drop = True)

    if (len(x)>=parameters['len_cutoff']):
        x=np.array(x)
        y=np.array(y)
        val,cov=optimize.curve_fit(f,x,y)
        deviation=np.sqrt(np.mean((y-f(x,val[0],val[1],val[2],val[3]))**2))
        if (deviation<parameters['threshold_poly3']):
            return True
        else:
            return False

def rotate_single_track(data):
    """Rotates trajectories horizontally.

    :param data: DataFrame containing x and y coordinates.
    :type data: DataFrame
    :return: DataFrame containing rotated x and y coordinates.
    :rtype: DataFrame
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
