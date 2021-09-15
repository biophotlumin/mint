"""Module containing trajectory processing and filtering functions.

    rejoining rejoins trajectories that were split up.
    SNR_spot_estimation estimates SNR for each feature.
    acceleration_minimization_norm1 applies an acceleration minimization algorithm to smooth trajectories with a lot of spatial noise.
    minimization prepares data for acceleration_minimization_norm1.
    SNR_threshold filters trajectories for which the average SNR is below threshold.
    _gaussian returns a gaussian function with the given parameters.
    _spot_moments returns the gaussian parameters of a 2D distribution by calculating its moments.
    _fit_spot_by_gaussian returns the gaussian parameters of a 2D distribution found by a fit.
    MSD_filtering computes Mean Square Displacement and filters trajectories accordingly.
    f returns a third-degree polynom.
    polynomial_fit filters trajectories based on how much they fit a third-degree polynom.
"""
#Imports
import numpy as np
from scipy import optimize
from scipy.fftpack import *
from scipy.signal import *
import pandas as pd
import trackpy as tp
import math
from sklearn import linear_model
import cvxpy as cp

def rejoining(tracks,threshold_t,threshold_r):
    """Rejoins split trajectories that are visually part of the same entity. Inputs a DataFrame and two integers, returns a DataFrame.

        Inputs a DataFrame of trajectories from trackpy.link as 'tracks'. 
        threshold_t is the maximum distance, in pixels, between two trajectories for them to be rejoined.
        threshold_r is  the maximum amout of time, in number of frames, between two trajectories for them to be rejoined.
        Returns a DataFrame with rejoined trajectories.
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
    """Estimates SNR for each feature. Inputs a NumPy array, a DataFrame, and an integer. Returns a DataFrame.

        tracks is a DataFrame containing particle, frame, x, y. 'base_level' is the base level of the detector (419 for FNDs experiment).
        Returns a one-column dataframe containing SNR and one-column dataframe containing the signal=integral of the 2D gaussian.
        SNR is defined as the height of the gaussian fit divided by the noise from the background + shot noise.
        The signal is defined as : N=2pi(F-Fo)*sigmaX*sigmaY (volume of 2D Gaussian of standard deviation sigmaX and sigmaY in both directions)
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

def acceleration_minimization_norm1(measure, sigma0,px, nn = 0):
    """
    Parameters
    ----------
    measure : array (n, 2)
        measured data (probably noisy).
    sigma0 : number
        standard deviation caused by the measure (in nanometer)
    nn : int, optional
        number of data that aren't taken into account at the extremities of the solution.
        For some methods, the extreme values are less reliable
    
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
    
    prob.solve(solver='SCS') #alternatively, 'GUROBI' or 'MOSEK'
    solution = variable.value
    if nn == 0:
        return solution
    else:
        return solution[nn:n-nn]

def minimization(subdata,parameters):
    """Prepares data for acceleration_minimization_norm1. Inputs a DataFrame and a dictionary.

        subdata is a DataFrame containing data from a .csv file being processed.
        parameters is a dictionary of calculation parameters, as defined in script.py.
    """
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

    subdata['x'] = subdata['x']/px
    subdata['y'] = subdata['y']/px

    subdata.drop(subdata.head(1).index,inplace=True)
    subdata.drop(subdata.tail(1).index,inplace=True)
    subdata = subdata.reset_index(drop = True)

    return subdata

#Filtering functions

def SNR_threshold(tracks):
    """Returns a new tracks without trajectories for which the average SNR is below threshold. Inputs and returns a DataFrame.
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
    """Returns a gaussian function with the given parameters. Inputs integers, returns a function.
    """
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: feet + height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def _spot_moments(data):
    """Returns the gaussian parameters of a 2D distribution by calculating its moments. Inputs a NumPy array, returns integers.

        feet is calculated by averaging pixels on the edge of 'data'.
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
    """Returns the gaussian parameters of a 2D distribution found by a fit. Inputs a NumPy array, returns integers. 
    """
    parameters = _spot_moments(data)
    feet = parameters[0]
    params = parameters[1:6]
    errorfunction = lambda p: np.ravel(_gaussian(feet,*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(func=errorfunction, x0=params, maxfev=1200)

    return feet, p

def MSD_filtering(tracks,threshold):
    """Computes Mean Square Displacement and filters trajectories accordingly. Inputs a DataFrame and an integer, returns a DataFrame.
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
    """Third-degree polynom. Inputs and returns floats.
        """
    return a*x**3+b*x**2+c*x+d

def polynomial_fit(data,parameters):
    """Calculates wether or not a trajectory is within a specified range of a given third-degree polynom.

        Inputs a DataFrame with x and y coordinates. Returns a boolean.
        """
    x = data.x
    y = data.y
    n = len(x)
    nn = 10 #Supresses the first and last 10 points
    x = x.iloc[nn:n-nn]
    x = x.reset_index(drop = True)

    y = y.iloc[nn:n-nn]
    y = y.reset_index(drop = True)

    if (len(x)>=parameters['len_cutoff']):

        xm = x[:,np.newaxis]
        ym = y[:,np.newaxis]

        model = linear_model.LinearRegression()
        model.fit(xm, ym)

        x=np.array(x)
        y=np.array(y)
        val,cov=optimize.curve_fit(f,x,y)
        deviation=np.sqrt(np.mean((y-f(x,val[0],val[1],val[2],val[3]))**2))

        if (deviation<parameters['threshold_poly3']):
            return True
        else:
            return False

