from scipy.signal import savgol_filter
from traj_calc import *

log_analysis = {
    'r_poly':0,
    'r_speed':0,
}

list_r_conf = []

parameters = {
# Data Extraction
'r_conf_cut' : 0.64,
'px' : 0.173, # in µm
'dt' : 0.05, # in s
'min_thr_prec' : 50, # in nm
'sliding_window':3,
'sigma':129,
'len_cutoff':30, # Number of points
'threshold_poly3':1.4, # Deviation from third-degree polynom
}   
settings = {
# Data Extraction
'polynomial_fit':True,
'minimization':True,
'antero_retro':True,
'theta':True,
'conf_list':False,
}

dt = parameters['dt']
sw = parameters['sliding_window']
conf_threshold = parameters['min_thr_prec']*1e-3/(dt*(sw-1)) # in µm/s, threshold defined by a minimal velocity of 10nm/dt (or 0.2µm/s)

def inst_velocity(x,y,dt):
    """Calculates instantaneous velocity.

    :param x: Array of x coordinates
    :type x: NumPy array
    :param y: Array of y coordinates
    :type y: NumPy array
    :param dt: Sampling period in seconds.
    :type dt: float
    :return: Array of point by point instantaneous velocity
    :rtype: NumPy array
    """    

    size = len(x)
    v_inst = np.zeros(size)
    v_inst[0] = (np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2))/dt # Instantaneous speed of the first point calculated from the first segment
    v_inst[size-1] = (np.sqrt((x[size-1]-x[size-2])**2+(y[size-1]-y[size-2])**2))/dt # Instantaneous speed of the last point calculated from the last segment

    for u in range(size-2):
        u = u+1
        v_inst[u] = (np.sqrt((x[u+1]-x[u-1])**2+(y[u+1]-y[u-1])**2))/(2*dt)

    return v_inst

def confinement(x,y,sw):
    """Calculates point by point confinement ratio

    :param x: Array of x coordinates
    :type x: NumPy array
    :param y: Array of y coordinates
    :type y: NumPy array
    :param sw: Confinement ratio sliding window in number of points.
    :type sw: int
    :return: List of point by point confinement ratio
    :rtype: list
    """    

    r_conf = []
    size=len(x)

    if (size>=2*sw+1):
        for i in range(sw): # Extrapolate the confinement ratio at the beginning of the trajectory
            d_net = np.sqrt((x[sw*2]-x[0])**2+(y[sw*2]-y[0])**2)
            d_total = 0
            for u in range(2*sw):
                d_total = d_total + np.sqrt((x[u+1]-x[u])**2+(y[u+1]-y[u])**2)
            r_conf.append((d_net/d_total))

        for i in range (sw,size-sw):
            if (i == size-sw-1): # Extrapolate a number of points equal to the sliding window at the end
                d_total = 0
                d_net = np.sqrt((x[i+sw]-x[i-sw])**2+(y[i+sw]-y[i-sw])**2)
                for j in range (2*sw):
                    d_total = d_total + np.sqrt((x[(i-sw)+j+1]-x[(i-sw)+j])**2+(y[(i-sw)+j+1]-y[(i-sw)+j])**2)
                if (d_total!=0):
                    for p in range(sw+1):
                        r_conf.append((d_net/d_total))
                else:
                    for p in range(sw+1):
                        r_conf.append(0)
            else:

                d_total = 0
                d_net = np.sqrt((x[i+sw]-x[i-sw])**2+(y[i+sw]-y[i-sw])**2)
                for j in range (2*sw):
                    d_total = d_total + np.sqrt((x[(i-sw)+j+1]-x[(i-sw)+j])**2+(y[(i-sw)+j+1]-y[(i-sw)+j])**2)
                if (d_total!=0):
                        r_conf.append((d_net/d_total))
                else:
                        r_conf.append(0)
    else:
        pass
    return r_conf

def per_traj_calc(data, trajectory, settings, parameters, condition, slide, animal, name):

    subdata = data[data.particle==trajectory]
    subdata = subdata.reset_index(drop = True)
    m_data = pd.DataFrame()

    if settings['polynomial_fit']: # Fit trajectories to third-degree polynom, and discard trajectories that deviate too much
        if polynomial_fit(subdata,parameters['len_cutoff'],parameters['threshold_poly3']) == True:
            pass
        else:
            print('poly')
            log_analysis['r_poly'] += 1
            return

    if settings['minimization']: # Experimental trajectory denoising
        subdata = minimization(subdata,parameters['px'],parameters['sigma']) # Defined noise level for every trajectories
        #subdata = point_minimization(subdata,parameters['px']) # Point-by-point calculation of noise level based on signal intensity

    x = subdata.x
    x = x.dropna()
    
    y = subdata.y
    y = y.dropna()
    
    x = x.reset_index(drop = True)
    y = y.reset_index(drop = True)

    size = len(x)
    
    r_conf = confinement(x,y,sw) # Separate trajectory into phases

    if settings['conf_list']:
        list_r_conf.append(r_conf)

    # Switch from pixels to µm   data_dict = {}
    x = x*parameters['px']
    x = x.dropna()
    y = y*parameters['px']
    y = y.dropna()

    v_inst = inst_velocity(x,y,dt) # Get instantaneous velocity for each point

    phase = np.zeros(size)

    # Categorize each phase as either a GO or STOP phase
    if (size>=2*sw+1):
        for i in range(len(r_conf)):
            if (r_conf[i]>parameters['r_conf_cut']):
                    phase[i] = 2 # GO phase

    else:
        for i in range(sw,(size-sw)): # STOP phase refinment
            vel_list = []
            for j in range((i+(-1*(sw//2))),(i+(sw//2))):
                vel = (np.sqrt((x[j+1]-x[j])**2+(y[j+1]-y[j])**2))/dt
                vel_list.append(vel)

            if (phase[i]==2) & (np.mean(vel_list)<conf_threshold):
                phase[i] = 0
    
    subdata['Vinst'] = v_inst
    subdata = subdata.reset_index(drop = True)
    
    diff = []
    cut = []

    for i in range(size-1):
        diff.append(phase[i+1]-phase[i])

    for j in range(size-1):
        if (diff[j]==1 or diff[j]==-1 or diff[j]==2 or diff[j]==-2):
            cut.append(j+1)

    cut.append(len(subdata))

    min_x = x[cut[0]-1]
    max_x = x[(cut[len(cut)-1])-1]

    min_y = y[cut[0]-1]
    max_y = y[(cut[len(cut)-1])-1]

    for phase_number in range(len(cut)-1): # Per phase processing. '-1' : First and last phases are deleted
        
        start = cut[phase_number]
        stop = cut[phase_number+1]
        sub_phase = subdata.loc[start:stop-1]
        intensity = sub_phase.mass

        sub_v_inst = sub_phase.Vinst 
        phase_length = len(sub_phase)
        phase_duration = phase_length*dt

        variance = (np.std(intensity))**2
        mean_intensity = np.mean(intensity)

        if (phase_length==1):
            variance = 0

        if settings['theta']:
            # Calculate theta angle of each particle based on variation of intensity
            # Specific to nanoKTP or similarly behaving nanoparticles

            thetalist = []
            savgol = savgol_filter(intensity,window_length=9,polyorder=3,mode="nearest")

            for n in savgol:
                thetalist.append(np.arcsin(np.sqrt((n-np.min(savgol))/(np.max(savgol)-np.min(savgol))))\
                        if np.max(savgol) != np.min(savgol) else np.nan) # Prevent division by zero

            theta = np.array(thetalist)*180/np.pi
            theta_std = np.std(theta)

        curvilign_velocity = np.abs(np.mean(sub_v_inst))
        vectorial_velocity = np.abs((np.sqrt((x[stop-1]-x[start])**2+(y[stop-1]-y[start])**2))/(dt*phase_length))
        
        if settings['antero_retro']:
            # Check wether trajectory belongs to the right or left eye
            if slide == "oeil_droit":
                sign = 1
            else:
                sign = -1

            # Change the sign of the velocity accordingly
            if ((x[stop-1]-x[start])>0):
                curvilign_velocity = -sign * curvilign_velocity
                vectorial_velocity = -sign * vectorial_velocity
            else:
                curvilign_velocity = sign * curvilign_velocity
                vectorial_velocity = sign * vectorial_velocity
        
        if (phase[start]==0):
            phase_sign = 0
        if (phase[start]==2):
            phase_sign = 2

        run_length = curvilign_velocity*phase_duration
        curv_length = curvilign_velocity*dt

        data_dict = {'trajectory': trajectory, 
                    'phase':phase_sign,
                    'phase_number':phase_number,
                    'phase_length':phase_length,
                    'vectorial_velocity':vectorial_velocity,
                    'curvilign_velocity':curvilign_velocity,
                    'phase_duration':phase_duration,
                    'run_length':run_length,
                    'intensity':mean_intensity,
                    'variance':variance,
                    'condition':condition,
                    'slide':slide,
                    'curv_length':curv_length,
                    'rejoined_trajectory':subdata.rejoined_particle.unique()[0],
                    'animal':animal,
                    'file':name,
                    'min_x':min_x,
                    'max_x':max_x,
                    'min_y':min_y,
                    'max_y':max_y,
                    'n_particles':data.n_particles.unique()[0],
                    }
        
        if settings['theta']:
            temp_dict = {'theta_std':theta_std}
            data_dict = {**data_dict, **temp_dict}
        
        m_data = pd.concat((m_data, pd.DataFrame([data_dict])))

    return m_data


data_path = r'/media/baptiste/Windows/Users/LUMIN10/Documents/Données/video_benchmark_min Results - 20230605_141811/video_benchmark_min/Exp1_20190205_06_kif5a_nKTP/HET/larve3/oeil_droit/190205_nanoKTP_kif5a.lif - Series006.tif/190205_nanoKTP_kif5a.lif - Series006.tif_rejoined.csv'
data = pd.read_csv(data_path,sep='\t')
data.head()

f_phase_parameters = pd.DataFrame()

for trajectory in set(data.particle):
    
    data_dict = per_traj_calc(data, trajectory, settings, parameters, 'HET', 'oeil_droit', 'larve3', '190205_nanoKTP_kif5a.lif - Series006.tif')
    f_phase_parameters.reset_index(inplace=True, drop=True)
    f_phase_parameters = pd.concat((f_phase_parameters,data_dict))
    f_phase_parameters.dropna(axis=0,how='all',inplace=True)

print(f_phase_parameters)