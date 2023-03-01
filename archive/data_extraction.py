"""Functions used to extract transport parameters from trajectories."""

#Imports
import os
import numpy as np
import pandas as pd

from pathlib import Path
from output import dict_dump
from traj_calc import *
from scipy.signal import savgol_filter
from utils import folder_structure_creation, csv_sniffer

log_analysis = {
    'r_poly':0,
    'r_speed':0,
}

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
    v_inst[0] = (np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2))/dt #Instantaneous speed of the first point calculated from the first segment
    v_inst[size-1] = (np.sqrt((x[size-1]-x[size-2])**2+(y[size-1]-y[size-2])**2))/dt #Instantaneous speed of the last point calculated from the last segment

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
        for i in range(sw): #Extrapolate the confinement ratio at the beginning of the trajectory
            d_net = np.sqrt((x[sw*2]-x[0])**2+(y[sw*2]-y[0])**2)
            d_total = 0
            for u in range(2*sw):
                d_total = d_total + np.sqrt((x[u+1]-x[u])**2+(y[u+1]-y[u])**2)
            r_conf.append((d_net/d_total))

        for i in range (sw,size-sw):
            if (i == size-sw-1): #Extrapolate a number of points equal to the sliding window at the end
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

def phase_calculations(parameters,data,settings,condition,slide,name,animal):
    """Calculates per phase variables of interest.

    :param parameters: Dictionary containing calculation parameters
    :type parameters: dict
    :param data: DataFrame containing trajectory coordinates
    :type data: DataFrame
    :param settings: Dictionary containing calculation settings
    :type settings: dict
    :param condition, slide, name, animal: Elements of the file path of the current .csv. Used for statistical analysis downstream.
    :type condition, slide, name, animal: str
    :return: DataFrame containing variables of interest.
    :rtype: DataFrame
    """    

    dt = parameters['dt']
    sw = parameters['sliding_window']
    conf_threshold = parameters['min_thr_prec']*1e-3/(dt*(sw-1)) # in µm/s, threshold defined by a minimal velocity of 10nm/dt (or 0.2µm/s)
    f_phase_parameters = pd.DataFrame()

    for trajectory in set(data.particle):

        subdata = data[data.particle==trajectory]
        subdata = subdata.reset_index(drop = True)

        if settings['polynomial_fit']: #Fit trajectories to third-degree polynom, and discard trajectories that deviate too much
            if polynomial_fit(subdata,parameters['len_cutoff'],parameters['threshold_poly3']) == True:
                pass
            else:
                log_analysis['r_poly'] += 1
                continue

        if settings['minimization']: #Experimental trajectory denoising
            subdata = minimization(subdata,parameters['px'],parameters['sigma']) #Defined noise level for every trajectories
            #subdata = point_minimization(subdata,parameters['px']) #Point-by-point calculation of noise level based on signal intensity

        x = subdata.x
        x = x.dropna()
        
        y = subdata.y
        y = y.dropna()
        
        x = x.reset_index(drop = True)
        y = y.reset_index(drop = True)

        size = len(x)
        
        r_conf = confinement(x,y,sw) #Separate trajectory into phases

        #Switch from pixels to µm   
        x = x*parameters['px']
        x = x.dropna()
        y = y*parameters['px']
        y = y.dropna()

        v_inst = inst_velocity(x,y,dt) #Get instantaneous velocity for each point

        phase = np.zeros(size)

        #Categorize each phase as either a GO or STOP phase
        if (size>=2*sw+1):
            for i in range(len(r_conf)):
                if (r_conf[i]>parameters['r_conf_cut']):
                        phase[i] = 2 #GO phase

        else:
            for i in range(sw,(size-sw)): #STOP phase refinment
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

        for phase_number in range(len(cut)-1): #Per phase processing. '-1' : First and last phases are deleted
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
                #Calculates theta angle of each particle based on variation of intensity
                #Specific to nanoKTP or similarly behaving nanoparticles

                thetalist = []
                savgol = savgol_filter(intensity,window_length=9,polyorder=3,mode="nearest")

                for n in savgol:
                    thetalist.append(np.arcsin(np.sqrt((n-np.min(savgol))/(np.max(savgol)-np.min(savgol))))\
                         if np.max(savgol) != np.min(savgol) else np.nan) #Prevents division by zero

                theta = np.array(thetalist)*180/np.pi
                theta_std = np.std(theta)

            curvilign_velocity = np.abs(np.mean(sub_v_inst))
            vectorial_velocity = np.abs((np.sqrt((x[stop-1]-x[start])**2+(y[stop-1]-y[start])**2))/(dt*phase_length))
            
            if settings['antero_retro']:
                #Checks wether trajectory belongs to the right or left eye
                if slide == "oeil_droit":
                    sign = 1
                else:
                    sign = -1

                #Change the sign of the velocity accordingly
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

            f_phase_parameters.reset_index(inplace=True, drop=True)
            f_phase_parameters = pd.concat((f_phase_parameters,pd.DataFrame([data_dict])))

    return f_phase_parameters

def trajectory_calculations(phase_parameters,settings):
    """Calculates per trajectory variables of interest.

    :param phase_parameters: DataFrame containing per phase variables of interest. 
    :type phase_parameters: DataFrame
    :param settings: Dictionary containing calculation settings
    :type parameters: dict
    :return: DataFrame containing per trajectory variables of interest
    :rtype: DataFrame
    """

    iteration = []
    file_list=[]

    phase_parameters.sort_values(by=['file'],inplace=True)
    phase_parameters.sort_values(by=['rejoined_trajectory'],inplace=True)

    #Initializes lists

    intensity_GO,intensity_STOP,variance_GO,variance_STOP,curvilign_velocity,processivity,run_length,pausing_frequency,pausing_time,\
        diag_size,fraction_paused,moving_particles,number_stop,duration,curv_length,pausing_time = ([] for _ in range(16))

    if settings['antero_retro']:

        curvilign_velocity_antero,curvilign_velocity_retro,processivity_antero,processivity_retro,run_length_antero,run_length_retro,\
            directionality,switch,p_directionality_GO,switch_a_to_r,switch_r_to_a,switch_var_STOP,switch_normal,pausing_time_antero,\
                pausing_time_retro,pausing_time_switch = ([] for _ in range(16))

    if settings['theta']:

        theta_std_GO = []
        theta_std_STOP = []

    condition = []
    animal = []
    slide = []

    for file in set(phase_parameters.file.unique()):
        print(f'Per trajectory calculations of trajectories in file {file}')
        for item in set(phase_parameters[(phase_parameters.file==file)].rejoined_trajectory):
            
            data = phase_parameters[(phase_parameters.file==file) & (phase_parameters.rejoined_trajectory==item)]
            data = data.reset_index(drop = True)

            #STOP phases
            data_STOP=data.loc[data['phase']==0,:]
            data_STOP = data_STOP.reset_index(drop = True)

            #GO phases
            data_GO=data.loc[data['phase']==2,:]
            data_GO = data_GO.reset_index(drop = True)

            if len(data_GO)==0: #Check if trajectory contains at least one GO phase
                continue

            if settings['antero_retro']:
                #Antero
                data_GO_antero=data_GO.loc[data_GO['curvilign_velocity']>0,:]
                data_GO_antero = data_GO_antero.reset_index(drop = True)
                #Retro
                data_GO_retro=data_GO.loc[data_GO['curvilign_velocity']<0,:]
                data_GO_retro = data_GO_retro.reset_index(drop = True)

                #Discard spurrious trajectories
                if np.mean(np.abs(data_GO_antero.curvilign_velocity)) >= 4 or np.mean(np.abs(data_GO_retro.curvilign_velocity)) >= 4 :
                    log_analysis['r_speed'] += 1
                    continue

            else:
                if np.mean(np.abs(data_GO.curvilign_velocity)) >= 4 :
                    log_analysis['r_speed'] += 1
                    continue

            iteration.append(item)
            file_list.append(file)

            condition.append(data.loc[data.rejoined_trajectory==item,'condition'].unique()[0])
            animal.append(data.loc[data.rejoined_trajectory==item,'animal'].unique()[0])
            slide.append(data.loc[data.rejoined_trajectory==item,'slide'].unique()[0])

            #Intensity
            traj_intensity_GO = np.mean(data_GO.intensity)
            intensity_GO.append(traj_intensity_GO)
            traj_intensity_STOP = np.mean(data_STOP.intensity)
            intensity_STOP.append(traj_intensity_STOP)

            #Variance
            traj_variance_GO = np.mean(data_GO.variance)
            variance_GO.append(traj_variance_GO)
            traj_variance_STOP = np.mean(data_STOP.variance)
            variance_STOP.append(traj_variance_STOP)

            #Diagonal size
            min_x = np.min(data.min_x)
            max_x = np.max(data.max_x)

            min_y = np.min(data.min_y)
            max_y = np.max(data.max_y)

            delta_x = max_x-min_x
            delta_y = max_y-min_y

            traj_diag_size =np.sqrt(delta_x**2+delta_y**2)
            diag_size.append(traj_diag_size)

            #Curvilign size
            curv_size = np.sum(data.run_length)
            curv_length.append(np.abs(curv_size))

            #Pausing frequency
            trajectory_time = np.sum(data.phase_duration)
            n_stop = len(data[data.phase==0])
            number_stop.append(n_stop)

            if (n_stop==0):
                traj_pausing_frequency=0
            else:
                traj_pausing_frequency = 60*n_stop/trajectory_time

            pausing_frequency.append(traj_pausing_frequency)

            #Trajectory duration
            t_duration = (np.sum(data.phase_duration)*0.05)
            duration.append(t_duration)
            
            #Pausing time
            if len(data_STOP)==0:
                traj_pausing_time=0
            else:
                traj_pausing_time = np.mean(data_STOP.phase_duration)

            pausing_time.append(traj_pausing_time)

            #Ratio of moving particles
            n_particles = data.n_particles.tolist()
            moving_particles.append((phase_parameters[(phase_parameters.file==file)].rejoined_trajectory.nunique()/np.unique(n_particles))[0])

            if settings['antero_retro']:

                #Curvilign velocity
                traj_curvilign_velocity_antero = np.mean(data_GO_antero.curvilign_velocity)
                curvilign_velocity_antero.append(traj_curvilign_velocity_antero)
                traj_curvilign_velocity_retro = np.mean(data_GO_retro.curvilign_velocity)
                curvilign_velocity_retro.append(traj_curvilign_velocity_retro)

                #Run length
                traj_run_length_antero = np.mean(data_GO_antero.run_length)
                run_length_antero.append(traj_run_length_antero)
                traj_run_length_retro = np.mean(data_GO_retro.run_length)
                run_length_retro.append(traj_run_length_retro)

                #Processivity
                traj_processivity_antero = np.mean(data_GO_antero.phase_duration)
                processivity_antero.append(traj_processivity_antero)
                traj_processivity_retro = np.mean(data_GO_retro.phase_duration)
                processivity_retro.append(traj_processivity_retro)

                #Fraction of time paused
                if len(data_STOP)==0 or len(data_GO_antero)==0 or len(data_GO_retro)==0:
                    traj_fraction_paused=0
                else:
                    traj_fraction_paused = traj_pausing_time/(traj_pausing_time+traj_processivity_antero+traj_processivity_retro)

                fraction_paused.append(traj_fraction_paused)

                #Directionality
                if (len(data_GO_antero)==0):
                    distance_antero=0
                else:
                    distance_antero = np.sum(data_GO_antero.run_length)
                if (len(data_GO_retro))==0:
                    distance_retro=0
                else:
                    distance_retro = np.sum(data_GO_retro.run_length)

                distance_total = np.abs(distance_antero)+np.abs(distance_retro)
                
                if distance_total==0:
                    continue

                #Phase directionality
                t_p_directionality_GO = len(data_GO_retro)/len(data_GO) 
                p_directionality_GO.append(t_p_directionality_GO)

                #directionality.append(np.abs(distance_antero)/np.abs(distance_total)) #Fraction of anterograde transport
                directionality.append(np.abs(distance_retro)/np.abs(distance_total)) #Fraction of retrograde transport

                #Directionality reversal
                t_switch = 0 #Overall number of reversals
                t_switch_a_to_r = 0 #Anterograde to retrograde reversals
                t_switch_r_to_a = 0 #Retrograde to anterograde reversals
                t_switch_var_STOP = [] #Variance of signal intensity in STOP phases
                t_pausing_time_antero = [] #Pausing time between anterograde phases
                t_pausing_time_retro = [] #Pausing time between retrograde phases
                t_pausing_time_switch = [] #Pausing time between phases of opposite direction

                for p in set(data.trajectory):
                    p_data = data.loc[data.trajectory==p]
                    if len(p_data) > 2:
                        for i in range(len(p_data.phase_number)-2):
                                if (p_data[p_data.phase_number==i].phase.values==2) & (p_data[p_data.phase_number==i+1].phase.values==0) & (p_data[p_data.phase_number==i+2].phase.values==2):
                                    if (p_data[p_data.phase_number==i].curvilign_velocity.values > 0) & (p_data[p_data.phase_number==i+2].curvilign_velocity.values <0):
                                        t_switch += 1
                                        t_switch_a_to_r += 1
                                        t_switch_var_STOP.append(p_data[p_data.phase_number==i+1].variance.values)
                                        t_pausing_time_switch.append(float(p_data[p_data.phase_number==i+1].phase_duration))
                                    elif (p_data[p_data.phase_number==i].curvilign_velocity.values < 0) & (p_data[p_data.phase_number==i+2].curvilign_velocity.values >0):
                                        t_switch += 1
                                        t_switch_r_to_a += 1
                                        t_switch_var_STOP.append(p_data[p_data.phase_number==i+1].variance.values)
                                        t_pausing_time_switch.append(float(p_data[p_data.phase_number==i+1].phase_duration))
                                    elif (p_data[p_data.phase_number==i].curvilign_velocity.values > 0) & (p_data[p_data.phase_number==i+2].curvilign_velocity.values > 0):
                                        t_pausing_time_antero.append(float(p_data[p_data.phase_number==i+1].phase_duration))
                                    elif (p_data[p_data.phase_number==i].curvilign_velocity.values < 0) & (p_data[p_data.phase_number==i+2].curvilign_velocity.values < 0):
                                        t_pausing_time_retro.append(float(p_data[p_data.phase_number==i+1].phase_duration))

                switch.append(t_switch)
                switch_normal.append(t_switch/t_duration)
                switch_a_to_r.append(t_switch_a_to_r)
                switch_r_to_a.append(t_switch_r_to_a)
                switch_var_STOP.append(np.mean(t_switch_var_STOP)\
                     if len(t_switch_var_STOP) > 0 else np.nan)
                pausing_time_antero.append(np.mean(t_pausing_time_antero)\
                     if len(t_pausing_time_antero) > 0 else np.nan)
                pausing_time_retro.append(np.mean(t_pausing_time_retro)\
                     if len(t_pausing_time_retro) > 0 else np.nan)
                pausing_time_switch.append(np.mean(t_pausing_time_switch)\
                     if len(t_pausing_time_switch) > 0 else np.nan)

            else:
                #Curvilign velocity
                traj_curvilign_velocity = np.mean(data_GO.curvilign_velocity)
                curvilign_velocity.append(traj_curvilign_velocity)

                #Run length
                traj_run_length = np.mean(data_GO.run_length)
                run_length.append(traj_run_length)

                #Processivity
                traj_processivity = np.mean(data_GO.phase_duration)
                processivity.append(traj_processivity)

                #Fraction of time paused
                if len(data_STOP)==0 or len(data_GO)==0 :
                    traj_fraction_paused=0
                else:
                    traj_fraction_paused = traj_pausing_time/(traj_pausing_time+traj_processivity)

                fraction_paused.append(traj_fraction_paused)

            if settings['theta']:

                #Theta standard deviation
                t_theta_std_GO = np.mean(data_GO.theta_std)
                theta_std_GO.append(t_theta_std_GO)
                t_theta_std_STOP = np.mean(data_STOP.theta_std)
                theta_std_STOP.append(t_theta_std_STOP)

    data_dict = {'condition':condition,
                'animal':animal,
                'slide':slide,
                'trajectory':iteration,
                'file':file_list,
                'n_stop':number_stop,
                'fraction_moving':moving_particles,
                'pausing_time':pausing_time, 
                'pausing_frequency':pausing_frequency,
                'diag_size':diag_size,
                'duration':duration,
                'intensity_go':intensity_GO,
                'intensity_stop':intensity_STOP,
                'variance_go':variance_GO,
                'variance_stop':variance_STOP,
                'fraction_paused':fraction_paused,
                'curvilign_length':curv_length,}


    if settings['antero_retro']:

        temp_dict = {'curv_velocity_antero':curvilign_velocity_antero,
                    'curv_velocity_retro':curvilign_velocity_retro,
                    'processivity_antero':processivity_antero,
                    'processivity_retro':processivity_retro,
                    'switch':switch,
                    'run_length_antero':run_length_antero,
                    'run_length_retro':run_length_retro,
                    'phase_dir_go':p_directionality_GO,
                    'switch_var_stop':switch_var_STOP,
                    'switch_a_to_r':switch_a_to_r,
                    'switch_r_to_a':switch_r_to_a,
                    'switch_normal':switch_normal,
                    'fraction_moving':moving_particles,
                    'directionality':directionality,
                    'pausing_time_antero':pausing_time_antero,
                    'pausing_time_retro':pausing_time_retro,
                    'pausing_time_switch':pausing_time_switch,}

        data_dict = {**data_dict, **temp_dict}
        
    else:
        temp_dict = {'curvilign_velocity':curvilign_velocity,
                    'processivity':processivity,
                    'run_length':run_length,}
        data_dict = {**data_dict, **temp_dict}
    
    if settings['theta']:
        temp_dict = {'theta_std_go':theta_std_GO,
                    'theta_std_stop':theta_std_STOP,}
        data_dict = {**data_dict, **temp_dict}

    trajectory_parameters = pd.DataFrame(data_dict)

    return trajectory_parameters

def data_extraction(input_folder,parameters,settings):
    """Extracts transport parameters from trajectories contained in a folder.
    
    :param input_folder: Path to a folder of .csv files containing trajectories.
    :type input_folder: str or Path
    :param parameters: Dictionary containing calculation parameters
    :type parameters: dict
    :param settings: Dictionary containing calculation settings
    :type settings: dict
    """        

    output_folder = folder_structure_creation(input_folder)[0]

    if len(input_folder.parents) < 3 : #Guard against os.walk running on an empty folder, if the input folder is placed at the root of a drive
        input_folder = input_folder.parent

    os.makedirs(output_folder)

    phase_parameters = pd.DataFrame()

    for path, subfolder, files in os.walk(input_folder): #Scan entire folder structure for files
        for name in files:
            if name.endswith('_rejoined.csv') == False:  #Check for correct file
                continue

            #Build output file path
            file_path = os.path.join(path, name)
            file_folder_path = os.path.split(path)[0]
            slide_path,slide = os.path.split(file_folder_path)
            animal_path,animal = os.path.split(slide_path)
            condition = os.path.split(animal_path)[1]

            #Define output file names
            phase_parameters_output = output_folder.joinpath("Per phase parameters.csv")
            traj_parameters_output = output_folder.joinpath("Trajectory average parameters.csv") 

            data = pd.read_csv(file_path,sep=csv_sniffer(file_path))
            
            print("Per phase calculations of "+name)
            phase_parameters = pd.concat((phase_parameters,phase_calculations(parameters,data,settings,condition,slide,name,animal)))

    if phase_parameters.empty:
        raise RuntimeError('No trajectories retained during analysis')
    else:
        phase_parameters.to_csv(phase_parameters_output, sep = '\t')

    trajectory_parameters = trajectory_calculations(phase_parameters,settings)

    #Writes results to .csv files
    trajectory_parameters.to_csv(traj_parameters_output ,sep='\t')
    phase_parameters.drop(['min_x','max_x','min_y','max_y'], axis='columns', inplace=True)
    phase_parameters.to_csv(phase_parameters_output, sep = '\t')

    dict_dump(Path(output_folder).parent,log_analysis,'log')


if __name__ == '__main__':
    from pathlib import Path
    import time

    parameters = {
    #Data Extraction
    'r_conf_cut' : 0.64,
    'px' : 0.173, #in µm
    'dt' : 0.05, #in s
    'min_thr_prec' : 50, # in nm
    'sliding_window':3,
    'sigma':129,
    'len_cutoff':30, #Number of points
    'threshold_poly3':1.4, #Deviation from third-degree polynom
    }   
    settings = {
    #Data Extraction
    'polynomial_fit':True,
    'minimization':True,
    'antero_retro':True,
    'theta':True,
    }

    start = time.time()
    input_folder = Path(r"/media/baptiste/Windows/Users/LUMIN10/Documents/video_benchmark_int Results - 20230213_182939/video_benchmark_int")
    data_extraction(input_folder,parameters,settings)
    end = time.time()
    duration = end - start
    print('%dh%s' % (int(duration//3600),f'{int((duration%3600)/60):02d}'))