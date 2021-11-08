"""Functions used to extract transport parameters from trajectories previously obtained.

    inst_velocity calculates the instantaneous velocity between each point of a trajectory.
    confinement calculates the confinement ratio of a trajectory to segment it into phases.
    phase_calculations calculates per phase variables of interest.
    phase_calculations_antero_retro behaves like phase_calculations but can distinguish between anterograde and retrograde transport.
    trajectory_calculations averages per phase data and calculates per trajectory variables of interest.
    trajectory_calculations_antero_retro behaves like trajectory_calculations but can distinguish between anterograde and retrograde transport.
    data_extraction scans the output folder for .csv files and calls the aforementioned functions.
"""
#Imports
import numpy as np
import pandas as pd
import os
from scipy import *
from scipy.fftpack import *
from scipy.signal import *
from matplotlib.pyplot import *
from utils import folder_structure_creation
from trajectory_calculations import *

def inst_velocity(x,y,parameters):
    """Calculates the instantaneous velocity between each point of a trajectory. Inputs two floats and a dictionary, returns a NumPy array.
    
        x and y are coordinates (in pixels).
        parameters is a dictionary of calculation parameters, as defined in script.py.
    """
    size = len(x)
    dt = parameters['dt']
    v_inst = np.zeros(size)
    v_inst[0] = (np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2))/dt #Instantaneous speed of the first point calculated from the first segment
    v_inst[size-1] = (np.sqrt((x[size-1]-x[size-2])**2+(y[size-1]-y[size-2])**2))/dt #Instantaneous speed of the last point calculated from the last segment

    for u in range(size-2):
        u = u+1
        v_inst[u] = (np.sqrt((x[u+1]-x[u-1])**2+(y[u+1]-y[u-1])**2))/(2*dt)

    return v_inst

def confinement(x,y,parameters):
    """Calculates the confinement ratio of a trajectory to segment it into phases. Inputs two floats and a dictionary, returns a list.
        x and y are coordinates (in pixels).
        parameters is a dictionary of calculation parameters, as defined in script.py.
    """

    r_conf = []
    size=len(x)
    sliding_window = parameters['sliding_window']
    if (size>=2*sliding_window+1):
        for i in range(sliding_window): #Extrapolate the first 5 points
            d_net = np.sqrt((x[sliding_window*2]-x[0])**2+(y[sliding_window*2]-y[0])**2)
            d_total = 0
            for u in range(2*sliding_window):
                d_total = d_total + np.sqrt((x[u+1]-x[u])**2+(y[u+1]-y[u])**2)
            r_conf.append((d_net/d_total)**2)

        for i in range (sliding_window,size-sliding_window):
            if (i == size-sliding_window-1): #Extrapolate the last 5 points
                d_total = 0
                d_net = np.sqrt((x[i+sliding_window]-x[i-sliding_window])**2+(y[i+sliding_window]-y[i-sliding_window])**2)
                for j in range (2*sliding_window):
                    d_total = d_total + np.sqrt((x[(i-sliding_window)+j+1]-x[(i-sliding_window)+j])**2+(y[(i-sliding_window)+j+1]-y[(i-sliding_window)+j])**2)
                if (d_total!=0):
                    for p in range(sliding_window+1):
                        r_conf.append((d_net/d_total)**2)
                else:
                    for p in range(sliding_window+1):
                        r_conf.append(0)
            else:

                d_total = 0
                d_net = np.sqrt((x[i+sliding_window]-x[i-sliding_window])**2+(y[i+sliding_window]-y[i-sliding_window])**2)
                for j in range (2*sliding_window):
                    d_total = d_total + np.sqrt((x[(i-sliding_window)+j+1]-x[(i-sliding_window)+j])**2+(y[(i-sliding_window)+j+1]-y[(i-sliding_window)+j])**2)
                if (d_total!=0):
                        r_conf.append((d_net/d_total)**2)
                else:
                        r_conf.append(0)
    else:
        pass
    return r_conf

def phase_calculations(parameters,data,condition,slide,name):
    """Calculates per phase variables of interest. Inputs a dictionary, a DataFrame, and three strings. Returns a DataFrame.

        parameters is a dictionary of calculation parameters, as defined in script.py.
        data is a DataFrame containing data from a .csv file being processed.
        condition, slide and name are different levels of the input folder structure. 
    """

    dt = parameters['dt']
    conf_threshold = parameters['min_theoretical_precision']*1e-3/dt # in µm/s, threshold defined by a minimal velocity of 10nm/dt (or 0.2µm/s)
    f_phase_parameters = pd.DataFrame()

    for trajectory in set(data.particle):
                
        subdata = data[data.particle==trajectory]
        subdata = subdata.reset_index(drop = True)

        x = subdata.x
        x = x.dropna()
        
        y = subdata.y
        y = y.dropna()
        
        x = x.reset_index(drop = True)
        y = y.reset_index(drop = True)

        size = len(x)
        
        r_conf = confinement(x,y,parameters) #Separate trajectory into phases

        #Switch from pixels to µm   
        x = x*parameters['px']
        x = x.dropna()
        y = y*parameters['px']
        y = y.dropna()

        v_inst = inst_velocity(x,y,parameters) #Get instantaneous velocity for each point

        phase = np.zeros(size)

        #Categorize each phase as either a GO or STOP phase
        if (size>=2*parameters['sliding_window']+1):
            for i in range(len(r_conf)):
                if (r_conf[i]>parameters['r_conf_cut']):
                        phase[i] = 2 #GO phase

            for i in range(size): #STOP phase refinment
                if (phase[i]==2) & (v_inst[i]<conf_threshold):
                    phase[i] = 0
        else:
            
            for i in range(size): #STOP phase refinment
                phase[i]==2 
                if (phase[i]==2) & (v_inst[i]<conf_threshold):
                    phase[i] = 0

        #Get min and max values for both X and Y to calculate length of rejoined trajectories later
        min_x = np.min(x)
        max_x = np.max(x)

        min_y = np.min(y)
        max_y = np.max(y)
        
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

            curvilign_velocity = np.abs(np.mean(sub_v_inst))
            vectorial_velocity = np.abs((np.sqrt((x[stop-1]-x[start])**2+(y[stop-1]-y[start])**2))/(dt*phase_length))
            
            if (phase[start]==0):
                phase_sign = 0
            if (phase[start]==2):
                phase_sign = 2

            run_length = curvilign_velocity*phase_duration
            f_phase_parameters.reset_index(inplace=True, drop=True)

            f_phase_parameters = f_phase_parameters.append([{'trajectory': trajectory, 'phase':phase_sign,'phase_number':phase_number,'phase_length':phase_length,\
            'vectorial_velocity':vectorial_velocity,'curvilign_velocity':curvilign_velocity,'phase_duration':phase_duration,'run_length':run_length,'intensity':mean_intensity,\
            'variance':variance,'condition':condition,'slide':slide,'rejoined_trajectory':subdata.rejoined_particle.unique()[0],\
                'file':name,'min_x':min_x,'max_x':max_x,'min_y':min_y,'max_y':max_y,'n_particles':data.n_particles.unique()[0]}])

    return f_phase_parameters

def phase_calculations_antero_retro(parameters,data,settings,condition,slide,name,animal):
    """Calculates per phase variables of interest. Inputs a dictionary, a DataFrame, and four strings. Returns a DataFrame.

        parameters is a dictionary of calculation parameters, as defined in script.py.
        data is a DataFrame containing data from a .csv file being processed.
        condition, slide, name and animal are different levels of the input folder structure. 
    """
    dt = parameters['dt']
    conf_threshold = parameters['min_theoretical_precision']*1e-3/dt # in µm/s, threshold defined by a minimal velocity of 10nm/dt (or 0.2µm/s)
    f_phase_parameters = pd.DataFrame()

    for trajectory in set(data.particle):
                
        subdata = data[data.particle==trajectory]
        subdata = subdata.reset_index(drop = True)

        if settings['polynomial_fit']==True:
            if polynomial_fit(subdata,parameters) == True:
                pass
            else:
                continue

        if settings['minimization']==True:
            subdata = minimization(subdata,parameters)

        x = subdata.x
        x = x.dropna()
        
        y = subdata.y
        y = y.dropna()
        
        x = x.reset_index(drop = True)
        y = y.reset_index(drop = True)

        size = len(x)
        
        r_conf = confinement(x,y,parameters) #Separate trajectory into phases

        #Switch from pixels to µm   
        x = x*parameters['px']
        x = x.dropna()
        y = y*parameters['px']
        y = y.dropna()

        v_inst = inst_velocity(x,y,parameters) #Get instantaneous velocity for each point

        phase = np.zeros(size)

        #Categorize each phase as either a GO or STOP phase
        if (size>=2*parameters['sliding_window']+1):
            for i in range(len(r_conf)):
                if (r_conf[i]>parameters['r_conf_cut']):
                        phase[i] = 2 #GO phase

            for i in range(size): #STOP phase refinment
                if (phase[i]==2) & (v_inst[i]<conf_threshold):
                    phase[i] = 0
        else:
            
            for i in range(size): #STOP phase refinment
                phase[i]==2 
                if (phase[i]==2) & (v_inst[i]<conf_threshold):
                    phase[i] = 0
        
        #Get min and max values for both X and Y to calculate length of rejoined trajectories later
        min_x = np.min(x)
        max_x = np.max(x)

        min_y = np.min(y)
        max_y = np.max(y)
        
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

        for phase_number in range(len(cut)-1): #Per phase processing. '-1' : First and last phases are deleted
            start = cut[phase_number]
            stop = cut[phase_number+1]
            sub_phase = subdata.loc[start:stop-1]
            intensity = sub_phase.mass

            sub_v_inst = sub_phase.Vinst 
            phase_length = len(sub_phase)
            phase_duration = phase_length*dt

            #if (np.max(sub_v_inst)>5):
                #break

            variance = (np.std(intensity))**2
            mean_intensity = np.mean(intensity)

            if (phase_length==1):
                variance = 0

            curvilign_velocity = np.abs(np.mean(sub_v_inst))
            vectorial_velocity = np.abs((np.sqrt((x[stop-1]-x[start])**2+(y[stop-1]-y[start])**2))/(dt*phase_length))
            
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
            f_phase_parameters.reset_index(inplace=True, drop=True)

            f_phase_parameters = f_phase_parameters.append([{'trajectory': trajectory, 'phase':phase_sign,'phase_number':phase_number,'phase_length':phase_length,\
            'vectorial_velocity':vectorial_velocity,'curvilign_velocity':curvilign_velocity,'phase_duration':phase_duration,'run_length':run_length,'intensity':mean_intensity,\
            'variance':variance,'condition':condition,'slide':slide,'rejoined_trajectory':subdata.rejoined_particle.unique()[0],'animal':animal,\
                'file':name,'min_x':min_x,'max_x':max_x,'min_y':min_y,'max_y':max_y,'n_particles':data.n_particles.unique()[0]}])
    return f_phase_parameters

def trajectory_calculations(phase_parameters):
    """Averages per phase data and calculates per trajectory variables of interest. Inputs a DataFrame, returns a DataFrame.

        phase_parameters is a DataFrame containing data as defined in phase_calculations.
    """

    iteration = []
    file_list=[]

    phase_parameters.sort_values(by=['file'],inplace=True)
    phase_parameters.sort_values(by=['rejoined_trajectory'],inplace=True)

    intensity,variance,curvilign_velocity,processivity,run_length,pausing_frequency,pausing_time,\
        diag_size,fraction_paused,moving_particles,number_stop = ([] for i in range(11)) #Initialize lists

    mean_condition = pd.DataFrame({'condition':[]},dtype=str)
    mean_slide = pd.DataFrame({'slide':[]},dtype=str)

    for file in set(phase_parameters.file.unique()):
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

            iteration.append(item)
            file_list.append(file)

            testcond = pd.Series(data=[data.loc[data.rejoined_trajectory==item,'condition'].values[0]], index=['condition'])
            testslide = pd.Series(data=[data.loc[data.rejoined_trajectory==item,'slide'].values[0]], index=['slide'])

            mean_condition = mean_condition.append(testcond,ignore_index=True)
            mean_slide = mean_slide.append(testslide,ignore_index=True)

            #Intensity
            traj_intensity = np.mean(data.intensity)
            intensity.append(traj_intensity)

            #Variance
            traj_variance = np.mean(data.variance)
            variance.append(traj_variance)
            
            #Diagonal size
            min_x = np.min(data.min_x)
            max_x = np.max(data.max_x)

            min_y = np.min(data.min_y)
            max_y = np.max(data.max_y)

            delta_x = max_x-min_x
            delta_y = max_y-min_y

            traj_diag_size =np.sqrt(delta_x**2+delta_y**2)
            diag_size.append(traj_diag_size)

            #Pausing frequency
            trajectory_time = np.sum(data.phase_duration)
            n_stop = len(data[data.phase==0])
            number_stop.append(n_stop)

            if (n_stop==0):
                traj_pausing_frequency=0
            else:
                traj_pausing_frequency = 60*n_stop/trajectory_time

            pausing_frequency.append(traj_pausing_frequency)

            #Curvilign velocity
            traj_curvilign_velocity = np.mean(data_GO.curvilign_velocity)
            curvilign_velocity.append(traj_curvilign_velocity)

            #Run length
            traj_run_length = np.mean(data_GO.run_length)
            run_length.append(traj_run_length)

            #Processivity
            traj_processivity = np.mean(data_GO.phase_duration)
            processivity.append(traj_processivity)
            
            #Pausing time
            if len(data_STOP)==0:
                traj_pausing_time=0
            else:
                traj_pausing_time = np.mean(data_STOP.phase_duration)

            pausing_time.append(traj_pausing_time)

            #Fraction of time paused
            if len(data_STOP)==0 or len(data_GO)==0 :
                traj_fraction_paused=0
            else:
                traj_fraction_paused = traj_pausing_time/(traj_pausing_time+traj_processivity)

            fraction_paused.append(traj_fraction_paused)

            #Ratio of moving particles
            n_particles = data.n_particles.tolist()
            moving_particles.append((phase_parameters[(phase_parameters.file==file)].rejoined_trajectory.nunique()/np.unique(n_particles))[0])

    cond_list = mean_condition['condition'].tolist()
    slide_list = mean_slide['slide'].tolist()

    trajectory_parameters = pd.DataFrame({'pausing_time':pausing_time, 'pausing_frequency':pausing_frequency,\
         'curvilign_velocity':curvilign_velocity,'intensity':intensity,'variance':variance,\
            'proc':processivity,'run_length':run_length, 'diag_size':diag_size,'fraction_paused':fraction_paused,\
                'trajectory':iteration,'n_stop':number_stop,'fraction_moving':moving_particles,\
                    'condition':cond_list,'slide':slide_list,'file':file_list})
    return trajectory_parameters

def trajectory_calculations_antero_retro(phase_parameters):
    """Averages per phase data and calculates per trajectory variables of interest. Inputs a DataFrame, returns a DataFrame.

        phase_parameters is a DataFrame containing data as defined in phase_calculations_antero_retro.
    """
    iteration = []
    file_list=[]

    phase_parameters.sort_values(by=['file'],inplace=True)
    phase_parameters.sort_values(by=['rejoined_trajectory'],inplace=True)
    intensity_GO,intensity_STOP,variance_GO,variance_STOP,curvilign_velocity_antero,curvilign_velocity_retro,\
        processivity_antero,processivity_retro,run_length_antero,run_length_retro,pausing_frequency,pausing_time,\
        diag_size,fraction_paused,directionality,moving_particles,number_stop,t_switch = ([] for i in range(18)) #Initializes lists

    mean_condition = pd.DataFrame({'condition':[]},dtype=str)
    mean_animal = pd.DataFrame({'animal':[]},dtype=str)
    mean_slide = pd.DataFrame({'slide':[]},dtype=str)

    for file in set(phase_parameters.file.unique()):
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

            #Antero
            data_GO_antero=data_GO.loc[data_GO['curvilign_velocity']>0,:]
            data_GO_antero = data_GO_antero.reset_index(drop = True)
            #Retro
            data_GO_retro=data_GO.loc[data_GO['curvilign_velocity']<0,:]
            data_GO_retro = data_GO_retro.reset_index(drop = True)

            if np.mean(np.abs(data_GO_antero.curvilign_velocity)) >= 4 or np.mean(np.abs(data_GO_retro.curvilign_velocity)) >= 4 :
                continue

            iteration.append(item)
            file_list.append(file)

            testcond = pd.Series(data=[data.loc[data.rejoined_trajectory==item,'condition'].values[0]], index=['condition'])
            testslide = pd.Series(data=[data.loc[data.rejoined_trajectory==item,'slide'].values[0]], index=['slide'])
            testanimal = pd.Series(data=[data.loc[data.rejoined_trajectory==item,'animal'].values[0]], index=['animal'])

            mean_condition = mean_condition.append(testcond,ignore_index=True)
            mean_slide = mean_slide.append(testslide,ignore_index=True)
            mean_animal = mean_animal.append(testanimal,ignore_index=True)

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

            #Pausing frequency
            trajectory_time = np.sum(data.phase_duration)
            n_stop = len(data[data.phase==0])
            number_stop.append(n_stop)

            if (n_stop==0):
                traj_pausing_frequency=0
            else:
                traj_pausing_frequency = 60*n_stop/trajectory_time

            pausing_frequency.append(traj_pausing_frequency)

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
            
            #Pausing time
            if len(data_STOP)==0:
                traj_pausing_time=0
            else:
                traj_pausing_time = np.mean(data_STOP.phase_duration)

            pausing_time.append(traj_pausing_time)

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

            #directionality.append(np.abs(distance_antero)/np.abs(distance_total))
            directionality.append(np.abs(distance_retro)/np.abs(distance_total))

            #Ratio of moving particles
            n_particles = data.n_particles.tolist()
            moving_particles.append((phase_parameters[(phase_parameters.file==file)].rejoined_trajectory.nunique()/np.unique(n_particles))[0])

            switch = 0
            for p in set(data.trajectory):
                p_data = data.loc[data.trajectory==p]
                if len(p_data) > 2:
                    for i in range(len(p_data.phase_number)-2):
                                if (p_data[p_data.phase_number==i].phase.values==2) & (p_data[p_data.phase_number==i+1].phase.values==0) & (p_data[p_data.phase_number==i+2].phase.values==2):
                                    if (p_data[p_data.phase_number==i].curvilign_velocity.values > 0) & (p_data[p_data.phase_number==i+2].curvilign_velocity.values <0):
                                        switch += 1
                                    elif (p_data[p_data.phase_number==i].curvilign_velocity.values < 0) & (p_data[p_data.phase_number==i+2].curvilign_velocity.values >0):
                                        switch += 1
                    
            t_switch.append(switch)

    cond_list = mean_condition['condition'].tolist()
    slide_list = mean_slide['slide'].tolist()
    animal_list = mean_animal['animal'].tolist()

    trajectory_parameters = pd.DataFrame({'pausing_time':pausing_time, 'pausing_frequency':pausing_frequency,\
         'curvilign_velocity_antero':curvilign_velocity_antero,'curvilign_velocity_retro':curvilign_velocity_retro,\
          'intensity_GO':intensity_GO,'intensity_STOP':intensity_STOP,'variance_GO':variance_GO,'variance_STOP':variance_STOP,\
           'processivity_antero':processivity_antero,'processivity_retro':processivity_retro,'switch':t_switch,\
            'run_length_antero':run_length_antero,'run_length_retro':run_length_retro,\
             'diag_size':diag_size,'fraction_paused':fraction_paused,'trajectory':iteration,\
              'n_stop':number_stop,'fraction_moving':moving_particles,'directionality':directionality,
                'condition':cond_list,'slide':slide_list,'animal':animal_list,'file':file_list})

    return trajectory_parameters

def data_extraction(parameters,input_folder,settings):
    """Runs through .csv files and calls the appropriate functions to extract transport parameters. Inputs a dictionary, a string and a dictionary.

        parameters is a dictionary of calculation parameters, as defined in script.py.
        input_folder is the path to the input folder.
        settings is a dictionary containing boolean values for optional data processing, as defined in script.py.
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

            if settings['antero_retro']==True:
                animal_path,animal = os.path.split(slide_path)
                condition = os.path.split(animal_path)[1]
            else:
                condition = os.path.split(slide_path)[1]

            #Define output file names
            phase_parameters_output = output_folder.joinpath("Per phase parameters.csv")
            traj_parameters_output = output_folder.joinpath("Trajectory average parameters.csv") 

            data = pd.read_csv(file_path,sep='\t')
            
            print("Per phase calculations of "+name)
            if settings['antero_retro']:
                if name == '210114_nKTP_dynapyrazole.lif - Series059.tif':
                    parameters['dt'] = 0.18
                    print(parameters['dt'])
                    phase_parameters = phase_parameters.append(phase_calculations_antero_retro(parameters,data,settings,condition,slide,name,animal))
                else:
                    parameters['dt'] = 0.1
                    phase_parameters = phase_parameters.append(phase_calculations_antero_retro(parameters,data,settings,condition,slide,name,animal))
            else:
                phase_parameters = phase_parameters.append(phase_calculations(parameters,data,condition,slide,name))
    phase_parameters.to_csv(phase_parameters_output, sep = '\t')
    print("Per trajectory calculations of "+name)
    if settings['antero_retro']==True:
        trajectory_parameters = trajectory_calculations_antero_retro(phase_parameters)
    else:
        trajectory_parameters = trajectory_calculations(phase_parameters)

    #Writes results to .csv files
    trajectory_parameters.to_csv(traj_parameters_output ,sep='\t')
    phase_parameters.drop(['min_x','max_x','min_y','max_y'], axis='columns', inplace=True)
    phase_parameters.to_csv(phase_parameters_output, sep = '\t')


if __name__ == '__main__':
    from pathlib import Path
    import time

    parameters = {
    #trackpy.batch
    'diameter':5,
    'minmass':180,
    'separation':10,
    #trackpy.link
    'search_range':19,
    'memory':5,
    'adaptive_stop':12,
    'adaptive_step':0.95,
    'stub_filtering':3,
    #trackpy.motion.msd
    'threshold':300,
    #SNR estimation
    'base_level':0,
    #Rejoining
    'threshold_t':10,
    'threshold_r':60,
    #Data Extraction
    'r_conf_cut' : 0.9**2,
    'px' : 0.173, #in µm
    'dt' : 0.05, #in s
    'min_theoretical_precision' : 50, # in nm
    'sliding_window':3,
    'sigma':129,
    'len_cutoff':30, #Number of points
    'threshold_poly3':1.4 #Deviation from third-degree polynom
    }   
    settings = {
    #Denoising
    'tophat':True,
    'wavelet':False,
    #Processing
    'stub_filtering':False,
    'MSD':True,
    'rejoining':True,
    'SNR_estimation':True,
    #Outputs
    'individual_images':True,
    'individual_txt':True,
    'group_image':True,
    #Data Extraction
    'polynomial_fit':True,
    'minimization':True,
    'antero_retro':True
    }

    start = time.time()
    input_folder = Path(r"/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20211104_102617/124")
    data_extraction(parameters,input_folder,settings)
    end = time.time()
    print((end-start)/60)

    """data = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20210921_182942 line average 4/124 Results - 20211028_174342/Per phase parameters.csv',sep='\t')
    rev = trajectory_calculations_antero_retro(data)
    rev.to_csv(r'/home/baptiste/rev3.csv',sep='\t')"""