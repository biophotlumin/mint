"""Functions used to output calculation results into files.
    trajectory_output writes trajectories into a .csv file. 
    image_output plots trajectories onto the first frame of the current file.
    trajectory_separation extracts each trajectory from a DataFrame and writes them into individual .txt files.
    dict_dump dumps a dictionary into a .txt file.
"""
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
from pathlib import Path

def trajectory_output(log,name,process,output):
    """Writes trajectories into a .csv file.

    Uses pandas to_csv to write the DataFrame containing trajectories into a .csv file.

    :param log: Dictionary containing the current output file path.
    :type log: dict
    :param name: Name of the file currently being processed.
    :type name: string
    :param process: sed to differentiate between .csv files written before and after optional trajectory processing.
    :type process: string
    :param output: DataFrame that is being dumped.
    :type output: DataFrame
    """

    filecsv = Path(str(log['output_file_path'])).joinpath(name+process+".csv") #Creates file path and name
    output.to_csv(filecsv,sep='\t')

def image_output(log,name,frames,trajectory,item):
    """Plots trajectories onto the first frame of a file. Inputs a dictionary, a string, a NumPy array, and a DataFrame.

    :param log: Dictionary containing the current output file path.
    :type log: dict
    :param name: Name of the file currently being processed.
    :type name: string
    :param frames: Name of the file currently being processed.
    :type frames: string
    :param trajectory: DataFrame containing the trajectory to be plotted.
    :type trajectory: DataFrame
    :param item: Number of the trajectory to be plotted. Pass False to plot all trajectories.
    :type item: int, boolean
    """  

    #Initializes plot
    plt.imshow(frames[0])
    plt.gca().invert_yaxis() #Plotting the y axis inverts it by default, so it must be inverted again
    bbox_props = dict(boxstyle="round", fc="w", ec="w")
    arrow_props = dict(arrowstyle="simple")

    trajectories = []
    if isinstance(item, bool):
        plt.title('Trajectories from '+name)
        filepng = Path(str(log['output_file_path'])).joinpath(name+".png") #Creates file path and name
        for item in set(trajectory.rejoined_particle):
            trajectories.append(item)
    else:
        plt.title(f'Trajectory {item} from {name}')
        trajectories.append(item)
        filepng = Path(str(log['output_file_path'])).joinpath(name+"-"+str(item)+".png")

    for item in trajectories: #Loops for each sub trajectory
        x = trajectory[trajectory.rejoined_particle==item].x
        y = trajectory[trajectory.rejoined_particle==item].y
        plt.plot(x,y,lw=0.5)
        if len(trajectories)>1:
            x0,y0 = x.iloc[0],y.iloc[0]
            plt.annotate(str(item),xy=(x0,y0),xytext=(x0+30,y0+30),bbox=bbox_props,arrowprops=arrow_props, size = 4)

    plt.savefig(filepng,dpi=300)
    plt.close()

def trajectory_separation(log,name,trajectory,settings,sub):
    """Separates individual trajectories and writes them into .txt files. 

    For each trajectory contained in the DataFrame, removes unnecessary columns, reindexes the remaining ones, and writes it into a .txt file.

    :param log: Dictionary containing the current output file path.
    :type log: dict
    :param name: Name of the file currently being processed.
    :type name: string
    :param trajectory: Number of the trajectory being plotted.
    :type trajectory: int
    :param settings: Dictionary specifying wether or not a specific process should be applied.
    :type settings: dict
    :param sub: Subset of the DataFrame containing the trajectories.
    :type sub: DataFrame
    """    

    file = Path(str(log['output_file_path'])).joinpath(name+"-"+str(trajectory)+".txt") #Creates file path and name
    pd.options.mode.chained_assignment = None
    if settings['SNR_estimation'] == True: #Checks wether SNR_estimation was used
        sub.drop(sub.columns.difference(['y','x','mass','feet']), axis='columns', inplace=True) #If it was, drops the following columns
        sub.columns = ['Yraw','Xraw','mass','feet']
        sub = sub.reindex(columns = ['Xraw','Yraw','mass','feet'])
        sub = sub.rename(columns={"mass":"Signal" , "feet":"Noise"})
    else:
        sub.drop(sub.columns.difference(['y','x','mass']), axis='columns', inplace=True) #If not, drops only the columns that aren't added by SNR_estimation
        sub.columns = ['Yraw','Xraw','mass']
        sub = sub.reindex(columns = ['Xraw','Yraw','mass'])
        sub = sub.rename(columns={"mass":"Signal"})

    sub.to_csv(file, sep = '\t',index = False)

def dict_dump(log,dict,file_name):
    """Writes the content of a dictionary into a .txt file

    :param log: Dictionary containing the current output file path.
    :type log: dict
    :param dict: Dictionary to dump.
    :type dict: dict
    :param file_name: Name of the text file.
    :type file_name: string
    """    

    with open(Path(log['output_folder']).joinpath(str(file_name)+".txt"), 'w') as dict_txt:
        for k, v in dict.items():
            print(str(k)+" : "+str(v), file=dict_txt)

def generate_report():
    # TODO Generate report with parameters, metrics and figures using ReportLab
    pass