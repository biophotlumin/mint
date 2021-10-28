"""Functions used to output calculation results into files.

    trajectory_output writes trajectories into a .csv file. 
    image_output plots a trajectory onto the first frame of the current file.
    trajectory_separation extracts each trajectory from a DataFrame and writes them into individual .txt files.
    final_image_output plots all trajectories found in a specific file into the first frame of said file.
"""
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def trajectory_output(log,name,process,planned_output):
    """Writes trajectories found in a video file into a .csv file. Inputs a dictionary, two strings, and a DataFrame.

                Uses pandas to_csv to write the DataFrame containing trajectories into a .csv file.
                parameters is a dictionary containing the current output file path.
                process is a string used to differentiate between .csv files written before and after optional trajectory processing. 
                name is the name of the file currently being processed.
                planned_output is the DataFrame that is written.
    """ 
    #filecsv = str(log['output_file_path'])+"\\"+name+process+".csv" #Creates file path and name
    filecsv = Path(str(log['output_file_path'])).joinpath(name+process+".csv")
    planned_output.to_csv(filecsv,sep='\t')
    
def image_output(log,name,frames,processed_trajectory,item,trajectory_number):
    """Plots a trajectory onto the first frame of a file. Inputs a dictionary, a string, a NumPy array, a DataFrame, and two integers.

                Uses matplotlib to plot a trajectory onto the first frame of a file.
                parameters is a dictionary containing the current output file path.
                name is the name of the file currently being processed.
                frames is a NumPy array containing raw frames from the video file.
                processed_trajectory is a DataFrame containing the trajectory to be plotted.
                item is the number of the trajectory being plotted.
                trajectory_number is the number of the trajectory being plotted.
    """ 
    #filepng = str(log['output_file_path'])+"\\"+name+str(trajectory_number)+".png" #Creates file path and name
    filepng = Path(str(log['output_file_path'])).joinpath(name+"-"+str(trajectory_number)+".png")
    pd.set_option("mode.chained_assignment", None)
    
    #Initializes plot
    plt.title('Trajectory from '+name)
    plt.imshow(frames[0])
    plt.gca().invert_yaxis() #Plotting the y axis inverts it by default, so it must be inverted again        
    bbox_props = dict(boxstyle="round", fc="w", ec="w")
    arrow_props = dict(arrowstyle="simple")
    x = processed_trajectory[processed_trajectory.particle==item].x
    y = processed_trajectory[processed_trajectory.particle==item].y
    plt.plot(x,y,lw=0.5)
    x0,y0 = x.iloc[0],y.iloc[0]
    plt.annotate(trajectory_number,xy=(x0,y0),xytext=(x0+30,y0+30),bbox=bbox_props,arrowprops=arrow_props, size = 4) #Plots arrows with number pointing to each trajectory
    plt.savefig(filepng,dpi=300)
    plt.close()

def trajectory_separation(log,name,trajectory_number,settings,sub_trajectory):
    """Separates individual trajectories and writes them into .txt files. Inputs a dictionary, a string, an integer, and dictionary, and a DataFrame.

                For each trajectory contained in the DataFrame, removes unnecessary columns, reindexes the remaining ones, and writes it into a .txt file.
                parameters is a dictionary containing the current output file path.
                name is the name of the file currently being processed.
                trajectory_number is the number of the trajectory being plotted.
                settings is a dictionary, containing Boolean values regarding wether or not a specific process should be applied.
                sub_trajectory is a subset of the DataFrame containing the trajectories.
        """
    file = Path(str(log['output_file_path'])).joinpath(name+"-"+str(trajectory_number)+".txt") #Creates file path and name
    pd.options.mode.chained_assignment = None
    if settings['SNR_estimation'] == True: #Checks wether SNR_estimation was used
        sub_trajectory.drop(sub_trajectory.columns.difference(['y','x','mass','feet']), axis='columns', inplace=True) #If it was, drops the following columns
        sub_trajectory.columns = ['Yraw','Xraw','mass','feet']
        sub_trajectory = sub_trajectory.reindex(columns = ['Xraw','Yraw','mass','feet'])
        sub_trajectory = sub_trajectory.rename(columns={"mass":"Signal" , "feet":"Noise"})
    else:
        sub_trajectory.drop(sub_trajectory.columns.difference(['y','x','mass']), axis='columns', inplace=True) #If not, drops only the columns that aren't added by SNR_estimation
        sub_trajectory.columns = ['Yraw','Xraw','mass']
        sub_trajectory = sub_trajectory.reindex(columns = ['Xraw','Yraw','mass'])
        sub_trajectory = sub_trajectory.rename(columns={"mass":"Signal"})

    sub_trajectory.to_csv(file, sep = '\t',index = False)

def final_image_ouput(log,name,frames,processed_trajectory):
    """Plots all trajectories onto the first frame of a file. Inputs a dictionary, a string, a NumPy array, and a DataFrame.

                Functionally identical to image_output, with the exception of the sub trajectory for loop being included within the function itself.
                parameters is a dictionary containing the current output file path.
                name is the name of the file currently being processed.
                frames is a NumPy array containing raw frames from the video file.
                processed_trajectory is a DataFrame containing the trajectory to be plotted.
    """
    filepng = Path(str(log['output_file_path'])).joinpath(name+".png") #Creates file path and name

    #Initializes plot
    plt.title('Trajectories from '+name)
    plt.imshow(frames[0])
    plt.gca().invert_yaxis() #Plotting the y axis inverts it by default, so it must be inverted again
    bbox_props = dict(boxstyle="round", fc="w", ec="w")
    arrow_props = dict(arrowstyle="simple")

    for item in set(processed_trajectory.particle): #Loops for each sub trajectory
        x = processed_trajectory[processed_trajectory.particle==item].x
        y = processed_trajectory[processed_trajectory.particle==item].y
        plt.plot(x,y,lw=0.5)
        x0,y0 = x.iloc[0],y.iloc[0]
        plt.annotate(str(item),xy=(x0,y0),xytext=(x0+30,y0+30),bbox=bbox_props,arrowprops=arrow_props, size = 4)

    plt.savefig(filepng,dpi=300)
    plt.close()
