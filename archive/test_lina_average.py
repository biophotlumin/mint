#Imports
import numpy as np
import os
import imageio
import trackpy as tp
import trajectory_calculations as ft
from image_denoising import *
from output_files_creation import *
import pathlib as pl
import pandas as pd
"""
def tracking(input_folder):

    for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
        for name in files:
            if name.endswith('.tif') == False:  #Check for correct file extension
                continue #Skips to next file if not .tif           

            #Building output folder
            file_path = os.path.join(path, name) #Get file path of current file
            #Opening video file
            frames = imageio.volread(file_path)

            #mean = ((frames[len(frames)-1]+frames[len(frames)-2])/2)
            mean = ((frames[29]+frames[30])/2)
            print(frames.shape)
            frames = line_average(frames)
            print(frames.shape)
            #delta = mean - frames[30]
            #print(delta)"""

"""def line_average(frames):
    len_init = frames.shape[0]
    for i in range(frames.shape[0]-2):
        if i==(int(len_init/2)):
            break
        for j in range((frames.shape[1]-2)):
            frames[i-1,j,] = ((frames[i-1,j,]+frames[i,j,])/2)
        frames = np.delete(frames,obj=i+1,axis=0)
    return frames"""


def line_average(frames):
    avg_frames = np.empty((int(frames.shape[0]/2),frames.shape[1],frames.shape[2]))

    for i in range(0,(frames.shape[0]-2),2):
        for j in range(frames.shape[1]):
            avg_frames[(int(i/2)),j,] = ((frames[i,j,]+frames[i+1,j,])/2)
    return avg_frames

"""input_folder = r'/media/baptiste/SHG_tracking_data/Jitter nKTP'
tracking(input_folder)"""

"""for i in range(0,100):
    print('I'+str(i))
    for j in range(0,100):
        print(j)
        if j == 80:
    
            break"""

from matplotlib import colors
import matplotlib.pyplot as plt
import pandas as pd
import os

"""input_folder = r'/home/baptiste/Documents/Pamela/ResultatsTXT'
traj_n = 1
plt.gca().invert_yaxis()
ax = plt.axes()
ax.set_facecolor('black')
for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
    for name in files:
        if name.endswith('.txt') == False:  #Check for correct file extension
            continue #Skips to next file if not .tif           

        #Building output folder
        file_path = os.path.join(path, name) #Get file path of current file
        print(file_path)
        df = pd.read_csv(file_path,sep='\t')
        x = df.Xraw
        y = df.Yraw
        bbox_props = dict(boxstyle="round", fc="w", ec="w")
        arrow_props = dict(arrowstyle="simple",color='white')
        plt.plot(x,y,lw=2)
        x0,y0 = x.iloc[0],y.iloc[0]
        plt.annotate(traj_n,xy=(x0,y0),xytext=(x0+30,y0+30),bbox=bbox_props,arrowprops=arrow_props, size = 4)
        traj_n += 1
plt.savefig('/home/baptiste/fig2.png',dpi=400)
plt.show()"""

"""input_folder = r'/home/baptiste/Documents/Pamela/ResultatsTXT'
traj_n = 1
plt.gca().invert_yaxis()
ax = plt.axes()
ax.set_facecolor('black')
for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
    for name in files:
        if name.endswith('.txt') == False:  #Check for correct file extension
            continue #Skips to next file if not .tif           

        #Building output folder
        file_path = os.path.join(path, name) #Get file path of current file
        print(file_path)
        df = pd.read_csv(file_path,sep='\t')
        x = df.Xraw
        y = df.Yraw
        bbox_props = dict(boxstyle="round", fc="w", ec="w")
        arrow_props = dict(arrowstyle="simple",color='white')
        plt.plot(x,y,lw=2)
        x0,y0 = x.iloc[0],y.iloc[0]
        plt.annotate(traj_n,xy=(x0,y0),xytext=(x0+30,y0+30),bbox=bbox_props,arrowprops=arrow_props, size = 4)
        traj_n += 1
plt.savefig('/home/baptiste/fig2.png',dpi=400)
plt.show()
"""

def final_image_ouput(processed_trajectory):
    """Plots all trajectories onto the first frame of a file. Inputs a dictionary, a string, a NumPy array, and a DataFrame.

                Functionally identical to image_output, with the exception of the sub trajectory for loop being included within the function itself.
                parameters is a dictionary containing the current output file path.
                name is the name of the file currently being processed.
                frames is a NumPy array containing raw frames from the video file.
                processed_trajectory is a DataFrame containing the trajectory to be plotted.
    """
    #Initializes plot
    plt.title('Trajectories from ')
    plt.gca().invert_yaxis() #Plotting the y axis inverts it by default, so it must be inverted again
    plt.xlim(0, 88.576)
    plt.ylim(0, 60.896)
    plt.gca().set_aspect('equal', adjustable='box')

    for item in set(processed_trajectory.rejoined_particle): #Loops for each sub trajectory
        x = (processed_trajectory[processed_trajectory.rejoined_particle==item].x)*0.173
        y = (processed_trajectory[processed_trajectory.rejoined_particle==item].y)*0.173
        plt.plot(x,y,lw=0.5)
    
    plt.show()

processed_trajectory = pd.read_csv(r'/media/baptiste/Windows/Users/LUMIN10/Documents/wfh/124 Results - 20211104_102617/124/Exp_2_20190226_kif5a_nKTP/HET/larve2/oeil_gauche/190226-nanoKTP-kif5a.lif - Series010.tif/190226-nanoKTP-kif5a.lif - Series010.tif_rejoined.csv',sep='\t')
final_image_ouput(processed_trajectory)