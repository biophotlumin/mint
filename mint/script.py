""" Script used to track particles in video files and extract trajectories.

    This Python script performs frame by frame analysis of video files, localizes particles and then reconstructs trajectories. 
    It can optionnally denoise individual frames, filter by mean square displacement (MSD), rejoin trajectories, and estimate signal to noise ratio.

    Parameters : 
        For particle localization (trackpy.batch) : 

            'diameter': Feature’s extent in each dimension.
            'minmass': Minimum integrated brightness.
            'separation': Minimum separation between features.

        For trajectory reconstitution (trackpy.link) : 
            'search_range': Maximum distance features can move between frames.
            'memory': Maximum number of frames during which a feature can vanish, then reappear nearby, and be considered the same particle.
            'adaptive_stop': When encountering an oversize subnet, retry by progressively reducing search_range until the subnet is solvable.
            'adaptive_step': Reduce search_range by multiplying it by this factor.

        For MSD filtering : 
            'threshold': Level against which the computed Mean Square Displacement is compared.

        For SNR estimation : 
            'BaseLevel': Base level of the detector.

    Settings : 
        'tophat' : Applies white top-hat transform to each individual frame.
        'wavelet' : Applies wavelet denoising to each individual frame.
        'MSD' : Filters trajectories based on computed Mean Square Displacement.
        'rejoining' : Rejoins trajectories that appear to be disjointed.
        'SNR_estimation' : Adds columns to the DataFrame containing the signal to noise ratio, the volume of the 2D Gaussian as well as its feet.

    input_folder is the root folder containing video files to be analyzed. For now, only .tif files are supported.
    The output is a clone of the root folder, identified by the time and date of the start of the run. Each .tif file creates a matching folder.
    Inside each file-specific folder is a .csv file containing the raw coordinates of trajectories found in that file. 
    Optionally, this folder will also contain a .csv file containing the coordinates of rejoined trajectories.
    For each trajectory, this script will optionally output a .txt file containing the coordinates of each particle, as well as an image of that trajectory plotted on the first frame of the video.
    It will also optionally output an image of all found trajectories plotted on the first frame of the video.
    If MSD filtering is enabled and the threshold is too high, no viable trajectories will be found and only the first .csv file will be created. 
"""
#Imports modules
from utils import *
from image_denoising import *
from output_files_creation import *
from data_extraction import *
from statistical_analysis import *
from tracking import *
import time

#Calculation parameters

parameters = {
    #trackpy.batch
    'diameter':9,
    'minmass':300,
    'separation':12,
    #trackpy.link
    'search_range':6,
    'memory':5,
    'adaptive_stop':5,
    'adaptive_step':0.9,
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
    'px' : 0.175, #in µm
    'dt' : 0.05, #in s
    'min_theoretical_precision' : 30, # in nm
    'sliding_window':3,
    'sigma':175,
    'len_cutoff':30, #Number of points
    'threshold_poly3':1.4 #Deviation from third-degree polynom
}

#Optional image and trajectory processing

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

log = {
    #Rejoining
    'number_rejoined':0,
}

##
##


#Define root input folder

input_folder = Path(r"F:\INSERM")

start = time.time()

if __name__ == '__main__':
    #Output folder initialization
    log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(input_folder)
    print(log['output_folder'])
    print(Path(log['output_folder']).joinpath(input_folder.name))
    #Calling main functions
    tracking(input_folder,parameters,settings,log)

    data_extraction(parameters,Path(log['output_folder']).joinpath(Path(log['output_folder']).joinpath(input_folder.name)),settings)

    statistical_analysis(settings,log['output_folder'])


end = time.time()
print((end-start)/60)

