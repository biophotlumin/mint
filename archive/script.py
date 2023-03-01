""" Script used to track particles in videos and extract trajectories.

    This Python script performs frame by frame analysis of videos, localizes particles and then reconstructs trajectories. 

    Please refer to https://github.com/biophotlumin/mint and https://lumin-mint.readthedocs.io/en/latest/ for more information.
"""

#Imports modules
from utils import *
from image_denoising import *
from output import *
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
    'stub_filtering':3, # Minimum length of trajectories, in points
    #trackpy.motion.msd
    'msd':300, # Threshold for MSD filtering
    #SNR estimation
    'base_level':0, #Base level for SNR estimation
    #Rejoining
    'threshold_t':10, # Temporal threshold for trajectory rejoining, in frames
    'threshold_r':40, # Spatial threshold for trajectory rejoining, in pixels
    #Data Extraction
    'r_conf_cut' : 0.64, # Confinement ratio cutoff between GO and STOP phases
    'px' : 0.173, # Pixel size, in µm
    'dt' : 0.05, # Time interaval between frames, in s
    'min_theoretical_precision' : 50, # Minimum theoretical precision, in nm
    'sliding_window':3, # Sliding windows for confinement ratio calculation
    'sigma':129, # Estimated precision of localization
    'len_cutoff':30, # Number of points
    'threshold_poly3':1.4, # Deviation from third-degree polynom
    #Stats
    'order':['WT','HET','HOM'], # Order of conditions in tables and graph
    'extension':'svg', # File extension of saved graphs
    'dpi':300, # DPI of saved graphs if they're not vectorial
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
    'individual_images':False,
    'individual_txt':False,
    'group_image':True,
    #Data Extraction
    'polynomial_fit':True,
    'minimization':True,
    'theta':True,
    'antero_retro':True,
    #Stats
    'ordering':True,
    'clean_up':False,
}

##
##

#Define root input folder

input_folder = Path(r"/media/baptiste/SHG_tracking_data/Zebrafish data/video_benchmark")

if __name__ == '__main__':

    start = time.time()
    
    #Output folder initialization
    log = {}
    log['input_folder'] = input_folder
    log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(input_folder)
    print(f'\nAnalyzing {input_folder}\n')
    print(f'Results stored in {Path(log["output_folder"].name)}')

    #Calling main functions
    tracking(input_folder,parameters,settings,log)

    data_extraction(Path(log['output_folder']).joinpath(input_folder.name),parameters,settings)

    statistical_analysis(settings,parameters,log['output_folder'])

    end = time.time()

    duration = end - start
    f_duration = '%dh%s' % (int(duration//3600),f'{int((duration%3600)/60):02d}')
    log['duration'] = f_duration
    print('Total runtime : ' )
    print(f_duration)

    #Dumps parameters and settings into .txt files at the root of the output folder
    dict_dump(log['output_folder'],parameters,'parameters')
    dict_dump(log['output_folder'],settings,'settings')
    dict_dump(log['output_folder'],log,'log')

    generate_report(log['output_folder'])

# TODO Declare arg types in functions