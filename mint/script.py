""" Script used to track particles in video files and extract trajectories.

    This Python script performs frame by frame analysis of video files, localizes particles and then reconstructs trajectories. 

    Please refer to https://github.com/biophotlumin/mint and https://lumin-mint.readthedocs.io/en/latest/ for more information.
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
    'threshold_r':40,
    #Data Extraction
    'r_conf_cut' : 0.64,
    'px' : 0.173, #in µm
    'dt' : 0.05, #in s
    'min_theoretical_precision' : 50, # in nm
    'sliding_window':3,
    'sigma':129,
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
    'individual_images':False,
    'individual_txt':False,
    'group_image':True,
    #Data Extraction
    'polynomial_fit':True,
    'minimization':True,
    'antero_retro':True,
}

log = {
    #Rejoining
    'number_rejoined':0,
}

##
##


#Define root input folder

input_folder = Path(r"/home/baptiste/Documents/vis")

start = time.time()

if __name__ == '__main__':
    
    #Output folder initialization
    log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(input_folder)
    print(log['output_folder'])
    print(Path(log['output_folder']).joinpath(input_folder.name))

    #Calling main functions
    tracking(input_folder,parameters,settings,log)

    data_extraction(parameters,Path(log['output_folder']).joinpath(Path(log['output_folder']).joinpath(input_folder.name)),settings)

    # statistical_analysis(settings,log['output_folder'])


end = time.time()

duration = end - start
print('%dh%s' % (int(duration//3600),f'{int((duration%3600)/60):02d}'))

