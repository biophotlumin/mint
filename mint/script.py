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
    'stub_filtering':3,
    #trackpy.motion.msd
    'msd':300,
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
    'threshold_poly3':1.4, #Deviation from third-degree polynom
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
    'order':['WT','HET','HOM'],
    'extension':'svg',
    'dpi':300,
}

log = {
    #Rejoining
    'number_rejoined':0,
}

##
##


#Define root input folder

input_folder = Path(r"/media/baptiste/SHG_tracking_data/Zebrafish data/video_benchmark")

start = time.time()

if __name__ == '__main__':
    
    #Output folder initialization
    log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(input_folder)
    print(f'\nAnalyzing {input_folder}\n')
    print(f'Results stored in {Path(log["output_folder"].name)}')

    #Calling main functions
    tracking(input_folder,parameters,settings,log)

    data_extraction(Path(log['output_folder']).joinpath(input_folder.name),parameters,settings)

    # statistical_analysis(settings,log['output_folder'])


end = time.time()

duration = end - start
print('Total runtime : ')
print('%dh%s' % (int(duration//3600),f'{int((duration%3600)/60):02d}'))

# TODO Declare arg types in functions
# TODO Generate report through ReportLab
# TODO Read parameters and settings dict from a JSON file to use in command line