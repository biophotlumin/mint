""" Script used to track particles in videos and extract trajectories.

    This Python script performs frame by frame analysis of videos, localizes particles and then reconstructs trajectories. 

    Please refer to https://github.com/biophotlumin/mint and https://lumin-mint.readthedocs.io/en/latest/ for more information.
"""

#Imports modules
import os
import yaml
import time
import argparse

from pathlib import Path
from tracking import tracking
from data_extraction import data_extraction
from utils import  folder_structure_creation
from output import dict_dump, generate_report
from statistical_analysis import statistical_analysis

#Calculation parameters

parameters = {
    'extension_in':'tif',
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
    'msd':9, # Threshold for MSD filtering
    #SNR estimation
    'base_level':0, #Base level for SNR estimation
    #Rejoining
    'threshold_t':10, # Temporal threshold for trajectory rejoining, in frames
    'threshold_r':40, # Spatial threshold for trajectory rejoining, in pixels
    #Data Extraction
    'r_conf_cut' : 0.64, # Confinement ratio cutoff between GO and STOP phases
    'px' : 0.173, # Pixel size, in µm
    'dt' : 0.05, # Time interaval between frames, in s
    'min_thr_prec' : 50, # Minimum theoretical precision, in nm
    'sliding_window':3, # Sliding windows for confinement ratio calculation
    'sigma':129, # Estimated precision of localization
    'len_cutoff':30, # Number of points
    'threshold_poly3':1.4, # Deviation from third-degree polynom
    #Stats
    'order':['WT','HET','HOM'], # Order of conditions in tables and graph
    'extension_out':'svg', # File extension of saved graphs
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

input_folder = r'/media/baptiste/Windows/Users/LUMIN10/Documents/video_benchmark_int'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='M.I.N.T',description='Intraneuronal nanoparticle tracking')
    parser.add_argument('-f','--folder',default=input_folder if input_folder else os.getcwd(),help='Path to data folder')
    parser.add_argument('-p', '--params',default='',help='Path to config file')
    args = parser.parse_args()

    if args.params:
        config = yaml.safe_load(open(args.params))
        parameters = config['parameters']
        settings = config['settings']
    
    if args.folder:
        input_folder = args.folder

    input_folder = Path(input_folder)

    start = time.time()
    
    #Output folder initialization
    log = {}
    log['input_folder'] = input_folder
    log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(input_folder)
    os.makedirs(log['output_folder'])

    dict_dump(log['output_folder'],parameters,'parameters')
    dict_dump(log['output_folder'],settings,'settings')
    dict_dump(log['output_folder'],log,'log')

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
    print(f'Total runtime : {f_duration}' )

    dict_dump(log['output_folder'],log,'log')

    generate_report(log['output_folder'])

# TODO Iter over files then iter over paths