""" Script used to track particles in videos and extract trajectories.

    This Python script performs automated analysis of videos, including tracking,
    data analysis and statistical tests.

    Please refer to https: //github.com/biophotlumin/mint and
    https: //lumin-mint.readthedocs.io/en/latest/ for more information.
"""

#Imports modules
import os
import sys
import shutil
import warnings
import time
import argparse

from pathlib import Path
from .tracking import tracking, p_tracking
from .data_extraction import data_extraction
from .utils import  folder_structure_creation, load_params
from .output import dict_dump, generate_report
from .stat_analysis import statistical_analysis

def main():

    #Calculation parameters

    parameters = {
        'extension_in': 'tif',
        # trackpy.batch
        'diameter': 9,
        'minmass': 300,
        'separation': 12,
        # trackpy.link
        'search_range': 6,
        'memory': 5,
        'adaptive_stop': 5,
        'adaptive_step': 0.9,
        'stub_filtering': 3, # Minimum length of trajectories, in points
        # trackpy.motion.msd
        'msd': 9, # Threshold for MSD filtering
        # SNR estimation
        'base_level': 0, # Base level for SNR estimation
        # Rejoining
        'threshold_t': 10, # Temporal threshold for trajectory rejoining, in frames
        'threshold_r': 40, # Spatial threshold for trajectory rejoining, in pixels
        # Data Extraction
        'r_conf_cut':  0.64, # Confinement ratio cutoff between GO and STOP phases
        'px':  0.173, # Pixel size, in µm
        'dt':  0.05, # Time interaval between frames, in s
        'min_thr_prec':  50, # Minimum theoretical precision, in nm
        'sliding_window': 3, # Sliding window for confinement ratio calculation
        'sigma': 129, # Estimated precision of localization, in nm
        'len_cutoff': 30, # Number of points
        'threshold_poly3': 1.4, # Deviation from third-degree polynom
        # Stats
        'order': ['WT', 'HET', 'HOM'], # Order of conditions in tables and graphs
        'extension_out': 'svg', # File extension of saved graphs
        'dpi': 300, # DPI of saved graphs if they're not vectorial
    }

    # Optional image and trajectory processing

    settings = {
        'parallel': True,
        # Denoising
        'tophat': True,
        'wavelet': False,
        # Processing
        'stub_filtering': False,
        'MSD': True,
        'rejoining': True,
        'SNR_estimation': False,
        # Outputs
        'individual_images': False,
        'individual_txt': False,
        'group_image': True,
        # Data Extraction
        'polynomial_fit': True,
        'minimization': True,
        'theta': True,
        'antero_retro': True,
        'conf_list': False,
        # Stats
        'ordering': True,
        'clean_up': False,
        #
        'parallel_tracking': False,
    }

    ##
    ##

    #Define root input folder

    input_folder = r''

    parser = argparse.ArgumentParser(prog='M.I.N.T',
                                     description='Intraneuronal nanoparticle tracking')
    parser.add_argument('-f', '--folder',
                        default=os.getcwd(),
                        help='Path to data folder')
    parser.add_argument('-p', '--params',
                        default='', help='Path to config file')
    parser.add_argument('-l', '--locate',
                        default=False, help='Locate and link particles',
                        action='store_true')
    parser.add_argument('-e', '--extract',
                        default=False, help='Extract transport parameters',
                        action='store_true')
    parser.add_argument('-s', '--stats',
                        default=False, help='Statistical analysis',
                        action='store_true')

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    print('\n')
    print(' ######################################################################### ')
    print(' ######################################################################### ')
    print('\n')
    print('                    ███╗   ███╗   ██╗   ███╗  ██╗████████╗                 ')
    print('                    ████╗ ████║   ██║   ████╗ ██║╚══██╔══╝                 ')
    print('                    ██╔████╔██║   ██║   ██╔██╗██║   ██║                    ')
    print('                    ██║╚██╔╝██║   ██║   ██║╚████║   ██║                    ')
    print('                    ██║ ╚═╝ ██║██╗██║██╗██║ ╚███║██╗██║                    ')
    print('                    ╚═╝     ╚═╝╚═╝╚═╝╚═╝╚═╝  ╚══╝╚═╝╚═╝                    ')
    print('\n')
    print('                                   v0.3.0                                  ')
    print('\n')
    print(' ######################################################################### ')
    print(' ######################################################################### ')

    if args.params:
        config = load_params(args.params)
        input_folder = config.get('input_folder', os.getcwd())
        settings = config['settings']
        parameters = config['parameters']
        if args.folder != os.getcwd():
            input_folder = args.folder
    else:
        if args.folder:
            input_folder = args.folder
            config = load_params(Path(input_folder).joinpath('config.yml'))
            if config is not None:
                settings = config['settings']
                parameters = config['parameters']
            else:
                warnings.warn('Running with default parameters')

    input_folder = Path(input_folder)

    start = time.time()

    #Output folder initialization
    log = {}
    log['input_folder'] = input_folder
    (log['output_folder'],
     log['identifier'],
     log['root_input_folder']
     ) = folder_structure_creation(input_folder)

    os.makedirs(log['output_folder'])

    dict_dump(log['output_folder'], parameters, 'parameters')
    dict_dump(log['output_folder'], settings, 'settings')
    dict_dump(log['output_folder'], log, 'log')

    print(f'\nAnalyzing {input_folder.name}\n')
    print(f'Results stored in {Path(log["output_folder"].name)}')

    if settings['parallel_tracking']:
        tracker = p_tracking
    else:
        tracker = tracking

    if args.locate:
        tracker(input_folder, parameters, settings, log)
    if args.extract:
        if args.locate:
            data_extraction(Path(log['output_folder']).joinpath(input_folder.name),
                            parameters, settings)
        else:
            data_extraction(input_folder, parameters, settings)
    if args.stats:
        if args.extract:
            statistical_analysis(settings, parameters, Path(log['output_folder']))
        else:
            # Rerunning stats overwrites previous results for now
            shutil.rmtree(log['output_folder'])
            log['output_folder'] = input_folder
            statistical_analysis(settings, parameters, input_folder)

    if True not in vars(args).values():
        tracker(input_folder, parameters, settings, log)
        data_extraction(Path(log['output_folder']).joinpath(input_folder.name),
                        parameters, settings)
        statistical_analysis(settings, parameters, log['output_folder'])

    end = time.time()

    duration = end - start
    log['duration_s'] = duration
    f_duration = f'{int(duration//3600)}h{int((duration%3600)/60):02d}'
    log['duration'] = f_duration
    print(f'Total runtime : {f_duration}')

    dict_dump(log['output_folder'], log, 'log')

    if True not in vars(args).values():
        generate_report(log['output_folder'])

if __name__ == '__main__':
    main()

# TODO GUI
# TODO osqp and qdldl wheels not existing for Python-3.12