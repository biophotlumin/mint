"""Provides an interactive command line interface (CLI), as an alternative to the GUI or script.py

    The CLI essentially recapitulates the parameters and settings available in script.py.
    It should also raise an error if the parameters provided are not in the right format.
"""
import os
import time
from pathlib import Path

from utils import *
from image_denoising import *
from output_files_creation import *
from data_extraction import *
from statistical_analysis import *
from tracking import *

from prompt_toolkit.validation import Validator, ValidationError
from PyInquirer import prompt,style_from_dict, Token

class NumberValidator(Validator):
    def validate(self, document):
        try:
            float(document.text)
        except ValueError:
            raise ValidationError(
                message='Please enter a number',
                cursor_position=len(document.text))  # Move cursor to end

class FolderValidator(Validator):
    def validate(self, document):
        if os.path.exists(document.text)==True:
            pass
        else :
            raise ValidationError(
                message='Folder does not exist',
                cursor_position=len(document.text))  # Move cursor to end

if __name__ == '__main__':

    print('\n')
    print(' ███╗   ███╗██╗███╗   ██╗████████╗')
    print(' ████╗ ████║██║████╗  ██║╚══██╔══╝')
    print(' ██╔████╔██║██║██╔██╗ ██║   ██║')
    print(' ██║╚██╔╝██║██║██║╚██╗██║   ██║')
    print(' ██║ ╚═╝ ██║██║██║ ╚████║   ██║')
    print(' ╚═╝     ╚═╝╚═╝╚═╝  ╚═══╝   ╚═╝')
    print('\n')
    print(' ================================')
    print('              0.1.3')
    print(' ================================')
    print('\n')
    print(' ################################')
    print('         Please refer to')
    print('    lumin-mint.readthedocs.io')
    print('                or')
    print('    github.com/biophotlumin/mint')
    print('      for further information')
    print(' ################################')
    print('\n')
    questions = [
        {
            'type': 'input',
            'name': 'input_folder',
            'message': 'Select folder with input data',
            'default': '',
            'filter': lambda val: Path(val),
            'validate': FolderValidator
        },
        {
            'type': 'checkbox',
            'message': 'Preprocessing',
            'name': 'preprocess',
            'choices': [
                {
                    'name': 'Tophat'
                },
                {
                    'name': 'Wavelet'
                }
            ]},         
        {
            'type': 'input',
            'name': 'diameter',
            'message': '[trackpy] Feature diameter (in pixels)',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'minmass',
            'message': '[trackpy] Minimum integrated brightness',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'separation',
            'message': '[trackpy] Feature separation (in pixels)',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'search_range',
            'message': '[trackpy] Search range (in pixels)',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'memory',
            'message': '[trackpy] Memory (in number of frames)',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'adaptive_stop',
            'message': '[trackpy] Adaptive stop (in pixels)',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'adaptive_step',
            'message': '[trackpy] Adaptive step',
            'validate': NumberValidator,
            'filter': lambda val: float(val)
        },
        {
            'type': 'checkbox',
            'message': 'Tracking options',
            'name': 'options',
            'choices': [
                {
                    'name': 'Stub filtering'
                },
                {
                    'name': 'MSD filtering'
                },
                {
                    'name': 'SNR estimation'
                },
                {
                    'name': 'Rejoining'
                },
            ]},
        {
            'type': 'input',
            'name': 'stub_filtering',
            'message': '[trackpy] Stub filtering (in number of points)',
            'validate': NumberValidator,
            'filter': lambda val: int(val),
            'when': lambda answers: 'Stub filtering' in answers['options'],
            
        },
        {
            'type': 'input',
            'name': 'threshold',
            'message': 'MSD threshold',
            'validate': NumberValidator,
            'filter': lambda val: int(val),
            'when': lambda answers: 'MSD filtering' in answers['options'] 
        },
        {
            'type': 'input',
            'name': 'base_level',
            'message': 'Base level for SNR estimation',
            'validate': NumberValidator,
            'filter': lambda val: int(val),
            'when': lambda answers: 'SNR estimation' in answers['options']
        },
        {
            'type': 'input',
            'name': 'threshold_t',
            'message': '[Rejoining] Temporal threshold (in number of frames)',
            'validate': NumberValidator,
            'filter': lambda val: int(val),
            'when': lambda answers: 'Rejoining' in answers['options']
        },
        {
            'type': 'input',
            'name': 'threshold_r',
            'message': '[Rejoining] Spatial threshold (in pixels)',
            'validate': NumberValidator,
            'filter': lambda val: int(val),
            'when': lambda answers: 'Rejoining' in answers['options']
        },
        {
            'type': 'checkbox',
            'message': 'Output options',
            'name': 'output',
            'choices': [
                {
                    'name': 'Individual images'
                },
                {
                    'name': 'Individual text files'
                },
                {
                    'name': 'Group images'
                },
            ]},
        {
            'type': 'checkbox',
            'message': 'Data extraction options',
            'name': 'data_extraction',
            'choices': [
                {
                    'name': 'Polynomial fit'
                },
                {
                    'name': 'Convex minimization (trajectory denoising)'
                },
                {
                    'name': 'Anterograde / retrograde transport'
                },
            ]},
        {
            'type': 'input',
            'name': 'r_conf_cut',
            'message': 'Confinement ratio',
            'validate': NumberValidator,
            'filter': lambda val: float(val)
        },
        {
            'type': 'input',
            'name': 'px',
            'message': 'Pixel size (in µm)',
            'validate': NumberValidator,
            'filter': lambda val: float(val)
        },
        {
            'type': 'input',
            'name': 'dt',
            'message': 'Sampling rate (in seconds)',
            'validate': NumberValidator,
            'filter': lambda val: float(val)
        },
        {
            'type': 'input',
            'name': 'min_theoretical_precision',
            'message': 'Minimum theoretical precision (in nm)',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'sliding_window',
            'message': 'Sliding window',
            'validate': NumberValidator,
            'filter': lambda val: int(val)
        },
        {
            'type': 'input',
            'name': 'sigma',
            'message': 'Estimated noise (in nm)',
            'validate': NumberValidator,
            'filter': lambda val: int(val),
            'when': lambda answers: 'Convex minimization (trajectory denoising)' in answers['data_extraction'] 
        },
        {
            'type': 'input',
            'name': 'len_cutoff',
            'message': 'Cutoff length (in points)',
            'validate': NumberValidator,
            'filter': lambda val: int(val),
            'when': lambda answers: 'Polynomial fit' in answers['data_extraction']
        },
        {
            'type': 'input',
            'name': 'threshold_poly3',
            'message': 'Deviation from 3rd degree polynom',
            'validate': NumberValidator,
            'filter': lambda val: float(val),
            'when': lambda answers: 'Polynomial fit' in answers['data_extraction']
        },
    ]
    style = style_from_dict({
    Token.Separator: '#69e5ae',
    Token.QuestionMark: '#61C2A2 bold',
    Token.Selected: '#69e5ae',  # default
    Token.Pointer: '#61C2A2 bold',
    Token.Instruction: '',  # default
    Token.Answer: '#69e5ae bold',
    Token.Question: '',
    })

    answers = prompt(questions, style=style)
    parameters = {
        #trackpy.batch
        'diameter':answers['diameter'],
        'minmass':answers['minmass'],
        'separation':answers['separation'],
        #trackpy.link
        'search_range':answers['search_range'],
        'memory':answers['memory'],
        'adaptive_stop':answers['adaptive_stop'],
        'adaptive_step':answers['adaptive_step'],
        'stub_filtering':answers.get('stub_filtering'),
        #trackpy.motion.msd
        'threshold':answers.get('threshold'),
        #SNR estimation
        'base_level':answers.get('base_level'),
        #Rejoining
        'threshold_t':answers.get('threshold_t'),
        'threshold_r':answers.get('threshold_r'),
        #Data Extraction
        'r_conf_cut' : answers['r_conf_cut'],
        'px' : answers['px'], #in µm
        'dt' : answers['dt'], #in s
        'min_theoretical_precision' : answers['min_theoretical_precision'], # in nm
        'sliding_window':answers['sliding_window'],
        'sigma':answers.get('sigma'),
        'len_cutoff':answers.get('len_cutoff'), #Number of points
        'threshold_poly3':answers.get('threshold_poly3') #Deviation from third-degree polynom
    }

    settings = {
        #Denoising
        'tophat':'Tophat' in answers['preprocess'],
        'wavelet':'Wavelet' in answers['preprocess'],
        #Processing
        'stub_filtering':'Stub filtering' in answers['options'],
        'MSD':'MSD filtering' in answers['options'],
        'rejoining':'Rejoining' in answers['options'],
        'SNR_estimation':'MSD filtering' in answers['options'],
        #Outputs
        'individual_images':'Individual images' in answers['output'],
        'individual_txt':'Individual text files' in answers['output'],
        'group_image':'Group images' in answers['output'],
        #Data Extraction
        'polynomial_fit':'Polynomial fit' in answers['data_extraction'],
        'minimization':'Convex minimization (trajectory denoising)' in answers['data_extraction'],
        'antero_retro':'Anterograde / retrograde transport' in answers['data_extraction']
    }

    log = {
        #Rejoining
        'number_rejoined':0,
    }


    print('Starting up...')
    print('\n')

    input_folder = answers['input_folder']

    start = time.time()
    
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
