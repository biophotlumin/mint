import PySimpleGUI as sg
from utils import folder_structure_creation
from image_denoising import *
from output_files_creation import *
from tracking import *
from data_extraction import *
from statistical_analysis import *
from test_locate import *
import time
import ast
import os

if os.path.isfile('default.ini')==True:
    with open('default.ini', 'r') as default:
        default_dict = ast.literal_eval(default.read())
else:
    default_dict = {
    'input_folder':'',
    'diameter':int(),
    'minmass':int(),
    'separation':int(),
    'search_range':int(),
    'memory':int(),
    'adaptive_stop':int(),
    'adaptive_step':float(),
    'threshold':int(),
    'BaseLevel':int(),
    'threshold_t':int(),
    'threshold_r':int(),
    'stub_length':int(),
    'r_conf_cut':float(),
    'px':float(),
    'dt':float(),
    'min_theoretical_precision':int(),
    'sliding_window':int(),
    'sigma':int(),
    'len_cutoff':int(),
    'threshold_poly3':float(),
    'tophat':True,
    'wavelet':True,
    'MSD':True,
    'rejoining':True,
    'SNR_estimation':True,
    'stub_filtering':True,
    'individual_images':True,
    'individual_txt':True,
    'group_image':True,
    'polynomial_fit':True,
    'minimization':True,
    'antero_retro':True
    }

sg.theme('DarkBlue')   

tab_trackpy_layout =  [          
            
            [sg.Text('Input folder', size=(8, 1)), sg.Input(default_text=default_dict['input_folder'],key='input_folder'), sg.FolderBrowse()],
            [sg.Frame(layout=[
            [sg.Text('Diameter'), sg.InputText(default_text=default_dict['diameter'],key='diameter', size=(5,1))],
            [sg.Text('Minmass'), sg.InputText(default_text=default_dict['minmass'],key='minmass', size=(5,1))],
            [sg.Text('Separation'), sg.InputText(default_text=default_dict['separation'],key='separation', size=(5,1))],
            ],\
                title='Locate',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Frame(layout=[
            [sg.Button('Test Locate')],
            [sg.Radio('First image','test_locate',key='first_image',enable_events=True)],
            [sg.Radio('Whole file','test_locate',key='whole_file',enable_events=True)],
            ],\
                title='Locate',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Frame(layout=[
            [sg.Text('Search Range'), sg.InputText(default_text=default_dict['search_range'],key='search_range', size=(5,1))],
            [sg.Text('Memory         '), sg.InputText(default_text=default_dict['memory'],key='memory', size=(5,1))],
            [sg.Text('Adaptative Stop'), sg.InputText(default_text=default_dict['adaptive_stop'],key='adaptive_stop', size=(5,1))],
            [sg.Text('Adaptative Step'), sg.InputText(default_text=default_dict['adaptive_step'],key='adaptive_step', size=(5,1))]],\
                title='Link',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Button('Run tracking')]]



tab_settings = [            
            [sg.Frame(layout=[      
            [sg.Checkbox('Tophat', size=(10,1), key='tophat',enable_events=True, default=default_dict['tophat'])],
            [sg.Checkbox('Wavelet', size=(10,1), key='wavelet',enable_events=True,default=default_dict['wavelet'])]],\
                title='Denoising',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Frame(layout=[
            [sg.Checkbox('MSD', size=(10,1), key='MSD',enable_events=True,default=default_dict['MSD']),\
                sg.Text('Threshold'), sg.InputText(default_text=default_dict['threshold'],disabled = not default_dict['MSD'],key='threshold', size=(5,1))]],\
                    title='',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Frame(layout=[
            [sg.Checkbox('Rejoining', size=(10,1), key='rejoining',enable_events=True, default=default_dict['rejoining']),\
                sg.Text('Temporal Threshold'),sg.InputText(default_text=default_dict['threshold_t'],key='threshold_t', size=(5,1)),\
                    sg.Text('Spatial Threshold'),sg.InputText(default_text=default_dict['threshold_r'],key='threshold_r', size=(5,1))]],\
                        title='',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Frame(layout=[  
            [sg.Checkbox('Stub filtering', size=(11,1), key='stub_filtering',enable_events=True,default=default_dict['stub_filtering']),\
                sg.Text('Stub length'), sg.InputText(default_text=default_dict['stub_length'],key='stub_length', disabled=not default_dict['stub_filtering'],size=(5,1))]],\
                    title='',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Frame(layout=[  
            [sg.Checkbox('SNR estimation', size=(11,1), key='SNR_estimation',enable_events=True,default=default_dict['SNR_estimation']),\
                sg.Text('Base Level'), sg.InputText(default_text=default_dict['BaseLevel'],key='BaseLevel', disabled=not default_dict['SNR_estimation'],size=(5,1))]],\
                    title='',title_color='white',relief = sg.RELIEF_SUNKEN)]]
                         
tab_extraction = [
            [sg.Frame(layout=[  
                [sg.Text('Pixel size (in µm)'),sg.InputText(default_text=default_dict['px'],key='px', size=(5,1))],
                [sg.Text('Sampling time (in s)'),sg.InputText(default_text=default_dict['dt'],key='dt', size=(5,1))]],
                    title='Video parameters',title_color='white',relief = sg.RELIEF_SUNKEN)],

            [sg.Frame(layout=[
                [sg.Text('Cutoff'),sg.InputText(default_text=default_dict['r_conf_cut'],key='r_conf_cut', size=(5,1))],  
                [sg.Text('Minimum theoretical precision (in nm)'),sg.InputText(default_text=default_dict['min_theoretical_precision'],key='min_theoretical_precision', size=(5,1))],
                [sg.Text('Sliding window'),sg.InputText(default_text=default_dict['sliding_window'],key='sliding_window', size=(5,1))]],
                    title='Confinement ratio',title_color='white',relief = sg.RELIEF_SUNKEN)],

            [sg.Frame(layout=[  
                [sg.Checkbox('Anterograde / retrograde', size=(20,1), key='antero_retro',enable_events=True,default=default_dict['antero_retro'])]],
                    title='Anterograde / retrograde',title_color='white',relief = sg.RELIEF_SUNKEN)],

            [sg.Frame(layout=[  
                [sg.Checkbox('Polynomial fit', size=(11,1), key='polynomial_fit',enable_events=True,default=default_dict['polynomial_fit']),\
                    sg.Text('Deviation'), sg.InputText(default_text=default_dict['threshold_poly3'],key='threshold_poly3', disabled=not default_dict['polynomial_fit'],size=(5,1)),\
                        sg.Text('Length cutoff'), sg.InputText(default_text=default_dict['len_cutoff'],key='len_cutoff', disabled=not default_dict['polynomial_fit'],size=(5,1))]],
                    title='Polynomial fit',title_color='white',relief = sg.RELIEF_SUNKEN)],
            
            [sg.Frame(layout=[  
                [sg.Checkbox('Minimization', size=(10,1), key='minimization',enable_events=True,default=default_dict['minimization']),\
                    sg.Text('Sigma'), sg.InputText(default_text=default_dict['sigma'],key='sigma', disabled=not default_dict['minimization'],size=(5,1))]],
                    title='Minimization',title_color='white',relief = sg.RELIEF_SUNKEN)],

            [sg.Button('Run data extraction')]
            ]


tab_stats = [
            [sg.Frame(layout=[  
                [sg.Checkbox('Individual images', size=(20,1), key='individual_images',enable_events=True,default=default_dict['individual_images'])],
                    [sg.Checkbox('Individual text files', size=(20,1), key='individual_txt',enable_events=True,default=default_dict['individual_txt'])],
                        [sg.Checkbox('Group images', size=(20,1), key='group_image',enable_events=True,default=default_dict['group_image'])]],
                    title='Minimization',title_color='white',relief = sg.RELIEF_SUNKEN)],
    
            [sg.Button('Run statistical analysis')]] 

layout = [  
        [sg.TabGroup([[sg.Tab('TrackPy', tab_trackpy_layout), sg.Tab('Settings', tab_settings),sg.Tab('Data Extraction', tab_extraction),sg.Tab('Output and stats',tab_stats)]])],    
        [sg.Button('Run'), sg.Button('Exit')] ]

# Create the Window
window = sg.Window('MINT', layout)
# Event Loop to process "events" and get the "values" of the inputs

while True:
    
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Exit': # if user closes window or clicks cancel
        break

    parameters = {
    #trackpy.batch
    'diameter':int(values['diameter']),
    'minmass':int(values['minmass']),
    'separation':int(values['separation']),
    #trackpy.link
    'search_range':int(values['search_range']),
    'memory':int(values['memory']),
    'adaptive_stop':int(values['adaptive_stop']),
    'adaptive_step':float(values['adaptive_step']),
    'stub_filtering':int(values['stub_length']),
    #trackpy.motion.msd
    'threshold':int(values['threshold']),
    #SNR estimation
    'BaseLevel':int(values['BaseLevel']),
    #rejoining
    'threshold_t':int(values['threshold_t']),
    'threshold_r':int(values['threshold_r']),
    #Data Extraction
    'r_conf_cut':float(values['r_conf_cut']),
    'px':float(values['px']),
    'dt':float(values['dt']),
    'min_theoretical_precision':int(values['min_theoretical_precision']),
    'sliding_window':int(values['sliding_window']),
    'sigma':int(values['sigma']),
    'len_cutoff':int(values['len_cutoff']),
    'threshold_poly3':float(values['threshold_poly3']),
    }

    settings = {
    'tophat':values['tophat'],
    'wavelet':values['wavelet'],
    'stub_filtering':values['stub_filtering'],
    'MSD':values['MSD'],
    'rejoining':values['rejoining'],
    'SNR_estimation':values['SNR_estimation'],
    'individual_images':values['individual_images'],
    'individual_txt':values['individual_txt'],
    'group_image':values['group_image'],
    'polynomial_fit':values['polynomial_fit'],
    'minimization':values['minimization'],
    'antero_retro':values['antero_retro']
    }

    log = {
    #Rejoining
    'number_rejoined':0,
}
    if values['MSD']==True or values['MSD']==False:
        window['threshold'].update(disabled=(not values['MSD'])) 
    
    if values['rejoining']==True or values['rejoining']==False:
        window['threshold_t'].update(disabled=(not values['rejoining'])) 
        window['threshold_r'].update(disabled=(not values['rejoining'])) 

    if values['stub_filtering']==True or values['stub_filtering']==False:
        window['stub_length'].update(disabled=(not values['stub_filtering']))

    if values['SNR_estimation']==True or values['SNR_estimation']==False:
        window['BaseLevel'].update(disabled=(not values['SNR_estimation']))

    if values['polynomial_fit']==True or values['polynomial_fit']==False:
        window['len_cutoff'].update(disabled=(not values['polynomial_fit'])) 
        window['threshold_poly3'].update(disabled=(not values['polynomial_fit']))  

    if values['minimization']==True or values['minimization']==False:
        window['sigma'].update(disabled=(not values['minimization']))   
        
    if event == 'Run':
        log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(values['input_folder'])
        start = time.time()

        tracking(values['input_folder'],parameters,settings,log)
        data_extraction(parameters,Path(log['output_folder']).joinpath(Path(log['output_folder']).joinpath(input_folder.name)),settings)
        statistical_analysis(settings,log['output_folder'])

        end = time.time()
        print((end-start)/60)

        parameters['input_folder'] = values['input_folder']
        default_dict = {**parameters, **settings}
        with open("default.ini", 'w') as default_ini:
            print(default_dict, file=default_ini)

    if event == 'Run tracking':
        log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(values['input_folder'])
        start = time.time()
        tracking(values['input_folder'],parameters,settings,log)
        end = time.time()
        print((end-start)/60)

        parameters['input_folder'] = values['input_folder']
        default_dict = {**parameters, **settings}
        with open("default.ini", 'w') as default_ini:
            print(default_dict, file=default_ini)
    
    if event == 'Test Locate':
        if values['whole_file'] == True:
            whole = True
        else:
            whole = False
        test_locate(values['input_folder'],parameters,whole,settings)
    

window.close()