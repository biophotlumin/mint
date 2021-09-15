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
    'tophat':True,
    'wavelet':True,
    'MSD':True,
    'rejoining':True,
    'SNR_estimation':True
    }

sg.theme('DarkBlue')   

tab_trackpy_layout =  [          
            
            [sg.Text('Input folder', size=(8, 1)), sg.Input(default_text=default_dict['input_folder'],key='input_folder'), sg.FolderBrowse()],
            [sg.Frame(layout=[
            [sg.Text('Diameter'), sg.InputText(default_text=default_dict['diameter'],key='diameter', size=(5,1)),sg.Button('Test Locate')],
            [sg.Text('Minmass'), sg.InputText(default_text=default_dict['minmass'],key='minmass', size=(5,1)),\
                sg.Radio('First image','test_locate',key='first_image',enable_events=True)],
            [sg.Text('Separation'), sg.InputText(default_text=default_dict['separation'],key='separation', size=(5,1)),\
                sg.Radio('Whole file','test_locate',key='whole_file',enable_events=True)],
            ],\
                title='Locate',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Frame(layout=[
            [sg.Text('Search Range'), sg.InputText(default_text=default_dict['search_range'],key='search_range', size=(5,1))],
            [sg.Text('Memory         '), sg.InputText(default_text=default_dict['memory'],key='memory', size=(5,1))],
            [sg.Text('Adaptative Stop'), sg.InputText(default_text=default_dict['adaptive_stop'],key='adaptive_stop', size=(5,1))],
            [sg.Text('Adaptative Step'), sg.InputText(default_text=default_dict['adaptive_step'],key='adaptive_step', size=(5,1))]],\
                title='Link',title_color='white',relief = sg.RELIEF_SUNKEN)],
            [sg.Button('Run TrackPy')]]



tab_settings = [            
            #[sg.Frame(layout=[      
            [sg.Checkbox('Tophat', size=(10,1), key='tophat',enable_events=True, default=default_dict['tophat'])],
            [sg.Checkbox('Wavelet', size=(10,1), key='wavelet',enable_events=True,default=default_dict['wavelet'])],
            [sg.Checkbox('MSD', size=(10,1), key='MSD',enable_events=True,default=default_dict['MSD']),\
                sg.Text('Threshold'), sg.InputText(default_text=default_dict['threshold'],disabled = not default_dict['MSD'],key='threshold', size=(5,1))],
            [sg.Checkbox('Rejoining', size=(10,1), key='rejoining',enable_events=True, default=default_dict['rejoining']),\
                sg.Text('Temporal Threshold'),sg.InputText(default_text=default_dict['threshold_t'],key='threshold_t', size=(5,1)),\
                    sg.Text('Spatial Threshold'),sg.InputText(default_text=default_dict['threshold_r'],key='threshold_r', size=(5,1))],
            [sg.Checkbox('SNR estimation', size=(11,1), key='SNR_estimation',enable_events=True,default=default_dict['SNR_estimation']),\
                sg.Text('Base Level'), sg.InputText(default_text=default_dict['BaseLevel'],key='BaseLevel', disabled=not default_dict['SNR_estimation'],\
                     size=(5,1))]]
                         #title='Options',title_color='white', relief=sg.RELIEF_SUNKEN, tooltip='Use these to set flags')]]  
                         # 
tab_extraction = [[sg.Button('Run')]]
tab_stats = [[sg.Button('Run')]] 

layout = [  
        [sg.TabGroup([[sg.Tab('TrackPy', tab_trackpy_layout), sg.Tab('Settings', tab_settings),sg.Tab('Data Extraction', tab_extraction),sg.Tab('Statistics',tab_stats)]])],    
        [sg.Button('Run'), sg.Button('Exit')] ]

# Create the Window
window = sg.Window('TrackPy', layout)
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
    #trackpy.motion.msd
    'threshold':int(values['threshold']),
    #SNR estimation
    'BaseLevel':int(values['BaseLevel']),
    #rejoining
    'threshold_t':int(values['threshold_t']),
    'threshold_r':int(values['threshold_r'])
    }

    settings = {
    'tophat':values['tophat'],
    'wavelet':values['wavelet'],
    'MSD':values['MSD'],
    'rejoining':values['rejoining'],
    'SNR_estimation':values['SNR_estimation']
    }

    log = {
    #Rejoining
    'number_rejoined':0,
}
    if values['MSD']==True or values['MSD']==False:
        window['threshold'].update(disabled=(not values['MSD'])) 

    if values['SNR_estimation']==True or values['SNR_estimation']==False:
        window['BaseLevel'].update(disabled=(not values['SNR_estimation']))    
        
    if event == 'Run':
        log['output_folder'],log['identifier'],log['root_input_folder'] = folder_structure_creation(values['input_folder'])
        start = time.time()

        tracking(values['input_folder'],parameters,settings,log)
        data_extraction(log['output_folder']+"\\"+values['input_folder'].replace(log['root_input_folder'],''))
        statistical_analysis(log['output_folder'])

        end = time.time()
        print((end-start)/60)

        parameters['input_folder'] = values['input_folder']
        default_dict = {**parameters, **settings}
        with open("default.ini", 'w') as default_ini:
            print(default_dict, file=default_ini)

    if event == 'Run TrackPy':
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