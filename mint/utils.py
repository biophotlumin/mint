"""Module regrouping several small utilitary functions.
    extraction_csv automatically extracts .csv file from an output folder containing multiple .txt and .png files.
    folder_structure_creation returns the time and date of the start of the run, as well as the path to the output and root input folder.
"""
#Imports
import os
import shutil
from datetime import datetime
from pathlib import Path

def extraction_csv(input_folder):
    """Automatically extracts .csv file from an output folder containing multiple .txt and .png files.
        """
    print(input_folder)
    shutil.copytree(input_folder, input_folder+" CSV ONLY",symlinks=False,ignore=shutil.ignore_patterns('*.png','*.txt'))

def folder_structure_creation(input_folder):
    """Returns the time and date of the start of the run, as well as the path to the output and root input folder, all as strings.

            input_folder is the root folder containing video files to be analyzed.
            identifier is a string containing the suffix 'Results' as well as the time and date of the start of the run. 
            It is added to output file paths so that the results are stored in a folder different than root.
            output_folder is the input folder with the identifier added, and is thus the root ouput folder.
            root_input_folder is the folder above input_folder in the folder structure. It is used to generate proper output file paths later on.
                """               
    identifier = " Results - "+str(datetime.now().strftime('%Y%m%d_%H%M%S')) #Create a string that is added to output file paths 

    output_folder = Path(input_folder.parent).joinpath(input_folder.name + identifier) #Set path for root output folder
    root_input_folder = os.path.dirname(input_folder) #Get folder directly before root input folder

    if os.path.dirname(root_input_folder) == root_input_folder: #Prevents conflict in case input folder is placed at the root of a drive
        root_input_folder = input_folder

    if input_folder =='':
        input_folder = root_input_folder

    return output_folder,identifier,root_input_folder