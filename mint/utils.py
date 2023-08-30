"""Module regrouping several small utility functions.
"""

#Imports
import os
import csv
import nd2
import shutil
import imageio
import warnings

from datetime import datetime
from pathlib import Path


def extraction_csv(input_folder):
    """Extracts .csv files from a folder into another.

    :param input_folder: Path to the input folder.
    :type input_folder: string or Path
    """

    print(input_folder)
    shutil.copytree(input_folder, input_folder+" CSV ONLY",symlinks=False,ignore=shutil.ignore_patterns('*.png','*.txt'))

def folder_structure_creation(input_folder):
    """Generates paths to output folders. 

    Adds "Results" and the time and date of the start of the run to the input folder.

    :param input_folder: Path to the input folder.
    :type input_folder: string or Path
    :return: `output_folder` is the path to the output folder.

            `identifier` is "Results" + time and date.

            `root_input_folder` is the root of the input folder.
    :rtype: Path, string, Path
    """      

    identifier = " Results - "+str(datetime.now().strftime('%Y%m%d_%H%M%S')) #String added to output file path

    output_folder = Path(input_folder.parent).joinpath(input_folder.name + identifier) #Set path for root output folder
    root_input_folder = os.path.dirname(input_folder) #Get folder directly before root input folder

    if os.path.dirname(root_input_folder) == root_input_folder: #Prevent conflict in case input folder is placed at the root of a drive
        root_input_folder = input_folder

    if input_folder =='':
        input_folder = root_input_folder

    return output_folder,identifier,root_input_folder

def csv_sniffer(file_path):
    """Get delimiter from .csv file.

    :param file_path: Path to file.
    :type file_path: Path or string.
    :return: Delimiter.
    :rtype: string
    """

    with open(file_path, newline='') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1))
    return dialect.delimiter

def print_pb(text,i,max):
    """Print `text` into console on top of a progress bar.

    :param text: Text to be printed.
    :type text: str
    :param i: Progress bar index.
    :type i: int
    :param max: Maximum index of the progress bar.
    :type max: int
    """    
    
    m_max = max #(max-1) if max > 1 else 1
    print("\033[K",end="")
    print(text)
    print(f'| {("▊"*int(i/m_max*60))}{"_"*int((1-(i/m_max))*60)} | {i}/{max} | {round((i/m_max)*100,2)}%',end='\r',flush=True)

def get_file_list(input_folder, extension):
    """Build a list of file paths and names to loop over.

    :param input_folder: Path to folder or file
    :type input_folder: Path or string
    :param extension: Extension of files to loop over.
    :type extension: str
    :return: Lists of file paths and names.
    :rtype: list
    """ 

    path_list = []
    name_list = []

    if input_folder.endswith(f'.{extension}'):
        path_list.append(Path(input_folder))
        name_list.append(Path(input_folder).name)

    else:
        for path, subfolder, files in os.walk(input_folder): # Scan entire folder structure for files
                for name in files:
                    if name.endswith(f'.{extension}') == False:  # Check for correct file extension
                        continue # Skip to next file if not correct extension        

                    file_path = os.path.join(path, name) # Get file path of current file
                    path_list.append(file_path) # Append to file path list
                    name_list.append(name)
    
    return path_list, name_list

try:
    import imagej

    IJ_INSTALLED = True

    def pyimagej_reader(path):

        path = str(path) #jpype doesn't like PosixPath objects, apparently

        ij = imagej.init()
        jimage = ij.io().open(path)
        frames = ij.py.from_java(jimage)

        try:
            frames.shape
        except AttributeError:
            warnings.warn(f'Extension {path.split(".")[-1]} is not supported')
        if len(frames.shape) != 3:
            warnings.warn(f'File has {len(frames.shape)} dimension(s), needs 3')

        return frames
    
except ImportError:

    IJ_INSTALLED = False

    def pyimagej_reader(path):
        print('PyImageJ not installed, Bio-Formats reader not available') # Shouldn't print 

def get_frames(path):

    extension = Path(path).suffix

                                                # match extension:
    if extension == '.tif':     #     case 'tif':
        frames = imageio.volread(path)     #         frames = imageio.volread(file_path)
    elif extension == '.nd2':   #     case 'nd2':
        frames = nd2.imread(path)          #         frames = nd2.imread(file_path)
    elif IJ_INSTALLED == True:
        frames = pyimagej_reader(path)
    else:
        warnings.warn(f'Extension {extension} is not supported')

    return frames

