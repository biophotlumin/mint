"""
Module regrouping several small utility functions.
"""

#Imports
import csv
import os
import yaml
import shutil
import warnings
import numpy as np

from datetime import datetime
from pathlib import Path, PosixPath, WindowsPath

Path_type = str | Path | PosixPath | WindowsPath

def extraction_csv(
        input_folder: Path_type,
        ) -> None:
    """
    Extracts .csv files from a folder into another.

    Parameters
    ----------
    input_folder : str or Path
        Path to the input folder.
    """

    input_folder = str(input_folder)

    shutil.copytree(input_folder, input_folder+" CSV ONLY", symlinks=False,
                    ignore=shutil.ignore_patterns('*.png', '*.txt'))

def folder_structure_creation(
        input_folder: Path_type,
        ) -> tuple:
    """
    Generates paths to output folders.

    Parameters
    ----------
    input_folder : str or Path
        Path to the input folder.

    Returns
    -------
    output_folder : Path
        Path to the output folder.
    identifier : str
        "Results" + time and date.
    root_input_folder : Path
        Root of the input folder.
    """

    input_folder = Path(input_folder)

    # String added to output file path
    identifier = " Results - "+str(datetime.now().strftime('%Y%m%d_%H%M%S'))

    # Set path for root output folder
    output_folder = Path(input_folder.parent).joinpath(input_folder.name + identifier)
    # Get folder directly before root input folder
    root_input_folder = os.path.dirname(input_folder)

    # Prevent conflict in case input folder is placed at the root of a drive
    if os.path.dirname(root_input_folder) == root_input_folder:
        root_input_folder = input_folder

    if input_folder == '':
        input_folder = root_input_folder

    return output_folder, identifier, root_input_folder

def csv_sniffer(
        file_path: Path_type,
        ) -> str:
    """
    Detect the delimiter used in a CSV file.

    Parameters
    ----------
    file_path : Path or string
        Path to file.

    Returns
    -------
    delimiter : str
        Detected delimiter.
    """

    with open(file_path, newline='') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1))
    return dialect.delimiter

def print_pb(
        text: str,
        i: int,
        max: int,
        ) -> None:
    """
    Print text on top of a progress bar.

    Parameters
    ----------
    text : str
        Text to be printed.
    i : int
        Current index.
    max : int
        Maximum index.
    """

    m_max = max #(max-1) if max > 1 else 1
    print("\033[K", end="")
    print(text)
    print(f'| {("▊"*int(i/m_max*60))}{"_"*int((1-(i/m_max))*60)} | {i}/{max} |'
          f' {round((i/m_max)*100, 2)}%', end='\r', flush=True)

def get_file_list(
        input_folder: Path_type,
        extension: str | list,
        ) -> tuple[list, list]:
    """
    Build lists of file paths and names to loop over.

    Parameters
    ----------
    input_folder : Path or str
        Path to folder or file.
    extension : str
        Extension of files to loop over.

    Returns
    -------
    list
        List of file paths.
    list
        List of file names.
    """

    paths = (p for p in Path(input_folder).glob("**/*") if p.suffix in
             (set(extension) if isinstance(extension, list) else set([extension])))
    paths = [(str(path), path.name) for path in paths]
    parents, names = zip(*paths)

    return parents, names

def load_params(
        path: Path_type,
        ) -> dict:
    """
    Load parameters from a yaml file.

    This function loads parameters from a yaml file located at the
    given path. If the file does not exist, returns an empty dictionary.

    Parameters
    ----------
    path : Path or str
        Path to the yaml file.

    Returns
    -------
    dict
        Dictionary containing the loaded parameters.
    """

    try:
        config = yaml.safe_load(open(path))
    except FileNotFoundError:
        warnings.warn(f'File not found : {path}')
        return {}

    return config

def dict_dump(
        path,
        data: dict,
        file_name: str,
        overwrite: bool=False,
        ) -> None:
    """
    Dump a dictionary to a YAML file.

    Parameters
    ----------
    path : str or Path
        Path to the folder where the YAML file will be saved.
    data : dict
        Dictionary to dump.
    file_name : str
        Name of the YAML file.
    overwrite : bool, optional (default=False)
        Overwrite existing file if True, otherwise update the preexisting
        dictionary.

    """

    data = data.copy() # Just in case the dict has to be modified

    for k, v in data.items(): # YAML doesn't like PosixPath, convert to string
        if isinstance(v, PosixPath):
            data[k] = str(data[k])
        elif isinstance(v, np.floating):
            data[k] = float(data[k])

    file_path = Path(path).joinpath(f'{file_name}.yml')
    if file_path.is_file(): # Check for existing file
        if overwrite is True:
            with open(file_path, 'w') as f: # Overwrite if required
                yaml.dump(data, f)
        else:
            old_dict = yaml.safe_load(open(file_path))
            old_dict.update(data) # Otherwise update the preexisting dict
            with open(file_path, 'w') as f:
                yaml.dump(old_dict, f)
    else:
        with open(file_path, 'w') as f: # Otherwise create the file
            yaml.dump(data, f)

def dict_load(
        path: Path_type,
        name: str,
        ) -> dict:
    """
    Loads a YAML file from the specified path and returns the contents as a dictionary.

    Parameters
    ----------
    path : str
        The path to the directory containing the YAML file.
    name : str
        The name of the YAML file (without the extension).

    Returns
    -------
    dict
        The contents of the YAML file as a dictionary.
    """
    return yaml.safe_load(open(Path(path).joinpath(f'{name}.yml')))

class Logger():
    def __init__(self):
        self.logged = {}

    def log(self, k, v, type):
        if k not in self.logged.keys():
            if type == 'append':
                self.logged[k] = [v]
            else:
                self.logged[k] = v

        elif type == 'add':
            self.logged[k] += v
        elif type == 'append':
            self.logged[k].append(v)

    def dump(self, path):
        dict_dump(path=path, data=self.logged, file_name='log2')

    def get(self, k):
        return self.logged.get(k, [])

# global logger
logger = Logger()
