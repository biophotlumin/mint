"""Module regrouping several small utility functions.
"""

#Imports
import csv
import os
import yaml
import shutil
import warnings

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
        extension: str,
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
# TODO Add support for multiple extensions ?
# TODO Replace with glob
    path_list = []
    name_list = []

    if str(input_folder).endswith(f'.{extension}'): # TODO Change to Path suffix ?
        path_list.append(Path(input_folder))
        name_list.append(Path(input_folder).name)

    else:
        for path, subfolder, files in os.walk(input_folder):
                for name in files:
                    if name.endswith(f'.{extension}') is False:
                        continue

                    file_path = os.path.join(path, name)
                    path_list.append(file_path)
                    name_list.append(name)

    return path_list, name_list

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