"""Module regrouping several small utility functions.
"""

#Imports
import csv
import os
import shutil
from datetime import datetime
from pathlib import Path, PosixPath, WindowsPath

Path_type = str | Path | PosixPath | WindowsPath

def extraction_csv(input_folder: Path_type) -> None:
    """Extracts .csv files from a folder into another.

    :param input_folder: Path to the input folder.
    :type input_folder: string or Path
    """

    print(input_folder)
    shutil.copytree(input_folder, input_folder+" CSV ONLY", symlinks=False,
                    ignore=shutil.ignore_patterns('*.png', '*.txt'))

def folder_structure_creation(input_folder: Path_type) -> None:
    """Generates paths to output folders.

    Adds "Results" and the time and date of the start of the run to the input folder.

    :param input_folder: Path to the input folder.
    :type input_folder: string or Path
    :return: `output_folder` is the path to the output folder.

            `identifier` is "Results" + time and date.

            `root_input_folder` is the root of the input folder.
    :rtype: Path, string, Path
    """

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

def csv_sniffer(file_path: Path_type) -> str:
    """Get delimiter from .csv file.

    :param file_path: Path to file.
    :type file_path: Path or string.
    :return: Delimiter.
    :rtype: string
    """

    with open(file_path, newline='') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1))
    return dialect.delimiter

def print_pb(text: str, i: int, max: int) -> None:
    """Print `text` into console on top of a progress bar.

    :param text: Text to be printed.
    :type text: str
    :param i: Progress bar index.
    :type i: int
    :param max: Maximum index of the progress bar.
    :type max: int
    """

    m_max = max #(max-1) if max > 1 else 1
    print("\033[K", end="")
    print(text)
    print(f'| {("▊"*int(i/m_max*60))}{"_"*int((1-(i/m_max))*60)} | {i}/{max} |'
          f' {round((i/m_max)*100, 2)}%', end='\r', flush=True)

def get_file_list(input_folder: Path_type, extension: str) -> tuple[list, list]:
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
        for path, subfolder, files in os.walk(input_folder):
                for name in files:
                    if name.endswith(f'.{extension}') is False:
                        continue

                    file_path = os.path.join(path, name)
                    path_list.append(file_path)
                    name_list.append(name)

    return path_list, name_list