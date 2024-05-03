"""Functions used to output calculation results into files.
"""
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path, PosixPath

from .utils import Path_type

plt.switch_backend('agg')

def trajectory_output(
        output_file_path: Path_type,
        name: str,
        process: str,
        output: pd.DataFrame,
        ) -> None:
    """
    Write trajectories into a .csv file.

    Parameters
    ----------
    output_file_path : str or Path
        Output file path.
    name : str
        Name of the file currently being processed.
    process : str
        Suffix appended to the file name.
    output : DataFrame
        DataFrame that is being dumped.
    """

    # Create file path and name
    filecsv = Path(str(output_file_path)).joinpath(name+process+".csv")
    output.to_csv(filecsv, sep='\t')

def image_output(
        output_file_path: Path_type,
        name: str,
        frames: np.ndarray,
        trajectory: pd.DataFrame,
        item: int | bool,
        ) -> None:
    """
    Plot trajectories onto the first frame of a file.

    Parameters
    ----------
    output_file_path : str or Path
        Output file path.
    name : str
        Name of the file currently being processed.
    frames : array_like
        Frames of the file currently being processed.
    trajectory : DataFrame
        DataFrame containing the trajectory to be plotted.
    item : int or bool
        Number of the trajectory to be plotted. Pass False to plot all
        trajectories.
    """

    # Initialize plot
    plt.imshow(frames[0])
    # Plotting the y axis inverts it by default, so it must be inverted again
    plt.gca().invert_yaxis()
    bbox_props = dict(boxstyle="round", fc="w", ec="w")
    arrow_props = dict(arrowstyle="simple")

    trajectories = []
    if isinstance(item, bool):
        plt.title('Trajectories from '+name)
        filepng = Path(str(output_file_path)).joinpath(name+".png")
        for item in set(trajectory.rejoined_particle):
            trajectories.append(item)
    else:
        plt.title(f'Trajectory {item} from {name}')
        trajectories.append(item)
        filepng = Path(str(output_file_path)).joinpath(name+"-"+str(item)+".png")

    for item in trajectories: # Loop for each sub trajectory
        x = trajectory[trajectory.rejoined_particle == item].x
        y = trajectory[trajectory.rejoined_particle == item].y
        plt.plot(x, y, lw=0.5)
        if len(trajectories) > 1:
            x0, y0 = x.iloc[0], y.iloc[0]
            plt.annotate(str(item), xy=(x0, y0), xytext=(x0+30, y0+30), bbox=bbox_props,
                         arrowprops=arrow_props, size=4)

    plt.savefig(filepng, dpi=300)
    plt.close()

def trajectory_separation(
        output_file_path: Path_type,
        name: str,
        trajectory: int,
        settings: dict,
        sub: pd.DataFrame,
        ) -> None:
    """
    Separates individual trajectories and writes them into .txt files.

    Parameters
    ----------
    output_file_path : str or Path
        Current output file path.
    name : str
        Name of the file currently being processed.
    trajectory : int
        Number of the trajectory being saved.
    settings : dict
        Dictionary containing settings.
    sub : DataFrame
        Subset of the DataFrame containing the trajectory.
    """

    file = Path(str(output_file_path)).joinpath(name+"-"+str(trajectory)+".txt")
    pd.options.mode.chained_assignment = None
    if settings['SNR_estimation'] is True: # Check wether SNR_estimation was used
        sub.drop(sub.columns.difference(['y', 'x', 'mass', 'feet']),
                 axis='columns', inplace=True) # If it was, drops the following columns
        sub.columns = ['Yraw', 'Xraw', 'mass', 'feet']
        sub = sub.reindex(columns=['Xraw', 'Yraw', 'mass', 'feet'])
        sub = sub.rename(columns={"mass": "Signal", "feet": "Noise"})
    else:
        # If not, drops only the columns that aren't added by SNR_estimation
        sub.drop(sub.columns.difference(['y', 'x', 'mass']),
                 axis='columns', inplace=True)
        sub.columns = ['Yraw', 'Xraw', 'mass']
        sub = sub.reindex(columns=['Xraw', 'Yraw', 'mass'])
        sub = sub.rename(columns={"mass": "Signal"})

    sub.to_csv(file, sep='\t', index=False)

def dict_dump_legacy(
        path: Path_type,
        dict_: dict,
        file_name: str,
        ) -> None:
    """
    Dump a dictionary into a text file.

    Parameters
    ----------
    path : str or Path
        Folder where the dictionary will be saved.
    dict : dict
        Dictionary to dump.
    file_name : str
        Name of the text file.
    """

    with open(Path(path).joinpath(str(file_name)+".txt"), 'a') as dict_txt:
        for k, v in dict_.items():
            dict_txt.write(f'{str(k)} : {str(v)}\n')

def dict_load_legacy(
        input_folder: Path_type,
        dict_: str,
        ) -> dict:
    """
    Load a dictionary from a text file.

    Parameters
    ----------
    input_folder : str or Path
        Folder where the dictionary is located.
    dict_ : str
        Name of text file to be loaded.

    Returns
    -------
    loaded_dict : dict
        Loaded dictionary.
    """

    loaded_dict = {}

    with open(Path(input_folder).joinpath(f'{str(dict_)}.txt')) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            k, v = line.split(' : ')
            loaded_dict[k] = v

    return loaded_dict

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

