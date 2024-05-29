"""
Functions used to output calculation results into files.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

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