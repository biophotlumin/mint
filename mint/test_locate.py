"""Tests `locate` parameters.
"""

#Imports
import os
import matplotlib
import numpy as np
import trackpy as tp
import matplotlib.pyplot as plt

from pathlib import Path
from .input import get_frames
from .utils import get_file_list
from .denoising import filtering

# plt.switch_backend('TkAgg')

def get_range(parameter: int | list, scaling: int) -> list | np.ndarray:
    """
    Create a sequence of numbers according to the parameter.

    Parameters
    ----------
    parameter : int, float, or list of two numbers
        If `parameter` is a list of two ints, returns a range between the two numbers,
        with `scaling` as the step size.
    scaling : float
        The step size between the start and end points of the range
        when a list is provided.

    Returns
    -------
    result : list or ndarray
        Parameter or list of parameters

    Raises
    ------
    ValueError
        Raised if the specified argument do not match the expected format.
    """

    if isinstance(parameter, int) or isinstance(parameter, float):
        return [parameter]
    elif isinstance(parameter, list) and len(parameter) == 1:
        return parameter
    elif isinstance(parameter, list) and len(parameter) == 2:
        return np.arange(parameter[0], parameter[1], scaling)
    else:
        raise ValueError(f'Invalid parameter : {parameter}')


def test_locate(input_folder, parameters, settings):
    """
    Test locate parameters.

    Parameters
    ----------
    input_folder : Path or str
        Path to input folder containing raw videos.
    parameters : dict
        Dictionary containing calculation parameters.
    settings : dict
        Dictionary containing calculation settings.

    """

    path_list, name_list = get_file_list(input_folder, parameters['extension_in'])

    for (path, name, j) in zip(path_list, name_list,
                               [j for j in range(len(path_list))]):

        print(f'\nProcessing {name}')
        frames = get_frames(path)

        if len(frames) == 0:
            return

        processed_frame = frames[0]

        # Prevent conflict in case of wavelet filtering
        processed_frames = processed_frame.astype('float64')

        processed_frames = filtering(processed_frames, settings, parameters)

        # Show filter results

        if settings['show_filter']:
            plt.imshow(processed_frame)
            plt.title("Filter test on "+name, fontsize=10)
            plt.xlabel("x (pixels)", fontsize=10)
            plt.ylabel("y (pixels)", fontsize=10)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.show()
            plt.close()

        # Locating particles and finding trajectories

        tp.quiet(suppress=True) # Silencing TrackPy messages

        print('\tLocating')

        diameter_range = get_range(parameters['diameter'],
                                   parameters['diameter_scaling'])
        minmass_range = get_range(parameters['minmass'],
                                  parameters['minmass_scaling'])
        separation_range = get_range(parameters['separation'],
                                     parameters['separation_scaling'])

        output_folder = Path(input_folder).parent.joinpath(
            f'{name} - Locate parameters test')
        if not Path.is_dir(output_folder):
            os.makedirs(output_folder)

        matplotlib.use('Agg')
        values_tup = [(diameter, minmass, separation) for diameter in diameter_range
                      for minmass in minmass_range for separation in separation_range]
        for diameter, minmass, separation in values_tup:
            raw_coordinates = tp.locate(processed_frame,
                                        minmass=minmass,
                                        diameter=diameter,
                                        separation=separation,
                                        preprocess=False,
                                        engine='numba')

            plt.imshow(processed_frame)
            plt.title("Locate parameters test on "+name, fontsize=10)
            plt.xlabel("x (pixels)", fontsize=10)
            plt.ylabel("y (pixels)", fontsize=10)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            tp.annotate(raw_coordinates,
                        processed_frame,
                        plot_style={'markersize': 10},
                        color='red')
            plt.figtext(0.5, 0.01,
                        f'Number of particles found :  {len(raw_coordinates)}',
                        ha='center')
            plt.savefig(Path(output_folder).joinpath(f'{name}_diameter_{diameter}_minmass_{minmass}_separation_{separation}.png'))
            plt.close()

if __name__ == '__main__':

    parameters = {
        #trackpy.batch
        'diameter': [7],
        'minmass': [200],
        'separation': [12],
        'filter_separation': 12,
        'diameter_scaling': 2,
        'minmass_scaling': 10,
        'separation_scaling': 1,
        'extension_in': 'nd2',
    }

    #Optional image processing

    settings = {
        'tophat': True,
        'wavelet': True,
        'show_filter': False,
    }

    input_folder = r''
    test_locate(input_folder, parameters, settings=settings)