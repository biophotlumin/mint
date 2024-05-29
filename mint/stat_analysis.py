"""
Functions used for the statistical analysis of extracted transport parameters.
"""

# Imports
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import seaborn as sns
from scipy import stats

from .utils import Path_type, csv_sniffer, dict_dump

plt.switch_backend('agg')

stats_vars = { # Variables of interest as keys, labels as values
        'pausing_time': 'Pausing time (s)',
        'pausing_frequency': 'Pausing frequency (events/min)',
        'curvilign_velocity': 'Curvilign velocity (µm/s)',
        'curv_velocity_antero': 'Anterograde curvilign velocity (µm/s)',
        'curv_velocity_retro': 'Retrograde curvilign velocity (µm/s)',
        'processivity': 'Processivity (s)',
        'processivity_antero': 'Anterograde processivity (s)',
        'processivity_retro': 'Retrograde processivity (s)',
        'run_length': 'Run length (µm)',
        'run_length_antero': 'Anterograde run length (µm)',
        'run_length_retro': 'Retrograde run length (µm)',
        'diag_length': 'Length of trajectory (µm)',
        'fraction_paused': 'Fraction of time paused',
        'directionality': 'Ratio of retrograde to anterograde transport',
        'switch': 'Directionality reversal',
        'variance_go': 'Variance of intensity in GO phases',
        'variance_stop': 'Variance of intensity in STOP phases',
        'duration': 'Trajectory duration (s)',
        'curvilign_length': 'Curvilign length (µm)',
        'switch_a_to_r': 'Anterograde to retrograde reversal',
        'switch_r_to_a': 'Retrograde to anterograde  reversal',
        'phase_dir_go': 'Number of retrograde phases over total number of GO phases',
        'phase_dir_total': 'Number of retrograde phases over total number of phases',
        'switch_normal': 'Directionality reversal per µm',
        'switch_var_stop': 'Variation of intensity in STOP phases of reversal',
        # Theta might cause issues on Windows despite being UTF-8
        'theta_std_go': 'Standard deviation of θ in GO phases',
        'theta_std_stop': 'Standard deviation of θ in STOP phases',
        'pausing_time_antero': 'Pausing time in anterograde motion (s)',
        'pausing_time_retro': 'Pausing time in retrograde motion (s)',
        'pausing_time_switch': 'Pausing time in bidirectional motion (s)',
        # Must be kept last in the dict !
        'fraction_moving': 'Fraction of moving particles',
        # 'fraction_moving_msd': 'Fraction of moving particles (MSD)'

}
# Function names as keys, label as values (with a trailing space)
statistical_tests = {'kruskal': 'Kruskal-Wallis ',
                     'ranksums': 'rank-sums ',
                     't_test': 't-',
                     }

def statistical_analysis(
        settings: dict,
        parameters: dict,
        input_folder: Path_type
        ) -> None:
    """
    Scans through .csv files and runs them through statistical analysis.

    Parameters
    ----------
    settings : dict
        Dictionary containing settings.
    parameters : dict
        Dictionary containing parameters.
    input_folder : str or Path
        Folder containing .csv file to be analysed.
    """

    for path, subfolder, files in os.walk(input_folder):
        for name in files:
            if name.endswith('rajectory_parameters.csv') is False:
                continue

            file_path = os.path.join(path, name)
            data = pd.read_csv(file_path, sep=csv_sniffer(file_path))

            act_variables = {} # Get a dict of variables actually contained in the .csv
            for k, v in stats_vars.items():
                if k in data.columns:
                    act_variables[k] = v

            run_stats(settings, parameters, act_variables, data, input_folder)

def run_stats(
        settings: dict,
        parameters: dict,
        act_variables: dict,
        data: pd.DataFrame,
        input_folder: Path_type,
        ) -> None:
    """
    Statistical tests and plotting functions for each variable of interest.

    Parameters
    ----------
    settings : dict
        Dictionary containing settings.
    parameters : dict
        Dictionary containing parameters.
    act_variables : dict
        Dictionary of variables to be analysed.
    data : pandas.DataFrame
        DataFrame containing trajectory parameters.
    input_folder : str or Path
        Folder containing .csv file to be analysed.

    Raises
    ------
    RuntimeError
        Raised if the conditions in the order list does not match the data.
    """


    # Uncomment for subpopulation analysis

    # Everything but purely anterograde trajectories
    # non_antero = data.loc[data['directionality']>0]
    # Everything but purely retrograde trajectories
    # non_retro = data.loc[data['directionality']<1]
    # Everything but unidirectional trajectories
    # mixed = non_antero.loc[non_antero['directionality']<1]
    # Purely retrograde trajectories
    # pure_retro = data.loc[data['directionality']==1]
    # Purely anterograde trajectories
    # pure_antero = data.loc[data['directionality']==0]
    # Left eye only
    # left = data.loc[data['slide']=='oeil_gauche']
    # Right eye only
    # right = data.loc[data['slide']=='oeil_droit']
    if settings['gfp']:
        gfp_pos = data.loc[data['gfp'] == "[ True]"]
        data = gfp_pos

    if settings['ordering']:
        order = parameters['order'] # Optionally order by condition

        if len(order) != data.condition.nunique():
            raise RuntimeError(
                'Length of order list does not match number of conditions')
        if any(c not in list(data.condition.unique()) for c in order):
            raise RuntimeError(
                'Item(s) in order list not found in conditions')
        if len(order) != len(set(order)):
            raise RuntimeError(
                'Order list contains duplicates')

        data['condition_sorted'] = data.condition.astype("category")
        data.condition_sorted = data.condition_sorted.cat.set_categories(order)
        data.sort_values(['condition_sorted'], inplace=True)

    log_stats = {}

    log_stats['n_a_animal'] = data.animal.nunique()
    log_stats['n_a_file'] = data.file.nunique()
    log_stats['n_a_traj'] = data.trajectory.nunique()
    log_stats['raw_order'] = [x for x in data.condition.unique()]

    for condition in data.condition.unique():
        subdata = data[data.condition == condition]
        log_stats[f'n_a_animal_{condition}'] = subdata.animal.nunique()
        log_stats[f'n_a_file_{condition}'] = subdata.file.nunique()
        log_stats[f'n_a_traj_{condition}'] = subdata.trajectory.nunique()

    dict_dump(input_folder, log_stats, 'log')
    dict_dump(input_folder, stats_vars, 'vars')

    if settings['antero_retro']:
        # Invert retrograde variables to display them with positive values
        data.run_length_retro = data.run_length_retro*-1
        data.curv_velocity_retro = data.curv_velocity_retro*-1

    results = []

    print('Calculating stats for : ')
    for k, v in act_variables.items():
        print(f'\t{k}')

        # Must be kept last in the stats_vars dict as to not interfere with the others
        if k == 'fraction_moving':
            data = data.drop_duplicates(subset='file')

        if data.condition.nunique() > 2:
            if normality(data, k) is False:  # Check for normality
                results.append(run_variable(k, False, 'kruskal', data, True,
                                            input_folder, parameters))
            elif normality(data, k) is True:
                results.append(run_variable(k, True, 'kruskal', data, True,
                                            input_folder, parameters))

        elif data.condition.nunique() == 2:
            if normality(data, k) is False:  # Check for normality
                results.append(run_variable(k, False, 'ranksums', data, False,
                                            input_folder, parameters))
            elif normality(data, k) is True:
                results.append(run_variable(k, True, 't_test', data, False,
                                            input_folder, parameters))

        else:
            results.append(run_variable(k, False, None, data, False,
                                        input_folder, parameters))

    results = [item for sublist in results for item in sublist] # Flatten list of list

    text_file = open((Path(input_folder).joinpath("Statistical test results.txt")), 'w')
    text_file.writelines(results)
    text_file.close()

def run_variable(
        var: str,
        normal: bool,
        test: str | None,
        data: pd.DataFrame,
        dunn_b: bool,
        input_folder: Path_type,
        parameters: dict,
        ) -> list:
    """
    Call appropriate statistical test and plotting function for a variable.

    Parameters
    ----------
    var : str
        Current variable.
    normal : bool
        Whether the distribution of `var` in `data` is normal.
    test : str
        Statistical test to be applied.
    data : DataFrame
        DataFrame containing trajectory parameters.
    dunn_b : bool
        Whether Dunn's test should be applied.
    input_folder : str or Path
        Folder where the files will be saved.
    parameters : dict
        Dictionary containing parameters

    Returns
    -------
    list
        List of results.
    """
    var_res = []
    var_res.append(f'{var.upper()}\n')

    if normal:
        var_res.append(f'Distribution of {str(var)} is normal \n')
    else:
        var_res.append(f'Distribution of {str(var)} is not normal \n')

    if test:
        p_value = eval(test)(data, var)
        var_res.append(f'p-value of {statistical_tests[test]}test for '
                       # Formatted to disable scientific notation
                       f'{var} is {p_value:.10f}\n')
    else:
        p_value = 1

    if dunn_b:
        var_res.append('\n'+dunn(data, var)+'\n\n')

    boxplot(data, var, input_folder, str(round(p_value, 6)), parameters)
    barplot(data, var, input_folder, str(round(p_value, 6)), parameters)

    var_res.append(f'Sample size of {var} is {round(means(data, var)[0], 6)}\n')
    var_res.append(f'Mean of {var} is {round(means(data, var)[1], 6)}\n')
    var_res.append(f'Median of {var} is {round(means(data, var)[2], 6)}\n')
    var_res.append(f'Standard deviation of {var} is {round(means(data, var)[3], 6)}\n')
    var_res.append('\n')

    var_df = pd.DataFrame(columns=['Condition', 'Sample size', 'Mean', 'Median',
                                   'Standard deviation'], dtype='object')
    var_df.set_index('Condition', inplace=True)

    for cond in data.condition.unique():
        subdata = data[data.condition == cond]
        var_df.loc[cond] = [round(i, 6) for i in [*means(subdata, var)]]

    var_res.append(var_df.to_string()+'\n\n')

    return var_res

def ranksums(
        data:  pd.DataFrame,
        variable:  str
        ) -> float:
    """
    Two-sided Mann-Whitney U test.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.

    Returns
    -------
    float
        p-value.
    """

    index = data.condition.unique()
    x = data.loc[data.condition == index[0], variable]
    y = data.loc[data.condition == index[1], variable]
    p = stats.mannwhitneyu(x, y, use_continuity=True, alternative='two-sided',
                           nan_policy='omit')[1] # type: ignore
                            # `nan_policy` isn't explicitely defined as an argument
    return p

def kruskal(
        data: pd.DataFrame,
        variable:  str
        ) -> float:
    """
    Calculates the Kruskal-Wallis H-test for one nominal variable and one or more
    numerical variables.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.

    Returns
    -------
    float
        p-value.
    """

    list_of_arrays = []
    for index in data.condition.unique():
        list_of_arrays.append(data.loc[data.condition == index, variable])

    p = stats.kruskal(*list_of_arrays, nan_policy='omit')[1]

    return p

def t_test(
        data:  pd.DataFrame,
        variable:  str
        ) -> float:
    """
    Calculates the Kruskal-Wallis H-test for one nominal variable and one or more
    numerical variables.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.

    Returns
    -------
    float
        p-value.
    """

    list_of_arrays = []
    for index in data.condition.unique():
        list_of_arrays.append(data.loc[data.condition == index, variable])

    p = stats.ttest_ind(*list_of_arrays, nan_policy='omit')[1]

    return p

def dunn(
        data:  pd.DataFrame,
        variable:  str
        ) -> str:
    """
    Dunn's test.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.

    Returns
    -------
    str
        Table of p-values.
    """

    list_of_arrays = []
    conditions_list = []
    for index in data.condition.unique():
        list_of_arrays.append(data.loc[data.condition == index, variable])
        conditions_list.append(index)

    p = sp.posthoc_dunn(list_of_arrays, p_adjust='fdr_bh')

    index = []
    for i in range(len(list_of_arrays)):  # Rename DataFrame with analysed conditions
        index.append(i+1)
    label_dict = dict(zip(index, conditions_list))
    p.rename(index=label_dict, columns=label_dict, inplace=True)

    return p.to_string()

def boxplot(
        data: pd.DataFrame,
        variable: str,
        input_folder: Path_type,
        p: str,
        parameters: dict,
        ) -> None:
    """
    Bar plot.

    Generate a boxplot with a notch at the median,
    and save it as a file under `extension_out`.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.
    input_folder : str or Path
        Folder where the file will be saved.
    p : str
        p-value.
    parameters : dict
        Dictionary containing parameters.
    """

    sns.set_theme(style="ticks", palette=sns.color_palette("deep")) # 'crest'
    sns.boxplot(x=data['condition'],
            y=data[variable], hue=data['condition'], width=0.35, notch=True,
            data=data, showfliers=False, linewidth=1)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(stats_vars[variable])
    plt.annotate((f'p-value : {p}'), xy=(195, 310), xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath((f'Boxplot '
                                             f'{str(data["condition"].unique())} '
                                             f'{variable}.{parameters["extension_out"]}')),
                                             dpi=parameters["dpi"])
    plt.close()

def barplot(
        data: pd.DataFrame,
        variable: str,
        input_folder: Path_type,
        p: str,
        parameters: dict,
        ) -> None:
    """
    Bar plot.

    Generate a bar plot with the mean of `variable` as bar height,
    and SEM as error bars, and save it as a file under `extension_out`.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.
    input_folder : str or Path
        Folder where the file will be saved.
    p : str
        p-value.
    parameters : dict
        Dictionary containing parameters.
    """
    error = []

    for i in data['condition'].unique():
        c_data = data[variable].loc[data['condition'] == i]
        error.append(stats.sem(c_data, nan_policy='omit')
                               if len(c_data) > 3 else 0)

    # Order isn't properly inferred from DataFrame for columns with missing values
    # since seaborn 13.0, need to pass order explicitly

    data = data.dropna(subset=variable)
    bars = sns.barplot(data=data, y=variable, x='condition', estimator=np.mean,
                       errorbar=None, hue='condition', order=parameters['order'],
                       err_kws={'elinewidth': 2, 'capsize': 4, 'capthick': 2})

    # Can't pass a list to yerr directly since seaborn 13.0,
    # need to plot errorbars independently
    # If `hue` isn't used above :
    bars = bars.containers # bars = bars.containers[0]
    # Remove [0] here
    x_coordinates = [bar[0].get_x() + bar[0].get_width() / 2 for bar in bars]
    # and here
    y_coordinates = [bar[0].get_height() for bar in bars]
    plt.errorbar(x_coordinates, y_coordinates, yerr=error, elinewidth=2, capsize=4,
                 capthick=2, fmt='None', ecolor='black')

    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(stats_vars[variable])
    plt.annotate((f'p-value :  {p}'), xy=(195, 310), xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath(f'Barplot '
                                            f'{str(data["condition"].unique())} '
                                            f'{variable}.{parameters["extension_out"]}'),
                                            dpi=parameters["dpi"])
    plt.close()

def violinplot(
        data: pd.DataFrame,
        variable: str,
        p: str,
        ) -> None:
    """
    Violin plot.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.
    p : str
        p-value.

    """

    sns.violinplot(x=data['condition'], y=data[variable], inner='box', cut=0)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(stats_vars[variable])
    plt.annotate(("p-value :  "+p), xy=(195, 310), xycoords='figure points')
    plt.show()

def normality(
        data: pd.DataFrame,
        variable: str,
        ) -> bool:
    """
    Check for normality.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.

    Returns
    -------
    bool
        True if the sample fits a normal distribution.
    """

    if len(data[variable]) >= 8:
        p = stats.normaltest(a=data[variable], nan_policy='omit')[1]
        if p > 0.05:
            return True
        else:
            return False
    else:
        return False

def means(
        data: pd.DataFrame,
        variable: str
        ) -> tuple:
    """
    Calculate the sample size, mean, median,
    and standard deviation of `variable` in `data`.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing trajectory parameters.
    variable : str
        Current variable.

    Returns
    -------
    Tuple[float, float, float, float]
        Sample size, mean, median, standard deviation.
    """

    subdata = data[variable]
    sample_size = len(subdata.dropna())
    var_mean = np.nanmean(subdata)
    var_median = np.nanmedian(subdata)
    var_std = np.nanstd(subdata)

    return sample_size, var_mean, var_median, var_std

if __name__ == '__main__':

    input_folder = r''
    parameters = {
                'extension_out': 'png',
                'dpi': 300,
                'order': ['Condition 1', 'Condition 2', 'Condition 3'],
                }
    settings = {'ordering': True,
                'antero_retro': True,
                }
    statistical_analysis(settings, parameters, input_folder)
