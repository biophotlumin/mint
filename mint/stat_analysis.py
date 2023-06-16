"""Functions used for the statistical analysis of extracted transport parameters.
"""

# Imports
import os
import numpy as np
import pandas as pd
import seaborn as sns
import scikit_posthocs as sp
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from utils import *
from scipy import stats
from pathlib import Path
from statistics import mean
from output import dict_dump

vars = { # Variables of interest as keys, labels as values
        'pausing_time':'Pausing time (s)',
        'pausing_frequency':'Pausing frequency (events/min)',
        'curvilign_velocity':'Curvilign velocity (µm/s)',
        'curv_velocity_antero':'Anterograde curvilign velocity (µm/s)',
        'curv_velocity_retro':'Retrograde curvilign velocity (µm/s)',
        'processivity':'Processivity (s)',
        'processivity_antero':'Anterograde processivity (s)',
        'processivity_retro':'Retrograde processivity (s)',
        'run_length':'Run length (µm)',
        'run_length_antero':'Anterograde run length (µm)',
        'run_length_retro':'Retrograde run length (µm)',
        'diag_size':'Length of trajectory (µm)',
        'fraction_paused':'Fraction of time paused',
        'directionality':'Ratio of retrograde to anterograde transport',
        'switch':'Directionality reversal',
        'variance_go':'Variance of intensity in GO phases',
        'variance_stop':'Variance of intensity in STOP phases',
        'duration':'Trajectory duration (s)',
        'curvilign_length':'Curvilign length (µm)',
        'switch_a_to_r':'Anterograde to retrograde reversal',
        'switch_r_to_a':'Retrograde to anterograde  reversal',
        'phase_dir_go':'Number of retrograde phases over total number of GO phases',
        'phase_dir_total':'Number of retrograde phases over total number of phases',
        'switch_normal':'Directionality reversal per µm',
        'switch_var_stop':'Variation of intensity in STOP phases of reversal',
        'theta_std_go':'Standard deviation of θ in GO phases', # Theta might cause issues on Windows
        'theta_std_stop':'Standard deviation of θ in STOP phases',
        'pausing_time_antero':'Pausing time in anterograde motion (s)',
        'pausing_time_retro':'Pausing time in retrograde motion (s)',
        'pausing_time_switch':'Pausing time in bidirectional motion (s)',
        'fraction_moving':'Fraction of moving particles',
        
}

statistical_tests = {'kruskal':'Kruskal-Wallis ', #Function names as keys, label as values (with a trailing space)
                     'ranksums':'rank-sums ',
                     't_test':'t-',}

def statistical_analysis(settings,parameters,input_folder):
    """Perform statistical analysis.

    Scans through .csv files and runs them through statistical analysis.

    :param settings: Dictionary containing settings.
    :type settings: dict
    :param parameters: Dictionary containing parameters.
    :type parameters: dict
    :param input_folder: Folder containing .csv file to be analysed.
    :type input_folder: str or Path
    """    

    for path, subfolder, files in os.walk(input_folder): # Scan entire folder structure for files
        for name in files:
            if name.endswith('rajectory average parameters.csv') == False:  # Check for correct file
                continue # Skip to next file if not correct .csv 

            file_path = os.path.join(path, name)
            data = pd.read_csv(file_path,sep=csv_sniffer(file_path))

            act_variables = {} # Get a dict of variables actually contained in the .csv
            for k, v in vars.items():
                if k in data.columns:
                    act_variables[k] = v

            run_stats(settings,parameters,act_variables,data,input_folder)

def run_stats(settings,parameters,act_variables,data,input_folder):
    """Statistical tests and plotting functions for each variable of interest.

    :param ordering: Optional order of experimental conditions in graphs.
    :type ordering: List of str
    :param parameters: Dictionary containing parameters.
    :type parameters: dict
    :param act_variables: Dictionary of variables to be analysed.
    :type act_variables: dict
    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param input_folder: Folder containing .csv file to be analysed.
    :type input_folder: str or Path
    :raises RuntimeError: Raises error in case the order list does not match actual data.
    """    

    # Uncomment for subpopulation analysis

    # non_antero = data.loc[data['directionality']>0] #Everything but purely anterograde trajectories
    # non_retro = data.loc[data['directionality']<1] #Everything but purely retrograde trajectories
    # mixed = non_antero.loc[non_antero['directionality']<1] #Everything but unidirectional trajectories
    # pure_retro = data.loc[data['directionality']==1] #Purely retrograde trajectories
    # pure_antero = data.loc[data['directionality']==0] #Purely anterograde trajectories
    # left = data.loc[data['slide']=='oeil_gauche'] #Left eye only
    # right = data.loc[data['slide']=='oeil_droit'] #Right eye only
    # data = ?

    if settings['ordering']:
        order = parameters['order'] # Optionally order by condition
        if len(order) != data.condition.nunique():
            raise RuntimeError('Length of order list does not match number of conditions')
        if any(c not in data.condition.unique() for c in order):
            raise RuntimeError('Item(s) in order list not found in conditions')
        if len(order) != len(set(order)):
            raise RuntimeError('Order list contains duplicates')

        data['condition_sorted'] = data.condition.astype("category")
        data.condition_sorted = data.condition_sorted.cat.set_categories(order)
        data.sort_values(['condition_sorted'],inplace=True)
    
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
    
    dict_dump(input_folder,log_stats,'log')
    dict_dump(input_folder,vars,'vars')

    if settings['antero_retro']:
        # Invert retrograde variables to display them with positive values
        data.run_length_retro = data.run_length_retro*-1
        data.curv_velocity_retro = data.curv_velocity_retro*-1

    results = []

    print(f'Calculating stats for :')
    for k,v in act_variables.items():
        print(f'\t{k}')

        if k == 'fraction_moving':
            data = data.drop_duplicates(subset='file') # Must be kept last in the vars dict as to not interfere with the others

        if data.condition.nunique() > 2:
            if normality(data,k) == False: # Check for normality
                results.append(run_variable(k,False,'kruskal',data,True,input_folder,parameters))
            elif normality(data,k) == True:
                results.append(run_variable(k,True,'kruskal',data,True,input_folder,parameters))

        elif data.condition.nunique() == 2:
            if normality(data,k) == False: # Check for normality
                results.append(run_variable(k,False,'ranksums',data,False,input_folder,parameters))
            elif normality(data,k) == True:
                results.append(run_variable(k,True,'t_test',data,False,input_folder,parameters))

        else:
            results.append(run_variable(k,False,None,data,False,input_folder,parameters))

    results = [item for sublist in results for item in sublist] # Flatten list of list

    text_file = open((Path(input_folder).joinpath("Statistical test results.txt")),'w')
    text_file.writelines(results)
    text_file.close()
    
def run_variable(var,normal,test,data,dunn_b,input_folder,parameters):
    """Call appropriate statistical test and plotting function for a variable.

    :param var: Current variable.
    :type var: str
    :param normal: Wether the distribution of `var` in `data` is normal.
    :type normal: bool
    :param test: Statistical test to be applied.
    :type test: str
    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param dunn_b: Wether Dunn's test should be applied.
    :type dunn_b: bool
    :param input_folder: Folder where the files will be saved.
    :type input_folder: str or Path
    :param parameters: Dictionary containing parameters
    :type parameters: dict
    :return: List of results.
    :rtype: list
    """    
    var_res = []
    var_res.append(f'{var.upper()}\n')

    if normal:
        var_res.append(f'Distribution of {str(var)} is normal \n')
    else:
        var_res.append(f'Distribution of {str(var)} is not normal \n')

    if test:
        p_value = eval(test)(data,var)
        var_res.append(f'p-value of {statistical_tests[test]}test for {var} is {str(round(p_value,6))}\n')
    else:
        p_value = 1

    if dunn_b:
        var_res.append('\n'+dunn(data,var)+'\n\n')

    boxplot(data,var,input_folder,str(round(p_value,6)),parameters)
    barplot(data,var,input_folder,str(round(p_value,6)),parameters)

    var_res.append(f'Sample size of {var} is {round(means(data,var)[0],6)}\n')
    var_res.append(f'Mean of {var} is {round(means(data,var)[1],6)}\n')
    var_res.append(f'Median of {var} is {round(means(data,var)[2],6)}\n')
    var_res.append(f'Standard deviation of {var} is {round(means(data,var)[3],6)}\n')
    var_res.append('\n')

    var_df = pd.DataFrame(columns=['Condition','Sample size','Mean','Median','Standard deviation'],dtype='object')
    var_df.set_index('Condition',inplace=True)

    for cond in data.condition.unique():
        subdata = data[data.condition == cond]
        var_df.loc[cond] = [round(i,6) for i in [*means(subdata,var)]]

    var_res.append(var_df.to_string()+'\n\n')

    return var_res

def ranksums(data,variable):
    """Two-sided Mann-Whitney U test

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :return: p-value.
    :rtype: float
    """    

    index = data.condition.unique()
    x = data.loc[data.condition==index[0], variable]
    y = data.loc[data.condition==index[1], variable]
    p = stats.mannwhitneyu(x,y,use_continuity=True,alternative='two-sided',nan_policy='omit')[1]

    return p

def kruskal(data,variable):
    """Kruskal-Wallis test.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :return: p-value.
    :rtype: float
    """    

    list_of_arrays = []
    for index in data.condition.unique():
        list_of_arrays.append(data.loc[data.condition==index, variable])

    p = stats.kruskal(*list_of_arrays,nan_policy='omit')[1]
    return p

def t_test(data,variable):
    """t-test.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :return: p-value.
    :rtype: float
    """  

    list_of_arrays = []
    for index in data.condition.unique():
        list_of_arrays.append(data.loc[data.condition==index, variable])

    p = stats.ttest_ind(*list_of_arrays,nan_policy='omit')[1]
    return p

def dunn(data,variable):
    """Dunn's test.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :return: Table of p-values.
    :rtype: str
    """  

    list_of_arrays = []
    conditions_list = []
    for index in data.condition.unique():
        list_of_arrays.append(data.loc[data.condition==index, variable])
        conditions_list.append(index)

    p = sp.posthoc_dunn(list_of_arrays)

    index = []
    for i in range(len(list_of_arrays)): # Rename DataFrame with analysed conditions
        index.append(i+1)
    label_dict = dict(zip(index, conditions_list))
    p.rename(index=label_dict,columns=label_dict,inplace=True)

    return p.to_string()

def boxplot(data,variable,input_folder,p,parameters):
    """Box plot.

    Generate a box plot, with a notch at the median, and saves it as a file under `extension_out`.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :param input_folder: Folder where the file will be saved.
    :type input_folder: str or Path
    :param p: p-value.
    :type p: float
    :param parameters: Dictionary containing parameters.
    :type parameters: dict
    """

    sns.set_theme(style="ticks", palette=sns.color_palette("crest"))
    sns.boxplot(x=data['condition'],
            y=data[variable],  width=0.35, notch=True,
            data=data, showfliers =False,linewidth=1)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(vars[variable])
    plt.annotate((f'p-value : {p}'),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath(f'Boxplot {str(data["condition"].unique())} {variable}.{parameters["extension_out"]}'),dpi=parameters["dpi"])
    plt.close()
    
def barplot(data,variable,input_folder,p,parameters):
    """Bar plot.

    Generate a bar plot, with SEM as error bars, and saves it as a file under `extension_out`.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :param input_folder: Folder where the file will be saved.
    :type input_folder: str or Path
    :param p: p-value.
    :type p: float
    :param parameters: Dictionary containing parameters.
    :type parameters: dict
    """

    error = []
    
    for i in data['condition'].unique():
        error.append(stats.sem(data[variable].loc[data['condition']==i],nan_policy='omit') if len(data)>3 else 0)
    sns.barplot(y=data[variable].dropna(),x=data['condition'],estimator=mean,yerr=error, errorbar=None,\
            error_kw={'elinewidth':2,'capsize':4,'capthick':2})
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(vars[variable])
    plt.annotate((f'p-value : {p}'),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath(f'Barplot {str(data["condition"].unique())} {variable}.{parameters["extension_out"]}'),dpi=parameters["dpi"])
    plt.close()

def violinplot(data,variable,p):
    """Violin plot.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :param p: p-value.
    :type p: float
    """
    
    sns.violinplot(x=data['condition'],y=data[variable],inner='box',cut=0)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(vars[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    plt.show()

def normality(data,variable):
    """Check for normality.

    Test wether or not a sample fits a normal distribution.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :return: Wether the sample fits a normal distribution.
    :rtype: bool
    """    

    if len(data) >= 8:
        p = stats.normaltest(a=data[variable],nan_policy='omit')[1]
        if p > 0.05:
            return True
        else:
            return False
    else:
        return False

def means(data,variable):
    """Calculate distribution parameters.

    Calculate the sample size, mean, median and standard deviation of `variable` in `data`.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :return: Sample size, mean, median, standard deviation.
    :rtype: float
    """    

    subdata = data[variable]
    sample_size = len(subdata.dropna())
    var_mean = np.nanmean(subdata)
    var_median = np.nanmedian(subdata)
    var_std = np.nanstd(subdata)
    
    return sample_size, var_mean, var_median, var_std     

if __name__ == '__main__':

    input_folder = r''
    parameters = {'antero_retro':True,
                'extension_out':'svg',
                'dpi':300,
                'order':['WT','HET','HOM'],
                }
    settings = {'ordering':True,}
    statistical_analysis(settings,parameters,input_folder)
