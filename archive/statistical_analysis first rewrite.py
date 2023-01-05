"""Functions used for the statistical analysis of extracted transport parameters.

    statistical_analysis scans through .csv files and runs them through statistical analysis.
    variables calls statistical tests and plotting functions for each variable of interest.
    variables_antero_retro behaves like variables but can distinguish between anterograde and retrograde transport.
    ranksums performs a two-sided Mann-Whitney U test.
    kruskal performs a Kruskal-Wallis test.
    dunn performs Dunn's test.
    boxplot generates a boxplot, with a notch at the median, and saves it as a .png.
    barplot generates a barplot, with SEM as error bars, and saves it as a .png.
    violinplot generates a violinplot.
    normality tests wether or not a sample fits a normal distribution.
"""

## TO DO : FIGURES

#Imports
import os
import numpy as np
import pandas as pd
import seaborn as sns
import scikit_posthocs as sp
import matplotlib.pyplot as plt

from utils import *
from scipy import stats
from pathlib import Path
from statistics import mean

vars = { #Variables of interest as keys, labels as values
        'pausing_time':'Pausing time (s)',
        'pausing_frequency':'Pausing frequency (events/min)',
        'curvilign_velocity':'Curvilign velocity (µm/s)',
        'curvilign_velocity_antero':'Anterograde curvilign velocity (µm/s)',
        'curvilign_velocity_retro':'Retrograde curvilign velocity (µm/s)',
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
        'variance_GO':'Variance of intensity in GO phases',
        'variance_STOP':'Variance of intensity in STOP phases',
        'duration':'Trajectory duration (s)',
        'curvilign_length':'Curvilign length (µm)',
        'switch_a_to_r':'Anterograde to retrograde reversal',
        'switch_r_to_a':'Retrograde to anterograde  reversal',
        'phase_dir_GO':'Directionality as number of retrograde phases over total number of GO phases',
        'phase_dir_total':'Directionality as number of retrograde phases over total number of phases',
        'switch_normal':'Directionality reversal per µm',
        'switch_var_STOP':'Variation of intensity in STOP phases of reversal',
        'theta_std_GO':'theta_std_GO',
        'theta_std_STOP':'theta_std_STOP',
        'pausing_time_antero':'pausing_time_antero',
        'pausing_time_retro':'pausing_time_retro',
        'pausing_time_switch':'pausing_time_switch',
        'directionality_normal':'direction normal',
        'fraction_moving':'fraction_moving',
        
}

def statistical_analysis(settings,input_folder):
    """Scans through .csv files and runs them through statistical analysis. Inputs a dictionary and a string.

        parameters is a dictionary of calculation parameters, as defined in script.py.
        input_folder is the path to the input folder.
    """
    for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
        for name in files:
            if name.endswith('rajectory average parameters.csv') == False:  #Check for correct file extension
                continue #Skips to next file if not correct .csv 

            file_path = os.path.join(path, name)
            # print(os.path.dirname(name))
            data = pd.read_csv(file_path,sep=csv_sniffer(file_path))

            act_variables = {} #Get a dict of variables actually contained in the .csv
            for k, v in vars.items():
                if k in data.columns:
                    act_variables[k] = v

            run_stats(settings,act_variables,data,input_folder)

def run_stats(settings,act_variables,data,input_folder):
    """Calls statistical tests and plotting functions for each variable of interest. Inputs a DataFrame and a string.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        input_folder is the path to the input folder.
    """

    results = []

    #Uncomment for subpopulation analysis

    # non_antero = data.loc[data['directionality']>0] #Everything but purely anterograde trajectories
    # non_retro = data.loc[data['directionality']<1] #Everything but purely retrograde trajectories
    # mixed = non_antero.loc[non_antero['directionality']<1] #Everything but unidirectional trajectories
    # pure_retro = data.loc[data['directionality']==1] #Purely retrograde trajectories
    # pure_antero = data.loc[data['directionality']==0] #Purely anterograde trajectories
    # left = data.loc[data['slide']=='oeil_gauche'] #Left eye only
    # right = data.loc[data['slide']=='oeil_droit'] #Right eye only
    # data = ?

    order = settings['order'] #Optionally order by condition
    if len(order):
        data['condition_sorted'] = data.condition.astype("category")
        data.condition_sorted = data.condition_sorted.cat.set_categories(order)
        data.sort_values(['condition_sorted'],inplace=True)

    #Invert retrograde variables to display them with positive values
    data.run_length_retro = data.run_length_retro*-1
    data.curvilign_velocity_retro = data.curvilign_velocity_retro*-1

    for k,v in act_variables.items():
        print(k)
        if k == 'fraction_moving':
            data = data.drop_duplicates(subset='file') #Must be kept last in the vars dict as to not interfere with the others

        results.append(f'{k.upper()}\n')
        if data.condition.nunique() > 2:
            if normality(data,k) == False: #Check for normality
                results.append(f'Distribution of {str(k)} is not normal \n')
            elif normality(data,k) == True:
                results.append(f'Distribution of {str(k)} is normal \n')

            p_value = kruskal(data,k)
            results.append(f'p-value of Kruskal-Wallis test for {k} is {str(round(p_value,6))}\n')
            results.append('\n'+dunn(data,k)+'\n\n')

            boxplot(data,k,input_folder,str(round(p_value,6)))
            barplot(data,k,input_folder,str(round(p_value,6)))

        elif data.condition.nunique() == 2:
            if normality(data,k) == False: #Check for normality
                results.append(f'Distribution of {str(k)} is not normal \n')
                p_value = ranksums(data,k)
                results.append(f'p-value of rank sums test for {k} is {str(round(p_value,6))}\n')
            elif normality(data,k) == True:
                results.append(f'Distribution of {str(k)} is normal \n')
                p_value = t_test(data,k)
                results.append(f'p-value of t-test for {k} is {str(round(p_value,6))}\n')

            boxplot(data,k,input_folder,str(round(p_value,6)))
            barplot(data,k,input_folder,str(round(p_value,6)))


        results.append(f'Sample size of {k} is {means(data,k)[0]}\n')
        results.append(f'Mean of {k} is {means(data,k)[1]}\n')
        results.append(f'Median of {k} is {means(data,k)[2]}\n')
        results.append(f'Standard deviation of {k} is {means(data,k)[3]}\n')
        results.append('\n')

    text_file = open((Path(input_folder).joinpath("Statistical test results.txt")),'w')
    text_file.writelines(results)
    text_file.close()

def ranksums(data,variable):
    """Performs a two-sided Mann-Whitney U test. Inputs a DataFrame and a string. Returns a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        p is the p-value of the resulting test.
    """
    index = data.condition.unique()
    x = data.loc[data.condition==index[0], variable]
    y = data.loc[data.condition==index[1], variable]
    p = stats.mannwhitneyu(x,y,use_continuity=True,alternative='two-sided')[1]

    return p

def kruskal(data,variable):
    """Performs a Kruskal-Wallis test. Inputs a DataFrame and a string. Returns a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        p is the p-value of the resulting test.
    """
    list_of_arrays = []
    for index in set(data.condition.unique()):
        list_of_arrays.append(data.loc[data.condition==index, variable])

    p = stats.kruskal(*list_of_arrays,nan_policy='omit')[1]
    return p

def t_test(data,variable):
    """Performs a Kruskal-Wallis test. Inputs a DataFrame and a string. Returns a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        p is the p-value of the resulting test.
    """
    list_of_arrays = []
    for index in set(data.condition.unique()):
        tmp = data.loc[data.condition==index, variable]
        # print(tmp)
        list_of_arrays.append(data.loc[data.condition==index, variable])

    p = stats.ttest_ind(*list_of_arrays,nan_policy='omit')[1]
    return p

def dunn(data,variable):
    """Performs Dunn's test. Inputs a DataFrame and a string.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        This function prints out a DataFrame of individual p-values between each conditions.
    """
    list_of_arrays = []
    conditions_list = []
    for index in set(data.condition.unique()):
        list_of_arrays.append(data.loc[data.condition==index, variable])
        conditions_list.append(index)

    p = sp.posthoc_dunn(list_of_arrays)

    index = []
    for i in range(len(list_of_arrays)): #Rename DataFrame with analysed conditions
        index.append(i+1)
    label_dict = dict(zip(index, conditions_list))
    p.rename(index=label_dict,columns=label_dict,inplace=True)

    return p.to_string()

def boxplot(data,variable,input_folder,p):
    """Generates a boxplot, with a notch at the median, and saves it as a .png. Inputs a DataFrame, two strings and a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        input_folder is the path to the input folder.
        p is the p-value of whichever test was performed before, to be displayed within the graph. 
    """

    sns.set_theme(style="ticks", palette=sns.color_palette("crest"))
    sns.boxplot(x=data['condition'],
            y=data[variable],  width=0.35, notch=True,
            data=data, showfliers =False,linewidth=1)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(vars[variable])
    plt.annotate((f'p-value : {p}'),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath(f'Boxplot {str(data["condition"].unique())} {variable}.{settings["extension"]}'),dpi=settings["dpi"])
    plt.close()
    
def barplot(data,variable,input_folder,p):
    """Generates a barplot, with SEM as error bars, and saves it as a .png. Inputs a DataFrame, two strings and a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        input_folder is the path to the input folder.
        p is the p-value of whichever test was performed before, to be displayed within the graph. 
    """
    error = []
    
    for i in data['condition'].unique():
        error.append(stats.sem(data[variable].loc[data['condition']==i],nan_policy='omit'))
    sns.barplot(y=data[variable],x=data['condition'],estimator=mean,yerr=error,ci=None,\
            error_kw={'elinewidth':2,'capsize':4,'capthick':2})
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(vars[variable])
    plt.annotate((f'p-value : {p}'),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath(f'Barplot {str(data["condition"].unique())} {variable}.{settings["extension"]}'),dpi=settings["dpi"])
    plt.close()

def violinplot(data,variable,p):
    """Generates a violinplot. Inputs a DataFrame, a string and a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        p is the p-value of whichever test was performed before, to be displayed within the graph.
    """
    
    sns.violinplot(x=data['condition'],y=data[variable],inner='box',cut=0)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(vars[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    plt.show()

def normality(data,variable):
    """Tests wether or not a sample fits a normal distribution. Inputs a DataFrame and a string. Returns a Boolean.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
    """
    if len(data) > 0:
        p = stats.normaltest(a=data[variable],nan_policy='omit')[1]
        if p > 0.05:
            return True
        else:
            return False

def means(data,variable):

    for cond in set(data.condition.unique()):
        subdata = data.loc[data.condition==cond, variable]
        sample_size = len(subdata)
        var_mean = np.nanmean(subdata)
        var_median = np.nanmedian(subdata)
        var_std = np.nanstd(subdata)
    
    return sample_size, var_mean, var_median, var_std     

if __name__ == '__main__':
    input_folder = r'/home/baptiste/Documents/test_stats_kif'
    settings = {'antero_retro':True,
                'extension':'svg',
                'dpi':300,
                'order':['WT','HET','HOM'],
                }
    statistical_analysis(settings,input_folder)