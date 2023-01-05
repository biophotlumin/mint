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
#Imports
from utils import *
from numpy.core.fromnumeric import mean
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import scikit_posthocs as sp
from pathlib import Path
import numpy as np

ylabels = {
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
            print(os.path.dirname(name))
            data = pd.read_csv(file_path,sep=csv_sniffer(file_path))

            if settings['antero_retro']==True:
                variables_antero_retro(data,input_folder)


def variables_antero_retro(data,input_folder):
    """Calls statistical tests and plotting functions for each variable of interest. Inputs a DataFrame and a string.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        input_folder is the path to the input folder.
    """
    voi = ['curvilign_velocity_antero'] #Variables of interest

    #pub_boxplot(data,'duration',input_folder,str(round((kruskal(data,'duration')),6)))
    #pub_barplot(data,'duration',input_folder,str(round((kruskal(data,'duration')),6)))
    for item in set(voi):
        if normality(data,item) == False: #Check for normality
            pub_boxplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
            pub_barplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
            #violinplot(data,item,str(round((kruskal(data,item)),6)))
            #means(data,item)

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

colorpal_dyna = {'CONTROL':'tab:blue','DYNAPYRAZOLE':'tab:cyan'}
colorpal_kif = {'WT':'darkgreen','HET':'seagreen','HOM':'lightgreen'}

def pub_boxplot(data, variable,input_folder,p):
    if 'CONTROL' in data['condition'].unique():
        colorpal = colorpal_dyna
    else:
        colorpal = colorpal_kif
    sns.set_theme(style="ticks", palette="pastel")

    
    sns.boxplot(x=data['condition'],
        y=data[variable],  width=0.35, notch=True, palette=colorpal,
        data=data, showfliers =False)

    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(ylabels[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    #plt.savefig((input_folder+"\Boxplot "+str(data['condition'].unique())+" "+variable+".svg"))
    plt.savefig(Path(input_folder).joinpath("Boxplot "+str(data['condition'].unique())+" "+variable+".svg"))
    plt.close()

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def pub_barplot(data, variable,input_folder,p):

    colorpal = colorpal_dyna
    error = []
    error.append(stats.sem(data[variable].loc[data['condition']=='DYNAPYRAZOLE'],nan_policy='omit'))
    error.append(stats.sem(data[variable].loc[data['condition']=='CONTROL'],nan_policy='omit'))

    
    plt.subplot(1,2,1)
    plt.xlabel("Condition")
    plt.ylabel(ylabels[variable])
    sns.barplot(y=data[variable],x=data['condition'],estimator=mean,yerr=error,ci=None,\
        error_kw={'elinewidth':2,'capsize':4,'capthick':2},palette=colorpal)
    sns.despine(trim=True)
    
    plt.subplot(1,2,2)
    sns.boxplot(x=data['condition'],
        y=data[variable],  width=0.35, notch=True, palette=colorpal,
        data=data, showfliers =False)
    
    
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
    #plt.savefig((input_folder+"\Barplot "+str(data['condition'].unique())+" "+variable+".svg"))
    
    plt.savefig(Path(input_folder).joinpath("Barplot "+str(data['condition'].unique())+" "+variable+".svg"))
    plt.close()

if __name__ == '__main__':
    input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/test graphs'
    settings = {'antero_retro':True}
    statistical_analysis(settings,input_folder)