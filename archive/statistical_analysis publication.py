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
        'switch_normal':'Directionality reversal per µm',
        'switch_var_STOP':'Variation of intensity in STOP phases of reversal',
        'theta_std_GO':'theta_std_GO',
        'theta_std_STOP':'theta_std_STOP',
        'pausing_time_antero':'pausing_time_antero',
        'pausing_time_retro':'pausing_time_retro',
        'pausing_time_switch':'pausing_time_switch',
        'directionality_normal':'direction normal'
        
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
            print(csv_sniffer(file_path))
            if settings['antero_retro']==True:
                variables_antero_retro(data,input_folder)
            else:
                variables(data,input_folder)

def variables(data,input_folder):
    """Calls statistical tests and plotting functions for each variable of interest. Inputs a DataFrame and a string.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        input_folder is the path to the input folder.
    """
    voi = ['curvilign_velocity','processivity','run_length','diag_size','directionality',\
            'fraction_paused','pausing_frequency','pausing_time'] #Variables of interest

    for condition in set(data.condition): #Loop for different experiments
        subset = data.loc[data.condition==condition]
        results = []
        if subset['slide'].nunique() == 2: #Applies tests for two samples
            for item in set(voi):
                if normality(data,item) == False: #Check for normality
                    results.append('Distribution of '+str(item)+' is not normal \n')
                    results.append("p-value of Wilcoxon rank sum test for "+item+" is "+str(round((ranksums(subset,item)),6))+"\n")
                    boxplot(subset,item,input_folder,str(round((ranksums(subset,item)),6)))
                    barplot(subset,item,input_folder,str(round((ranksums(subset,item)),6)))
                elif normality(data,item) == True:
                        results.append('Distribution of '+str(item)+' is normal \n')
                        #Yet to implement t-test
                        #results.append("p-value of Kruskal-Wallis test for "+item+" is "+str(round((kruskal(data,item)),6))+"\n")
                        #boxplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
                        #barplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
            
            text_file = open((Path(input_folder).joinpath("Kruskal-Wallis test results.txt")),'w')
            text_file.writelines(results)
            text_file.close()

        elif subset['slide'].nunique() > 2: #Applies tests for three or more samples
            for item in set(voi):
                if normality(data,item) == False: #Check for normality
                    results.append('Distribution of '+str(item)+' is not normal \n')
                    results.append("p-value of Kruskal-Wallis test for "+item+" is "+str(round((kruskal(subset,item)),6))+"\n")
                    dunn(subset,item)
                    boxplot(subset,item,input_folder,str(round((kruskal(subset,item)),6)))
                    barplot(subset,item,input_folder,str(round((kruskal(subset,item)),6)))
                elif normality(data,item) == True:
                        results.append('Distribution of '+str(item)+' is normal \n')
                        #Yet to implement t-test
                        #results.append("p-value of Kruskal-Wallis test for "+item+" is "+str(round((kruskal(data,item)),6))+"\n")
                        #boxplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
                        #barplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
            text_file = open((Path(input_folder).joinpath("Kruskal-Wallis test results.txt")),'w') 
            text_file.writelines(results)
            text_file.close()
        
def variables_antero_retro(data,input_folder):
    """Calls statistical tests and plotting functions for each variable of interest. Inputs a DataFrame and a string.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        input_folder is the path to the input folder.
    """
    voi = ['curvilign_velocity_antero','curvilign_velocity_retro','processivity_antero','processivity_retro','curvilign_length','switch_a_to_r',\
        'run_length_antero','run_length_retro','diag_size','directionality','duration','phase_dir_GO','switch_r_to_a',\
            'fraction_paused','pausing_frequency','pausing_time','switch','variance_GO','variance_STOP','theta_std_GO','theta_std_STOP',\
                'pausing_time_antero','pausing_time_retro','pausing_time_switch'] #Variables of interest
    #voi = ['curvilign_velocity_antero'] #Variables of interest
    results = []
    #print(data.head)
    # non_antero = data.loc[data['directionality']>0] #Tout sauf pûrement antérograde
    # non_retro = data.loc[data['directionality']<1] #Tout sauf pûrement rétrograde
    # mixed = non_antero.loc[non_antero['directionality']<1]
    # pure_retro = data.loc[data['directionality']==0]
    # left = data.loc[data['slide']=='oeil_droit']
    #print(mixed.directionality.min())
    #print(mixed.directionality.max())#data = left
    #print(len(data))
    #str(round((kruskal(non_antero,'diag_size')),6))
    #pub_boxplot(data,'duration',input_folder,str(round((kruskal(data,'duration')),6)))
    #pub_barplot(data,'duration',input_folder,str(round((kruskal(data,'duration')),6)))
    for item in set(voi):
        print(item)
        if normality(data,item) == False: #Check for normality
            results.append('Distribution of '+str(item)+' is not normal \n')
            results.append("p-value of Kruskal-Wallis test for "+item+" is "+str(round((kruskal(data,item)),6))+"\n")
            # fraction_moving(data)
            dunn(data,item)
            pub_boxplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
            pub_barplot(data,item,input_folder,str(round((kruskal(data,item)),6)))
            #violinplot(data,item,str(round((kruskal(data,item)),6)))
            means(data,item)
        elif normality(data,item) == True:
                results.append('Distribution of '+str(item)+' is normal \n')
                results.append("p-value of Kruskal-Wallis test for "+item+" is "+str(round((t_test(data,item)),6))+"\n")
                pub_boxplot(data,item,input_folder,str(round((t_test(data,item)),6)))
                pub_barplot(data,item,input_folder,str(round((t_test(data,item)),6)))
    text_file = open((Path(input_folder).joinpath("Kruskal-Wallis test results.txt")),'w')
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

    print(variable)
    print(p)

def boxplot(data,variable,input_folder,p):
    """Generates a boxplot, with a notch at the median, and saves it as a .png. Inputs a DataFrame, two strings and a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        input_folder is the path to the input folder.
        p is the p-value of whichever test was performed before, to be displayed within the graph. 
    """

    sns.set_theme(style="ticks", palette="pastel")
    sns.boxplot(x=data['condition'],
            y=data[variable],  width=0.35, notch=True, #hue=data['slide'],palette=["m", "g"],
            data=data, showfliers =False,linewidth=1)
    #sns.stripplot(x=data['condition'],
            #y=data[variable],edgecolor='gray',palette=("flare"))
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(ylabels[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath("Boxplot "+str(data['condition'].unique())+" "+variable+".png"))
    plt.close()
    
def barplot(data,variable,input_folder,p):
    """Generates a barplot, with SEM as error bars, and saves it as a .png. Inputs a DataFrame, two strings and a float.

        data is a DataFrame containing transport parameters, as defined in data_extraction.py.
        variable is the name of the variable of interest being tested.
        input_folder is the path to the input folder.
        p is the p-value of whichever test was performed before, to be displayed within the graph. 
    """
    
    error = []
    for i in set(data['condition'].unique()):
        error.append(stats.sem(data[variable].loc[data['condition']==i],nan_policy='omit'))
    sns.barplot(x=data['condition'],y=data[variable],capsize=0.02,estimator=mean,yerr=error,ci=None)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(ylabels[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath("Barplot "+str(data['condition'].unique())+" "+variable+".png"))
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
    plt.ylabel(ylabels[variable])
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

colorpal_dyna = {'CONTROL':'tab:blue','DYNAPYRAZOLE':'tab:cyan'}
colorpal_kif = {'WT':'darkgreen','HET':'seagreen','HOM':'lightgreen'}

def pub_boxplot(data, variable,input_folder,p):
    if 'CONTROL' in data['condition'].unique():
        colorpal = colorpal_dyna
    else:
        colorpal = colorpal_kif
    sns.set_theme(style="ticks", palette="pastel")
    
    if variable == 'run_length_retro' or variable =='curvilign_velocity_retro':
        sns.boxplot(x=data['condition'],
            y=(data[variable]*-1),  width=0.35, notch=True, palette=colorpal,
            data=data, showfliers =False)
    
    elif variable == 'directionality':
        #colorpal = {'CONTROL':'white','DYNAPYRAZOLE':'white'}
        colorpal = {'WT':'white','HET':'white','HOM':'white'}
        sns.boxplot(x=data['condition'],
            y=(data[variable]),  width=0.35, notch=True, palette=colorpal,
            data=data, showfliers =False)
        stripal = {'WT':'darkgreen','HET':'seagreen','HOM':'lightgreen'}
        #stripal = {'CONTROL':'tab:blue','DYNAPYRAZOLE':'tab:cyan'}
        sns.stripplot(x=data['condition'],
            y=data[variable],edgecolor='gray',palette=stripal)

    else:
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
    
def pub_barplot(data, variable,input_folder,p):
    if 'CONTROL' in data['condition'].unique():
        colorpal = colorpal_dyna
        error = []
        error.append(stats.sem(data[variable].loc[data['condition']=='CONTROL'],nan_policy='omit'))
        error.append(stats.sem(data[variable].loc[data['condition']=='DYNAPYRAZOLE'],nan_policy='omit'))
    else:
        colorpal = colorpal_kif
        error = []
        # error.append(stats.sem(data[variable].loc[data['condition']=='HET'],nan_policy='omit'))
        error.append(stats.sem(data[variable].loc[data['condition']=='HOM'],nan_policy='omit'))
        error.append(stats.sem(data[variable].loc[data['condition']=='WT'],nan_policy='omit'))

    if variable == 'run_length_retro' or variable == 'curvilign_velocity_retro':
        sns.barplot(y=(data[variable]*-1),x=data['condition'],estimator=mean,yerr=error,ci=None,\
            error_kw={'elinewidth':2,'capsize':4,'capthick':2},palette=colorpal)
    else:
        sns.barplot(y=data[variable],x=data['condition'],estimator=mean,yerr=error,ci=None,\
            error_kw={'elinewidth':2,'capsize':4,'capthick':2},palette=colorpal)
    sns.despine(trim=True)
    
    plt.xlabel("Condition")
    plt.ylabel(ylabels[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    #plt.savefig((input_folder+"\Barplot "+str(data['condition'].unique())+" "+variable+".svg"))
    plt.savefig(Path(input_folder).joinpath("Barplot "+str(data['condition'].unique())+" "+variable+".svg"))
    plt.close()

def means(data,variable):
    print(variable)
    for cond in set(data.condition.unique()):
        print(cond)
        subdata = data.loc[data.condition==cond, variable]
        #subdata.to_csv(r'/home/baptiste/dwhelper/'+variable+'_'+cond+'.csv')
        print(len(subdata))
        print(subdata.mean())
        print(np.nanmean(subdata))
        #print(subdata.median())
        print()

def fraction_moving(data):
    """
    Calculates the fraction of moving particles. Inputs a DataFrame.
    """
    list_of_arrays = []
    conditions_list = []
    for index in set(data.condition.unique()):
        arr = []
        conditions_list.append(index)
        subdata = data.loc[data.condition==index]
        for file in subdata.file.unique():
            arr.append(float((subdata.loc[data.file==file]).fraction_moving.unique()))
        list_of_arrays.append(arr)
    pdunn = sp.posthoc_dunn(list_of_arrays)
    index = []
    for i in range(len(list_of_arrays)): #Rename DataFrame with analysed conditions
        index.append(i+1)
    label_dict = dict(zip(index, conditions_list))
    pdunn.rename(index=label_dict,columns=label_dict,inplace=True)
    pkw = stats.kruskal(*list_of_arrays,nan_policy='omit')[1]

    print("Kruskal-Wallis "+str(pkw))
    print("Dunn test ")
    print(str(pdunn))
    data_per_file = data.drop_duplicates('file')
    for cond in set(data_per_file.condition.unique()):
        print(cond)
        subdata = data_per_file.loc[data_per_file.condition==cond].fraction_moving
        print(subdata.mean())
        print(subdata.median())
        print()

    sns.set_theme(style="ticks", palette="pastel")
    sns.boxplot(x=data_per_file['condition'],
            y=data_per_file['fraction_moving'],  width=0.35, notch=True,
            data=data, showfliers =False)
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel("Ratio of moving particles")
    plt.annotate(("p-value : "+str(pkw)),xy=(195,310),xycoords='figure points')
    plt.show()
    plt.close()

    print(stats.sem(data_per_file['fraction_moving'].loc[data['condition']=='WT'],nan_policy='omit'))
    error = []
    error.append(stats.sem(data_per_file['fraction_moving'].loc[data['condition']=='sans_KTP'],nan_policy='omit'))
    error.append(stats.sem(data_per_file['fraction_moving'].loc[data['condition']=='avec_KTP'],nan_policy='omit'))
    # error.append(stats.sem(data_per_file['fraction_moving'].loc[data['condition']=='HOM'],nan_policy='omit'))
    # sns.barplot(y=data_per_file['fraction_moving'],x=data_per_file['condition'],estimator=mean,yerr=stats.sem(data_per_file['fraction_moving'],nan_policy='omit'),ci=None,\
    #     error_kw={'elinewidth':2,'capsize':4,'capthick':2})
    sns.barplot(y=data_per_file['fraction_moving'],x=data_per_file['condition'],estimator=mean,yerr=error,ci=None,\
         error_kw={'elinewidth':2,'capsize':4,'capthick':2})
    sns.despine(trim=True)
    
    plt.xlabel("Condition")
    plt.ylabel('Ratio of moving particles')
    plt.annotate(("p-value : "+str(pkw)),xy=(195,310),xycoords='figure points')
    plt.show()
    plt.close()
        
colorpal_ktp = ['yellow','blue']

def pub_boxplot_ktp(data, variable,input_folder,p):
    sns.set_theme(style="ticks", palette="pastel")
    
    if variable == 'run_length_retro' or variable =='curvilign_velocity_retro':
        sns.boxplot(x=data['condition'],
            y=(data[variable]*-1),  width=0.35, notch=True, palette=colorpal_ktp,
            data=data, showfliers =False)
    else:
        sns.boxplot(x=data['condition'],
            y=data[variable],  width=0.35, notch=True, palette=colorpal_ktp,
            data=data, showfliers =False)

    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(ylabels[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath("Boxplot "+str(data['condition'].unique())+" "+variable+".svg"))
    plt.close()
    
def pub_barplot_ktp(data, variable,input_folder,p):

    error = []
    error.append(stats.sem(data[variable].loc[data['condition']=='sans_KTP'],nan_policy='omit'))
    error.append(stats.sem(data[variable].loc[data['condition']=='avec_KTP'],nan_policy='omit'))

    if variable == 'run_length_retro' or variable == 'curvilign_velocity_retro':
        sns.barplot(y=(data[variable]*-1),x=data['condition'],estimator=mean,yerr=error,ci=None,\
            error_kw={'elinewidth':2,'capsize':4,'capthick':2},palette=colorpal_ktp)
    else:
        sns.barplot(y=data[variable],x=data['condition'],estimator=mean,yerr=error,ci=None,\
            error_kw={'elinewidth':2,'capsize':4,'capthick':2},palette=colorpal_ktp)
    sns.despine(trim=True)
    
    plt.xlabel("Condition")
    plt.ylabel(ylabels[variable])
    plt.annotate(("p-value : "+p),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath("Barplot "+str(data['condition'].unique())+" "+variable+".svg"))
    plt.close()

if __name__ == '__main__':
    input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124 Results - 20221011_145854 savgol 9/WT HOM TRI'
    settings = {'antero_retro':True}
    statistical_analysis(settings,input_folder)