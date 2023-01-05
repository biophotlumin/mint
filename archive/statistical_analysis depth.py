import pandas as pd
from utils import csv_sniffer
import seaborn as sns
from scipy import stats
import scikit_posthocs as sp
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import mean
import numpy as np

file_path = r'/media/baptiste/SHG_tracking_data/Zebrafish data/Dyna_tri_complet Results - 20220927_155531 nouveau switchs sans carrés/Trajectory average parameters.csv'

data = pd.read_csv(file_path,sep=csv_sniffer(file_path))

list_ctrl = []
list_dyna = []

# for i in set(data.animal.unique()):
#     subdata = data[data.animal==i]
#     for j in set(subdata.slide.unique()):
#         eye = subdata[subdata.slide==j]
#         varr = eye.fraction_moving
#         varr = varr.drop_duplicates()
#         if 'CONTROL' in subdata['condition'].unique():
#             list_ctrl.append(varr.std())
#         else:
#             list_dyna.append(varr.std())


# pkw = stats.kruskal(list_ctrl,list_dyna,nan_policy='omit')[1]
# sns.barplot(data=[list_ctrl,list_dyna])
# plt.show()
# print(np.nanmean(list_ctrl))
# print(np.nanmean(list_dyna))
# print(pkw)
# print(len(list_ctrl))
# print(len(list_dyna))

# for i in set(data.condition.unique()):
#     subdata = data[data.condition==i]
#     print(i)
#     for j in set(subdata.animal.unique()):
#         print(j)
#         ani = subdata[subdata.animal==j]
#         varr = ani.fraction_moving.unique().std()
#         print(varr)
#         if 'CONTROL' in subdata['condition'].unique():
#             list_ctrl.append(varr)
#         else:
#             list_dyna.append(varr)

# list_ctrl = [x / 0.06665904584532828 for x in list_ctrl]
# list_dyna = [x / 0.04554647559388804 for x in list_dyna]

# pkw = stats.kruskal(list_ctrl,list_dyna,nan_policy='omit')[1]
# sns.barplot(data=[list_ctrl,list_dyna])
# plt.show()
# plt.close
# print(np.nanmean(list_ctrl))
# print(np.nanmean(list_dyna))
# print(pkw)
# print(len(list_ctrl))
# print(len(list_dyna))
# sns.boxplot(data=[list_ctrl,list_dyna],palette=['white','white'])
# sns.stripplot(data=[list_ctrl,list_dyna])
# plt.show()

data_per_file = data.drop_duplicates('file')
ctrl = data_per_file[data_per_file.condition=='CONTROL']
dyna = data_per_file[data_per_file.condition=='DYNAPYRAZOLE']
ctrl = ctrl.fraction_moving
ctrl = ctrl/0.06665904584532828
dyna = dyna.fraction_moving
dyna = dyna/0.04554647559388804
print(len(ctrl))
print(len(dyna))
sns.set_theme(style="ticks", palette="pastel")
sns.boxplot(data = [ctrl,dyna],  width=0.35, notch=True, showfliers =False)
sns.despine(trim=True)
plt.xlabel("Condition")
plt.ylabel("Ratio of moving particles")
plt.show()
plt.close()

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
    sns.barplot(y=data_per_file['fraction_moving'],x=data_per_file['condition'],estimator=mean,yerr=stats.sem(data_per_file['fraction_moving'],nan_policy='omit'),ci=None,\
        error_kw={'elinewidth':2,'capsize':4,'capthick':2})
    sns.despine(trim=True)
    
    plt.xlabel("Condition")
    plt.ylabel('Ratio of moving particles')
    plt.annotate(("p-value : "+str(pkw)),xy=(195,310),xycoords='figure points')
    plt.show()
    plt.close()