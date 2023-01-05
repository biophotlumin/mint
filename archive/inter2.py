
from cProfile import label
from scipy import stats
import pandas as pd
from utils import csv_sniffer
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.mstats import pearsonr
from scipy.stats import ks_2samp

"""file_path = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124 Results - 20220413_13041 0.64/tri/Trajectory average parameters.csv'
#file_path = r'/media/baptiste/SHG_tracking_data/Zebrafish data/Dyna_tri_complet Results - 20220413_114440/outlier/Trajectory average parameters.csv'
data = pd.read_csv(file_path,sep=csv_sniffer(file_path))
subdata = data[data.condition=='WT']

print(stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit'))
p = stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit')[1]
print('{:f}'.format(p))
sns.set_theme(style="ticks", palette="pastel")
sns.boxplot(data=[subdata.theta_std_GO,subdata.theta_std_STOP],  width=0.3, notch=True, showfliers =False,linewidth=1)

for i, row in subdata.iterrows():
    if str(i).endswith('0'):
        plt.plot([subdata.theta_std_GO[i],subdata.theta_std_STOP[i]], lw=2)
plt.annotate(("p-value : "+str(p)),xy=(19,310),xycoords='figure points')
sns.despine(trim=True)
#plt.show()
plt.close()
n_bins = 232
counts, bins, bars = plt.hist((subdata.theta_std_GO-subdata.theta_std_STOP),bins=n_bins)
#plt.show()
plt.close()

plt.scatter(subdata.theta_std_GO,subdata.theta_std_STOP)
#plt.show()
plt.close()

print(pearsonr(subdata.theta_std_GO,subdata.theta_std_STOP))

binsize = 10

tempHET = data[data.condition=='HET']
print(int(((np.max(tempHET.theta_std_GO-tempHET.theta_std_STOP)-np.min(tempHET.theta_std_GO-tempHET.theta_std_STOP))/2)))
counts, bins, bars = plt.hist((tempHET.theta_std_GO-tempHET.theta_std_STOP),\
    bins=int(((np.max(tempHET.theta_std_GO-tempHET.theta_std_STOP)-np.min(tempHET.theta_std_GO-tempHET.theta_std_STOP))/2)),alpha=0.8,label='HET')

tempHOM = data[data.condition=='HOM']
counts, bins, bars = plt.hist((tempHOM.theta_std_GO-tempHOM.theta_std_STOP),\
    bins=int(((np.max(tempHOM.theta_std_GO-tempHOM.theta_std_STOP)-np.min(tempHOM.theta_std_GO-tempHOM.theta_std_STOP))/2)),alpha=0.8,label='HOM')

tempWT = data[data.condition=='WT']
counts, bins, bars = plt.hist((tempWT.theta_std_GO-tempWT.theta_std_STOP),\
    bins=int(((np.max(tempWT.theta_std_GO-tempWT.theta_std_STOP)-np.min(tempWT.theta_std_GO-tempWT.theta_std_STOP))/2)),alpha=0.8,label='WT')

print(counts)
print(bins)
print(len(tempWT),len(tempHET),len(tempHOM))
plt.legend()
plt.show()

p = stats.kruskal((tempWT.theta_std_GO-tempWT.theta_std_STOP),(tempHET.theta_std_GO-tempHET.theta_std_STOP),(tempHOM.theta_std_GO-tempHOM.theta_std_STOP),nan_policy='omit')[1]
print('{:f}'.format(p))"""

file_path = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/WT HOM/theta/Trajectory average parameters.csv'
#file_path = r'/media/baptiste/SHG_tracking_data/Zebrafish data/Dyna_tri_complet Results - 20220413_114440/outlier/Trajectory average parameters.csv'
data = pd.read_csv(file_path,sep=csv_sniffer(file_path))
subdata = data[data.condition=='HOM']

print(stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit'))
p = stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit')[1]
print('{:f}'.format(p))
sns.set_theme(style="ticks", palette="pastel")
sns.boxplot(data=[subdata.theta_std_GO,subdata.theta_std_STOP],  width=0.3, notch=True, showfliers =False,linewidth=1)

for i, row in subdata.iterrows():
    if str(i).endswith('0'):
        plt.plot([subdata.theta_std_GO[i],subdata.theta_std_STOP[i]], lw=2)
plt.annotate(("p-value : "+str(p)),xy=(19,310),xycoords='figure points')
sns.despine(trim=True)
#plt.show()
plt.close()
n_bins = 100
counts, bins, bars = plt.hist((subdata.theta_std_GO-subdata.theta_std_STOP),bins=n_bins)
#plt.show()
plt.close()

plt.scatter(subdata.theta_std_GO,subdata.theta_std_STOP)
#plt.show()
plt.close()

print(pearsonr(subdata.theta_std_GO,subdata.theta_std_STOP))


"""tempHET = data[data.condition=='HET']
print(int(((np.max(tempHET.theta_std_GO-tempHET.theta_std_STOP)-np.min(tempHET.theta_std_GO-tempHET.theta_std_STOP))/1)))
counts, bins, bars = plt.hist((tempHET.theta_std_GO-tempHET.theta_std_STOP),\
    bins=int(((np.max(tempHET.theta_std_GO-tempHET.theta_std_STOP)-np.min(tempHET.theta_std_GO-tempHET.theta_std_STOP))/1)),alpha=0.8,label='HET')"""

tempHOM = data[data.condition=='HOM']
counts, bins, bars = plt.hist((tempHOM.theta_std_GO-tempHOM.theta_std_STOP),\
    bins=int(((np.max(tempHOM.theta_std_GO-tempHOM.theta_std_STOP)-np.min(tempHOM.theta_std_GO-tempHOM.theta_std_STOP))/1)),alpha=0.8,label='HOM')

tempWT = data[data.condition=='WT']
counts, bins, bars = plt.hist((tempWT.theta_std_GO-tempWT.theta_std_STOP),\
    bins=int(((np.max(tempWT.theta_std_GO-tempWT.theta_std_STOP)-np.min(tempWT.theta_std_GO-tempWT.theta_std_STOP))/1)),alpha=0.8,label='WT')


#print(len(tempWT),len(tempHET),len(tempHOM))
plt.legend()
#plt.show()
print("kruskal")
p = stats.kruskal((tempWT.theta_std_GO-tempWT.theta_std_STOP),(tempHOM.theta_std_GO-tempHOM.theta_std_STOP),nan_policy='omit')[1]
print(p)
print('{:f}'.format(p))

prob = np.cumsum((tempWT.theta_std_GO-tempWT.theta_std_STOP))
#print(prob.shape)
#print(np.max(prob),np.min(prob))
#plt.plot(prob)
#plt.show()

plt.close()
data_sorted = np.sort((tempWT.theta_std_GO-tempWT.theta_std_STOP))
p2 = 1. * np.arange(len((tempWT.theta_std_GO-tempWT.theta_std_STOP))) / (len((tempWT.theta_std_GO-tempWT.theta_std_STOP)) - 1)
pWT = p2
plt.plot(data_sorted, p2,label='WT')

"""data_sorted = np.sort((tempHET.theta_std_GO-tempHET.theta_std_STOP))
p2 = 1. * np.arange(len((tempHET.theta_std_GO-tempHET.theta_std_STOP))) / (len((tempHET.theta_std_GO-tempHET.theta_std_STOP)) - 1)
pHET = p2
plt.plot(data_sorted, p2,label='HET')"""

data_sorted = np.sort((tempHOM.theta_std_GO-tempHOM.theta_std_STOP))
p2 = 1. * np.arange(len((tempHOM.theta_std_GO-tempHOM.theta_std_STOP))) / (len((tempHOM.theta_std_GO-tempHOM.theta_std_STOP)) - 1)
pHOM = p2
plt.plot(data_sorted, p2,label='HOM')
plt.legend()
plt.show()

print(ks_2samp(pWT,pHOM))
print(ks_2samp((tempHOM.theta_std_GO-tempHOM.theta_std_STOP),(tempWT.theta_std_GO-tempWT.theta_std_STOP)))

colorpal = {'WT':'darkgreen','HOM':'lightgreen'}
colorpal = {'CONTROL':'tab:blue','DYNAPYRAZOLE':'tab:cyan'}
print(np.mean(subdata['theta_std_STOP']))
print(np.median(subdata['theta_std_STOP']))
p = stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit')[1]
print(round(p,35))
sns.boxplot(width=0.35, notch=True,
            data=[subdata['theta_std_GO'],subdata['theta_std_STOP']], showfliers =False)

sns.despine(trim=True)
plt.xlabel("Condition")
#plt.ylabel(ylabels[variable])
plt.annotate(("p-value : "+ str(stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit')[1])),xy=(195,310),xycoords='figure points')
#plt.savefig((input_folder+"\Boxplot "+str(data['condition'].unique())+" "+variable+".svg"))
#plt.savefig(Path(input_folder).joinpath("Boxplot "+str(data['condition'].unique())+" "+variable+".svg"))
plt.show()
plt.close()

"""NA=0.95
am = np.arcsin(NA/1.3)
Atot = 8*np.pi/3

C = (4*np.pi/3+2*np.pi/3*((np.cos(am))**3-3*np.cos(am)))/Atot
B = np.pi*np.cos(am)*(np.sin(am))**2/Atot

print(B,C)

def eff(theta):
    return C+B*(np.sin(theta)**2)

theta = np.linspace(0,(90*(np.pi/180)),100)

e = eff(theta)
print(eff(0))
print(eff(60*(np.pi/180)))
theta = theta*(180/np.pi)
plt.grid(color='lightgray', linestyle='-', linewidth=0.5)
plt.plot(theta,e)

plt.show()"""

"""input_file = r'/home/baptiste/Documents/210330_nktp_kif5a.lif - Series111.tif_rejoined.csv'
df = pd.read_csv(input_file,sep='\t')
df.x = df.x*0.173
df.y = df.y*0.173
plt.gca().set_aspect('equal', adjustable='box')
for i in set(df.rejoined_particle.unique()):
    subdf = df[df.rejoined_particle==i]
    plt.plot(subdf.x,subdf.y)
plt.show()"""