from cmath import exp
from turtle import shape
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from scipy import optimize

# data = pd.read_csv("/media/baptiste/SHG_tracking_data/Zebrafish data/phase/210330_nktp_kif5a.lif - Series007.tif_rejoined.csv_4718.csv")

# stop = data.loc[data.phase==0]
# go = data.loc[data.phase==2]
# #plt.subplot(1, 2, 1)
# #plt.axis('scaled')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.scatter(go.x,go.y,lw=0.1,c='blue')
# plt.scatter(stop.x,stop.y,lw=0.1,c='red')

# #plt.subplot(1, 2, 2)
# #plt.scatter(data.x,data.y,c=data.conf,cmap='Reds')
# #plt.colorbar()

# plt.show()
# conf = pd.DataFrame()
# clist = []
# input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/phase'
# for path, subfolder, files in os.walk(input_folder): #Scan entire folder structure for files
#         for name in files:
#             file_path = os.path.join(path, name)
#             data = pd.read_csv(file_path)
#             #print(data.conf.values.tolist())
#             clist += data.conf.values.tolist()

#             #conf = pd.concat((conf,data.conf),axis=1)
#             #print(conf.shape)



# #data = pd.read_csv("/media/baptiste/SHG_tracking_data/Zebrafish data/phase/210330_nktp_kif5a.lif - Series086.tif_rejoined.csv_3118.csv")
# print(len(clist))
# #sns.histplot(clist,bins=500)
# #plt.show()
# n_bins = 15
# #counts, bins, bars = plt.hist(clist,bins=n_bins)
# counts, bins, bars = plt.hist(clist,bins=n_bins)
# plt.show()
# plt.close()
# bins = bins[:-1]
# #plt.scatter(np.log(bins),np.log(counts))
# #plt.scatter((bins),(counts))
# #plt.show()



#Am = counts[0]
#AM = counts[n_bins-1]
#print(Am)
#print(AM)
#bins = bins[:-1]
#pd.DataFrame({'bins':bins,'counts':counts}).to_csv('/media/baptiste/SHG_tracking_data/Zebrafish data/Histogramme rconf sans carrés.csv')

"""print(bins[n_bins])
def fonc(R,Am,Rm,AM,R_M):
    return Am*np.exp((-1*R)/Rm)+AM*np.exp((R-1)/R_M)

#f = lambda R,Rm,R_M : Am*np.exp(-R/Rm)+AM*np.exp((R-1)/R_M)
bins = bins[:-1]
pd.DataFrame({'bins':bins,'counts':counts}).to_csv('/media/baptiste/SHG_tracking_data/Zebrafish data/Histogramme rconf.csv')
val,cov=optimize.curve_fit(f=fonc,xdata=bins,ydata=counts,bounds=([0,0,0,0],[500000,1,500000,1]))
print(val)
plt.plot(bins, fonc(bins, *val))
plt.show()"""

"""data = pd.read_csv("/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124 Results - 20220331_190108/Per phase parameters.csv",sep='\t')
sns.histplot(np.abs(data.curvilign_velocity))
plt.show()"""

data = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124 Results - 20220413_130415 0.64/Per phase parameters.csv',sep='\t')
go = data[data.phase==2]
print(len(data))
print(len(go))
go[go.curvilign_velocity < 0] = (go[go.curvilign_velocity < 0])*-1


wt = go[go.condition == 'WT']
hom = go[go.condition == 'HOM']
bin_wt = (2*(np.subtract(*np.percentile(wt.curvilign_velocity, [75, 25]))))/(len(wt.curvilign_velocity)**(1/3))
bin_hom = (2*(np.subtract(*np.percentile(hom.curvilign_velocity, [75, 25]))))/(len(hom.curvilign_velocity)**(1/3))
print(bin_wt,bin_hom)
plt.hist((wt.curvilign_velocity),bins=int(bin_wt),histtype="bar",alpha=0.5)
plt.hist((hom.curvilign_velocity),bins=int(bin_hom),histtype="bar",alpha=0.5)
plt.legend(('wt','hom'))
plt.show()



