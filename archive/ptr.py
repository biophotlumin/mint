
from cmath import nan
import pandas as pd
import matplotlib.pyplot as plt
import os
from utils import csv_sniffer
import numpy as np

# df = pd.read_csv(r"/media/baptiste/SHG_tracking_data/Zebrafish data/Dyna_tri_complet Results - 20220927_155531 nouveau switchs sans carrés/Trajectory average parameters.csv",sep="\t")
# print(len(df))
# print(df.head)

# dd = df[df.pausing_time_retro.isna() & df.pausing_time_antero.isna() & df.pausing_time_switch.isna()]

# print(dd.pausing_time.mean())

# # dr = df[df.pausing_time_switch.notna()==True & df.condition=='WT']
# # print(len(dr))

# da = df[df.pausing_time_antero.notna()==True]
# print(da)
# da = da[da.condition=='CONTROL']
# print(da)
# from pathlib import Path
# bureau = Path('/home/baptiste/Bureau')
# pre = pd.read_csv(bureau.joinpath('pre_tri.csv'))
# post = pd.read_csv(bureau.joinpath('post_tri.csv'))
# # print(pre.head)
# # print(post.head)

# merge = pd.merge(post,pre,how='inner',on='intensity_STOP')
# merge.to_csv(bureau.joinpath('merge.csv'))

# moving = pd.read_csv(r'/home/baptiste/theta/14480.csv',sep='\t')
# fixed = pd.read_csv(r'/home/baptiste/theta/roi9_2.txt',sep='\t')
# print(moving.head)
# print(fixed.head)

# plt.plot(moving.mass)
# plt.plot(fixed.Signal)
# plt.show()

# candidats = []
# trajs = pd.read_csv(r'/home/baptiste/theta/Trajectory average parameters.csv',sep='\t')
# for f in set(trajs.file.unique()):
#     files = trajs[trajs.file==f]
#     for t in set(files.trajectory.unique()):
#         tr = files[files.trajectory==t]
#         print(len(tr))
#         if len(tr) > 1100 and len(tr) < 1400:
#             candidats.append(f,t)
# print(candidats)

# input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124'
# candidats = []
# verif = pd.read_csv(r'/home/baptiste/theta/Trajectory average parameters.csv')
# for path, subfolder, files in os.walk(input_folder): #Scan entire folder structure for files
#     for name in files:
#         if name.endswith('_rejoined.csv') == False:  #Check for correct file
#             continue

#         #Build output file path
#         file_path = os.path.join(path, name)
#         file_folder_path = os.path.split(path)[0]
#         slide_path,slide = os.path.split(file_folder_path) 

#         data = pd.read_csv(file_path,sep=csv_sniffer(file_path))
#         for t in set(data.rejoined_particle.unique()):
#             subdata = data[data.rejoined_particle==t]
#             if len(subdata) > 1200 and len(subdata) < 1300:
#                 candidats.append((files[2],t))
# print(len(candidats))

# print(candidats)

moving = pd.read_csv(r'/home/baptiste/theta/14480.csv',sep='\t')
moving = pd.read_csv(r'/home/baptiste/theta/Series042.csv',sep='\t')
moving = moving[moving.rejoined_particle==1028]
moving = moving.reset_index()
fixed = pd.read_csv(r'/home/baptiste/theta/roi9_2.txt',sep='\t')

# fixed = pd.read_csv(r'/home/baptiste/theta/roi2_0.txt',sep='\t')

xticks = np.arange(0,len(moving),100)
xlabels = xticks*0.05
xlabels = xlabels.astype(int)
# plt.xticks(ticks=xticks,labels=xlabels)
# plt.plot((moving.mass/11.8),c="tab:blue")
# plt.plot((fixed.Signal/11.8),c="tab:orange")
# plt.show()
# plt.close()
# print(len(moving))
# plt.scatter(xticks,(moving.mass/11.8))
# plt.show()


from scipy import  signal

# NA=0.95 #Numerical aperture of the objective
# am = np.arcsin(NA/1.3) #1.3 is the refractive index of water, change accordingly
# Atot = 8*np.pi/3

# C = (4*np.pi/3+2*np.pi/3*((np.cos(am))**3-3*np.cos(am)))/Atot
# B = np.pi*np.cos(am)*(np.sin(am))**2/Atot

# N_ratio = moving.mass/max(moving.mass)*(B+C)
# n2 = []
# # for i in range(len(moving)-1):
# #     if i==0:
# #         n2.append(np.mean([moving.mass[i],moving.mass[i+1]]))
# #     elif i==(len(moving)-1):
# #         n2.append(np.mean([moving.mass[i-1],moving.mass[i]]))
# #     else:
# #         n2.append(np.mean([moving.mass[i-1],moving.mass[i],moving.mass[i+1]]))

# for i in range(len(moving)-1):
#     if i==0:
#         n2.append(np.mean([moving.mass[i],moving.mass[i+1],moving.mass[i+2]]))
#     elif i==1:
#         n2.append(np.mean([moving.mass[i-1],moving.mass[i],moving.mass[i+1],moving.mass[i+2]]))
#     elif i==(len(moving)-1):
#         n2.append(np.mean([moving.mass[i-2],moving.mass[i-1],moving.mass[i]]))
#     elif i==(len(moving)-2):
#         n2.append(np.mean([moving.mass[i-2],moving.mass[i-1],moving.mass[i],moving.mass[i+1]]))
#     else:
#         n2.append(np.mean([moving.mass[i-2],moving.mass[i-1],moving.mass[i],moving.mass[i+1],moving.mass[i+2]]))


# N2_ratio = n2/max(n2)*(B+C)
# N2_ratio2 = []
# for i in range(len(n2)-1):
#     N2_ratio2.append((n2[i]-min(n2))/(max(n2)-min(n2))*(B+C))
# # plt.xticks(ticks=xticks,labels=xlabels)
# # plt.plot((N2_ratio2),c="tab:blue")
# # plt.show()
# # plt.close()
# thetalist = []
# # for n in N_ratio:
# #     if n>C:
# #         thetalist.append(np.arcsin(np.sqrt((n-C)/B)))
# #     else:
# #         thetalist.append(0)
# import math

# savgol = signal.savgol_filter(moving.mass,window_length=25,polyorder=3,mode="nearest")

# savgol = savgol/11.8
# plt.xticks(ticks=xticks,labels=xlabels)
# plt.plot(moving.mass/11.8,c='tab:blue')
# plt.plot((savgol),c="tab:orange")
# plt.show()
# plt.close()
# N_savgol = savgol/max(savgol)*(B+C)
# # N_savgol = []
# # for i in range(len(savgol)-1):
# #     N_savgol.append((n2[i]-min(n2))/(max(n2)-min(n2))*(B+C))
# for n in N_savgol:
#     # if math.isnan(np.arcsin(np.sqrt((n-C)/B))) == True:
#     #     print(n)
#     if n > C:
#         thetalist.append(np.arcsin(np.sqrt((n-C)/B)))
#     # else:
#     #     thetalist.append(0)
#     else:
#         thetalist.append(np.arcsin(np.sqrt((n-C)/B*-1))*-1)
#     # thetalist.append(np.arcsin(np.emath.sqrt((n-C)/B)))

# theta = np.array(thetalist)*180/np.pi
# print(theta)
# plt.xticks(ticks=xticks,labels=xlabels)
# plt.plot((theta),c="tab:blue")
# plt.show()
# plt.close()

# plt.xticks(ticks=xticks,labels=xlabels)
# plt.plot(fixed.Signal/11.8)
# sav_fixed = signal.savgol_filter(fixed.Signal,window_length=25,polyorder=3,mode="nearest")
# sav_fixed = sav_fixed/11.8
# plt.plot(sav_fixed)
# plt.show()
# plt.close()

# savgol = signal.savgol_filter(moving.mass,window_length=11,polyorder=3,mode="nearest")

# savgol = savgol/11.8
# plt.xticks(ticks=xticks,labels=xlabels)
# plt.plot(moving.mass/11.8,c='tab:blue')
# plt.plot((savgol),c="tab:orange")
# plt.show()
# plt.close()

# thetalist = []

# # savgol = moving.mass/11.8
# for n in savgol:
#     thetalist.append(np.arcsin(np.sqrt((n-np.min(savgol))/(np.max(savgol)-np.min(savgol)))))

# print(np.min(savgol))
# print(np.max(savgol))
# theta = np.array(thetalist)*180/np.pi
# print(theta)
# plt.xticks(ticks=xticks,labels=xlabels)
# yticks = np.arange(0,100,10)
# plt.yticks(ticks=yticks)
# print(len(theta))
# plt.plot((theta),c="tab:green")
# # plt.scatter(x=np.arange(0,len(theta)),y=theta,c="tab:blue",s=1)
# plt.show()
# plt.close()

# from scipy.fftpack import *
# from scipy.interpolate import interp1d

# timestep = 0.05 #en secondes ?
# N = len(fixed)
# freq = fftfreq(N,timestep)
# TF2 = abs(fft(fixed.Signal/11.8))*timestep
# DSP = TF2**2 *N*timestep
# DSP_norm = DSP/DSP[0]

# plt.figure(figsize=(12,5))
# plt.subplot(121)
# plt.plot(freq[:500],TF2[:500])
# plt.xlabel("fréquence en Hz")
# plt.ylabel("Transformée de Fourier du signal")
# plt.subplot(122)
# plt.plot(freq[:50],DSP_norm[:50])
# plt.xlabel("fréquence en Hz")
# plt.ylabel("Densité Spectrale de Puissance normalisée")
# plt.show()

print(np.std(fixed.Signal/11.8))
print(np.mean(fixed.Signal/11.8))

# from numpy.random import randn
# T = 60 #durée en secondes
# N = T*20 + 1 #Nbs de points à 20Hz
# t = np.linspace(0,T,N) #temps échantillonné à 20Hz
# intensity_base = 200
# intensity_mes = intensity_base + np.sqrt(intensity_base)*randn(N)
# sav_gol = signal.savgol_filter(intensity_mes,window_length=25,polyorder=3,mode="nearest")

# plt.plot(t,intensity_mes)
# plt.plot(t,sav_gol)
# plt.ylim([0,300])
# plt.show()

from scipy import stats
import seaborn as sns

file_path = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124 Results - 20221010_170113 nouveau théta savgol 11/WT HOM/tri/Trajectory average parameters.csv'
data = pd.read_csv(file_path,sep=csv_sniffer(file_path))
subdata = data[data.condition=='DYNAPYRAZOLE']

# print(stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit'))
# p = stats.kruskal(subdata.theta_std_GO,subdata.theta_std_STOP,nan_policy='omit')[1]
# print('{:f}'.format(p))
# print("GO ",np.mean(subdata.theta_std_GO))
# print("STOP ",np.mean(subdata.theta_std_STOP))
# print(((np.mean(subdata.theta_std_STOP)-np.mean(subdata.theta_std_GO))/np.mean(subdata.theta_std_GO))*100)
# sns.set_theme(style="ticks", palette="pastel")
# sns.boxplot(data=[subdata.theta_std_GO,subdata.theta_std_STOP],  width=0.3, notch=True, showfliers =False,linewidth=1)
# plt.show()


list_of_arrays = []
for index in set(data.condition.unique()):
    list_of_arrays.append(data.loc[data.condition==index, 'theta_std_GO'])
print('kruskal')
p = stats.kruskal(*list_of_arrays,nan_policy='omit')[1]
subdata_CTRL = data[data.condition=='WT']
subdata_DYNA = data[data.condition=='HOM']
print(p)
print("CTRL ",np.mean(subdata_CTRL.theta_std_GO))
print("DYNA ",np.mean(subdata_DYNA.theta_std_GO))
print(((np.mean(subdata_DYNA.theta_std_GO)-np.mean(subdata_CTRL.theta_std_GO))/np.mean(subdata_CTRL.theta_std_GO))*100)

list_of_arrays = []
for index in set(data.condition.unique()):
    list_of_arrays.append(data.loc[data.condition==index, 'theta_std_GO'])

p = stats.ttest_ind(*list_of_arrays,nan_policy='omit')[1]
print("t-test")
print(p)