from inspect import CO_ITERABLE_COROUTINE
import os
from re import sub
from matplotlib import cm
import matplotlib
from matplotlib.cbook import print_cycles
import pandas as pd
import numpy as np
from scipy.sparse import coo
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy

def prec(N,offset,alpha):
    return(offset+(alpha/(np.sqrt(N-40))))

input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/preloc/Jitter nKTP/tri Results - 20211104_152247 minmass100 no tophat/tri'

sd_x = []
sd_y = []
intensity = []
length = []
i = 0
for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
    for name in files:
        
        if name.endswith('.txt') == False:  #Check for correct file
            continue
        
        file_path = os.path.join(path, name)

        data = pd.read_csv(file_path,sep="\t")
        #cutoff = 600
        """if name.startswith('211022_nKTP_comparison_averaging.lif - ROI1_LA2'):
            
            cutoff = 1160
        elif name.startswith('211022_nKTP_comparison_averaging.lif - ROI2') or name.startswith('211022_nKTP_comparison_averaging.lif - ROI3') or name.startswith('211022_nKTP_comparison_averaging.lif - ROI4'):
            cutoff = 610
        else : 
            cutoff = 1200"""
        length.append(len(data.Xraw))
        if len(data) < 600:
            continue
        if (np.mean(data.Signal)/11.4) > 500:
            continue
        #print(name)
        #print(data.head)
        sd_x.append(np.std(data.Xraw))
        sd_y.append(np.std(data.Yraw))
        intensity.append(np.mean(data.Signal))
        #print(np.std(data.Yraw))
        #center_x = float(data.Xraw.mean(axis=0))
        #center_y = float(data.Yraw.mean(axis=0))

        if (np.mean(data.Signal)/11.4) > 430:
            print(name)
            print(np.std(data.Xraw)*173)
            print(np.sqrt((np.std(data.Xraw)**2)+(np.std(data.Yraw)**2)))
        
        

sigma_x = (sum(sd_x)/len(sd_x))*173
sigma_y = (sum(sd_y)/len(sd_y))*173
sigma = np.sqrt((sigma_x**2)+(sigma_y**2))
print(sigma_x)
print(sigma_y)
print(sigma)
df = pd.DataFrame(data=sd_x,columns=['sd_x'])
df['sd_x']=df['sd_x']*173
df['sd_y']=sd_y
df['sd_y']=df['sd_y']*173
df['sig'] = np.sqrt((df['sd_y']**2)+(df['sd_y']**2))
df['intensity']=intensity
df['intensity'] = df['intensity']/11.8
sdx = np.array(sd_x)
# print(np.mean(length))
# print(np.median(length))
# print(np.amax(length))
# print(np.amin(length))
# print(len(sd_x))
# print(np.mean(df.sig))
epx = pd.DataFrame((df['intensity'],df['sd_x']))
epx.to_csv('valeurs.csv')
p,c,id,_,__ = scipy.optimize.curve_fit(prec,df['intensity'],df['sd_x'],full_output=True)
print("curvefit")
print(p)
intensités_fit=np.linspace(0,500,100)
sns.scatterplot(data=df,x='intensity',y='sig')
#sns.scatterplot(data=df,x='intensity',y='sd_y',color='orange')
plt.plot(intensités_fit,prec(intensités_fit,30,600),linestyle="--",color='black')
plt.ylim([0,400])
# plt.savefig(r'/home/baptiste/PL.svg')
plt.show()
plt.close()


# i_f = r'/media/baptiste/SHG_tracking_data/Jitter nKTP/troisième série triée Results - 20211103_113050/troisième série triée/211022_nKTP_comparison_averaging.lif - ROI1_LA2.tif/211022_nKTP_comparison_averaging.lif - ROI1_LA2.tif.csv'
# data = pd.read_csv(i_f,sep="\t")
# for i in set(data.particle.unique()):
#     subdata = data[data.particle == i]
#     x = subdata.x
#     y = subdata.y
#     t = subdata.frame
#     df = pd.DataFrame(data=x,columns=['x'])
#     df['t'] = t
#     df['y'] = y
#     sns.lineplot(data=df,x='t',y='y')
#     plt.show()
#     print(stats.linregress(t,y))


# intensités=[69,124,297,530,1028,1663,2206,2898,3412]
# sigmax=np.array([60,52,46,18,25,26,14,18,12])

# xerror=0.05*sigmax

# intensités_fit=np.linspace(50,3500,1000)
# popt_x,pcov_x=scipy.optimize.curve_fit(prec,intensités,sigmax)


# plt.errorbar(intensités,sigmax,yerr=xerror,linestyle='',marker='+',c='tab:orange',label='sigma_x')
# plt.plot(intensités_fit,prec(intensités_fit,popt_x[0],popt_x[1]),linestyle="--",color='black')

# plt.legend()
# plt.xlabel('N_photons')
# plt.ylabel('sigma (nm)')
# plt.show()


from matplotlib.patches import Circle
coords = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/preloc/Jitter nKTP/tri Results - 20211104_152247 minmass100 no tophat/tri/211022_nKTP_comparison_averaging.lif - ROI8_LA2.tif/211022_nKTP_comparison_averaging.lif - ROI8_LA2.tif-8.txt',sep="\t")
print(coords.head)
coords['Xraw'] = coords['Xraw']*173
coords['Yraw'] = coords['Yraw']*173
x0 = np.mean(coords.Xraw)
y0 = np.mean(coords.Yraw)
coords['Xraw'] = (coords['Xraw']-x0)
coords['Yraw'] = (coords['Yraw']-y0)
coords['Signal'] = coords['Signal']/11.8
plt.scatter(coords.Xraw,coords.Yraw,s=2,c=coords.Signal,cmap='Reds')
plt.colorbar()
plt.xlim(-300,300)
plt.ylim(-300,300)
print(x0,y0)
draw_circle = plt.Circle((0,0), 48.45,fill=None,edgecolor='black',lw=1,ls='--')
plt.gca().add_patch(draw_circle)
plt.show()