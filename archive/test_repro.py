"""import trackpy as tp
import numpy as np
data = []
data = np.load('frames.csv.npy',allow_pickle=True)
if __name__ == '__main__':
    raw_coordinates = tp.batch(data, minmass=300, diameter=9, \
                separation=12,preprocess=False,engine='numba',processes='auto')"""
"""from imageio import volread
import numpy as np
frames = volread(r'F:\INSERM\190226-nanoKTP-kif5a.lif - Series121.tif')
np.save('frames.csv',frames)"""
"""import pandas as pd
import pathlib
df = pd.DataFrame()
input_folder = r'/home/baptiste'
test = pathlib.Path(input_folder).joinpath('Documents','GitHub')
df.to_csv(test.joinpath('df.csv'),sep='\t')"""

"""import numpy as np
import pandas as pd
subdata = pd.DataFrame([(3,5)],columns=['x','y'])
print(subdata)
array_x = subdata['x'].to_numpy()
array_x = array_x[:, np.newaxis]
array_y = subdata['y'].to_numpy()
array_y = array_y[:, np.newaxis]
array = np.concatenate((array_x,array_y),axis=1)
print(array)
    
xm = array_x[:,np.newaxis]
ym = array_y[:,np.newaxis]"""
"""import os
import shutil
from datetime import datetime
from pathlib import Path

def folder_structure_creation_old(input_folder):
          
    identifier = " - Results - "+str(datetime.now().strftime('%Y%m%d_%H%M%S')) #Create a string that is added to output file paths 
    output_folder = (input_folder+identifier) #Set path for root output folder
    root_input_folder = os.path.dirname(input_folder) #Get folder directly before root input folder
    if os.path.dirname(root_input_folder) == root_input_folder: #Prevents conflict in case input folder is placed at the root of a drive
        root_input_folder = input_folder
    return output_folder,identifier,root_input_folder

def folder_structure_creation_new(input_folder):
            
    identifier = " - Results - "+str(datetime.now().strftime('%Y%m%d_%H%M%S')) #Create a string that is added to output file paths 

    output_folder = Path(input_folder.parent).joinpath(input_folder.name + identifier) #Set path for root output folder
    root_input_folder = os.path.dirname(input_folder) #Get folder directly before root input folder

    #if os.path.dirname(root_input_folder) == root_input_folder: #Prevents conflict in case input folder is placed at the root of a drive
        #root_input_folder = input_folder

    return output_folder,identifier,root_input_folder

input_folder_old = r'F:\INSERM'

input_folder_new = Path(r"F:\INSERM")

#print(folder_structure_creation_old(input_folder_old)[0])

#print(folder_structure_creation_new(input_folder_new)[0])

old_output_folder,old_identifier,old_root_input_folder = folder_structure_creation_old(input_folder_old)

new_output_folder,new_identifier,new_root_input_folder = folder_structure_creation_new(input_folder_new)

print(old_output_folder)

print(old_output_folder+input_folder_old.replace(old_root_input_folder,''))

print(new_output_folder)

print(Path(new_output_folder).joinpath(input_folder_new.name))

data_input_folder = Path(new_output_folder).joinpath(input_folder_new.name)
data_input_folder = data_input_folder.joinpath("aaaa")
print(data_input_folder)

print(len(data_input_folder.parents))
print(data_input_folder.parents[0])"""

"""from os import sep
import numpy as np
import pandas as pd
from sklearn import linear_model
from scipy import optimize"""

"""x = np.array([1,2,3])
y = np.array([3,4,5])
print(x,y)
data = pd.read_csv("/media/baptiste/SHG_tracking_data/INSERM - Results - 20210916_155034/INSERM/190226-nanoKTP-kif5a.lif - Series122.tif/190226-nanoKTP-kif5a.lif - Series122.tif.csv",sep="\t")
x = data.x
y = data.y
xm = x[:,np.newaxis]
ym = y[:,np.newaxis]
print(x.shape)
print(xm.shape)
parameters = {'len_cutoff':30,'threshold_poly3':1.4}
def f(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d"""



import pandas as pd
import os
import statistics
import numpy as np

def statistical_analysis(input_folder):

    for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
        for name in files:
            if name.endswith('rajectory average parameters.csv') == False:  #Check for correct file extension
                continue #Skips to next file if not correct .csv 
            file_path = os.path.join(path, name)
            print(os.path.dirname(name))
            data = pd.read_csv(file_path,sep=',')
            #data = data.dropna(axis=0,how='all')
            print('n traj')
            print(len(data))
            #print(data.shape)
            dyna = data.loc[data.condition=='WT']
            print(len(dyna))
            pure_retro = dyna.loc[dyna['directionality']==1]
            pure_antero = dyna.loc[dyna['directionality']==0] 
            mixed1 = dyna.loc[dyna['directionality']>0]
            mixed2 =  mixed1.loc[mixed1['directionality']<1]
            mixed_antero = mixed2.loc[mixed2['directionality']<0.5]
            mixed_retro = mixed2.loc[mixed2['directionality']>=0.5]
            print('pure antero')
            print(len(pure_antero))
            print('pure retro')
            print(len(pure_retro))
            #print('mixed antero')
            #print(len(mixed_antero))
            #print('mixed retro')
            #print(len(mixed_retro))
            print('mixed')
            print(len(mixed2))
            print('n traj cond')
            print(len(dyna))
            print('n file cond')
            print(len(dyna['file'].unique()))
            print('n file')
            print(len(data['file'].unique()))

            t = data.loc[data.file.str.startswith('200128')]
            #print('190205')
            #print(t['animal'].nunique())
            u = data.loc[data.file.str.startswith('200114')]
            #print('210330')
            #print(u['animal'].nunique())
            #v = data.loc[data.file.str.startswith('190226')]
            #print('190226')
            #print(v['animal'].nunique())
            nlarv = t['animal'].nunique() + u['animal'].nunique() #+ v['animal'].nunique()
            print('true n larva')
            print(nlarv)
            tc = t.loc[t.condition=='HET']
            #vc = v.loc[v.condition=='HOM']
            uc = u.loc[u.condition=='HET']
            print('true n larva cond')
            print(tc['animal'].nunique() + uc['animal'].nunique())# + vc['animal'].nunique())
            
            print('n larva')
            print(data['animal'].nunique())
            print('n larva cond')
            print(dyna['animal'].nunique())
            
            print(data.condition.unique())
            print('no STOPs')
            no_stop = dyna.loc[data.pausing_time == 0.]
            print(len(no_stop))
            """n= []
            for i in set(data.file.unique()):
                n.append(data.loc[data.file==i].trajectory.nunique())
            print(n)
            print(statistics.mean(n))"""
            print(mixed2.animal.nunique())
            print(len(mixed2))
            mc = mixed2.loc[mixed2.condition=='CONTROL']
            print(mc.animal.nunique())
            print(len(mc))
            print(len(dyna.loc[dyna['directionality']==0]))
            print(len(dyna.loc[dyna['directionality']>0 ]))
            
            
input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/WT HOM/mixed'
statistical_analysis(input_folder)

"""from pathlib import Path
import glob
import os

input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20210921_182942 line average 4/124'
file = '190226-nanoKTP-kif5a.lif - Series025.tif_rejoined.csv'
t = input_folder + "/" + file
print(glob.glob(t))"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from trajectory_calculations import minimization

"""input_folder = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20211103_105350 no LA r40 redressement/124'

ref = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20211103_105350 no LA r40 redressement/124 Results - 20211104_025722/Trajectory average parameters.csv',sep='\t')

parameters = {'px':0.175,'sigma':175}

for path, subfolder, files in os.walk(input_folder): #Scan entire folder structure for files
    for name in files:
        if name.endswith('_rejoined.csv') == False:  #Check for correct file
            continue

        #Build output file path
        file_path = os.path.join(path, name)
        print(type(ref['file'].unique()))
        
        if name in ref['file'].unique():
            f = pd.read_csv(file_path,sep='\t')
            subref = ref.loc[ref['file']==name]
            for traj in set(subref['trajectory'].unique()):
                subf = f.loc[f['rejoined_particle']==traj]
                #subf = minimization(subf,parameters)
                bbox_props = dict(boxstyle="round", fc="w", ec="w")
                arrow_props = dict(arrowstyle="simple")
                x = subf.x
                y = subf.y
                plt.plot(x,y,lw=0.5)
                x0,y0 = x.iloc[0],y.iloc[0]
                plt.annotate(str(traj),xy=(x0,y0),xytext=(x0+5,y0+5),bbox=bbox_props,arrowprops=arrow_props, size = 4)
        plt.savefig(Path(r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20211103_105350 no LA r40 redressement/ref').joinpath(""+str(name)+".png"),dpi=300)
        plt.close()"""

"""data = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/paper2 Results - 20211029_094442/paper2/210330_nktp_kif5a.lif - Series057.tif/210330_nktp_kif5a.lif - Series057.tif_rejoined.csv',sep="\t")
sub = data.loc[data.rejoined_particle==450]
x = sub.x
y = sub.y
plt.plot(x,y,lw=0.5)
plt.show()"""