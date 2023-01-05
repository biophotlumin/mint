"""from numpy.lib.type_check import imag
import pims
import numpy as np
import pandas as pd
import trackpy as tp
import numba

images = pims.TiffStack_pil(r"F:\INSERM\190226-nanoKTP-kif5a.lif - Series122.tif")
print(images)

if __name__ == '__main__':
    tp.batch(images, minmass=300, diameter=9,separation=12,preprocess=False,engine='numba',processes='auto')"""
folder_pre = r'/media/baptiste/Windows/Users/LUMIN10/Desktop/Programme Baptiste/20210805 CODE REPRO VERSIONS/PRE'
folder_gurobi = r'/media/baptiste/Windows/Users/LUMIN10/Desktop/Programme Baptiste/20210805 CODE REPRO VERSIONS/GUROBI'
folder_mosek = r'/media/baptiste/Windows/Users/LUMIN10/Desktop/Programme Baptiste/20210805 CODE REPRO VERSIONS/MOSEK'
folder_scs = r'/media/baptiste/Windows/Users/LUMIN10/Desktop/Programme Baptiste/20210805 CODE REPRO VERSIONS/SCS'
import os
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import statistics
import seaborn as sns

dx = []
dy = []    
for path, subfolder, files in os.walk(folder_scs): #Scans entire folder structure for files
        for name in files:
            if name.endswith('.csv') == False:  #Check for correct file
                continue

            file_path_scs = os.path.join(path, name)
            file_path_gurobi = Path(folder_gurobi).joinpath(name.replace('SCS','GUROBI'))
            file_path_mosek = Path(folder_mosek).joinpath(name.replace('SCS','MOSEK'))
            data_scs = pd.read_csv(file_path_scs,sep=',')
            data_gurobi = pd.read_csv(file_path_gurobi,sep=',')
            data_mosek = pd.read_csv(file_path_mosek,sep=',')

            x_scs = data_scs.x
            y_scs = data_scs.y

            x_gurobi = data_gurobi.x
            y_gurobi = data_gurobi.y

            x_mosek = data_mosek.x
            y_mosek = data_mosek.y

            deltax = x_mosek.subtract(x_gurobi)
            deltay = y_mosek.subtract(y_gurobi)

            #plt.plot(deltax,deltay,lw=0.5)
            #plt.show()
            
            dx.append(float(deltax.mean(axis=0)))
            dy.append(float(deltay.mean(axis=0)))

sns.histplot(dx)
plt.show()
"""sns.histplot(dy)
plt.show()"""

print("mean")
print(statistics.mean(dx))
print(statistics.mean(dy))

print("median")
print(statistics.median(dx))
print(statistics.median(dy))

print("max")
print(max(dx))
print(max(dy))

print(len(dx))
print(len(dy))