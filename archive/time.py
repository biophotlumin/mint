duration = 4000


if duration < 3600:
    print((duration)/60)
else:
    print(duration//3600)
    print(duration%3600)
    print('%dh%s' % (int(duration//3600),f'{int((duration%3600)/60):02d}'))
    print('%dh%s' % (int(duration//3600),int((duration%3600)/60)))

import os
from posixpath import dirname
from re import L
import pandas as pd
import numpy as np


input_folder = r'/media/baptiste/SHG_tracking_data/Jitter nKTP/tri Results - 20211104_152247 minmass100 no tophat/tri'

sd_x = []
sd_y = []
intensity = []
length = []
i = 0
lt = []
for path, subfolder, files in os.walk(input_folder): #Scans entire folder structure for files
    for name in files:
        
        if name.endswith('.txt') == False:  #Check for correct file
            continue
        
        file_path = os.path.join(path, name)
        #print(dirname(file_path))
        data = pd.read_csv(file_path,sep="\t")
        length.append(len(data.Xraw))
        if len(data) < 600:
            continue
        if (np.mean(data.Signal)/11.4) > 500:
            continue
        lt.append(dirname(file_path))
x = np.array(lt)
print(len(np.unique(x)))

