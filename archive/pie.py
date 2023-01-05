import matplotlib.pyplot as plt

data = [52.9,41.7]
labels = ['LysoTracker +','LysoTracker -']
plt.pie(x=data,labels=labels,colors=['#87dce5ff','#2eabb8','#3274a1'],autopct='%.2f%%',textprops={'size':18})
plt.show()

"""import numpy as np
import pandas as pd
from math import sqrt
dtype=np.float64
df = pd.read_csv(r'/home/baptiste/mass.csv')
#for index,data in df.iterrows():
    #print(data)
    #print((30+(600/(np.sqrt(data-40)))))

array_mass = df['mass'].to_numpy()
list_mass = array_mass.tolist()
pd.DataFrame(array_mass).to_csv(r'/home/baptiste/sqrt.csv')
#print((30+(600/(np.sqrt(list_mass-40)))))
for i in set(array_mass):
    a = (30+(600/(sqrt(i-40))))

signal_list = []"""
"""for i in set(array_mass):
    signal_list.append((30+(600/(np.sqrt(i-40)))))"""
"""Estimated_Noise = np.array(signal_list)
print(len(Estimated_Noise))
print(len(signal_list))"""

# import os
# from pathlib import Path
# input_folder = Path(r"/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124")
# import pandas as pd
# n_raw_traj = 0
# for path, subfolder, files in os.walk(input_folder): #Scan entire folder structure for files
#     for name in files:
#         if name.endswith('_rejoined.csv') == False:  #Check for correct file
#             continue

#         #Build output file path
#         file_path = os.path.join(path, name)
#         file_folder_path = os.path.split(path)[0]
#         slide_path,slide = os.path.split(file_folder_path)


#         data = pd.read_csv(file_path,sep="\t")
#         n_raw_traj += int(data.rejoined_particle.nunique())

# print(n_raw_traj)

