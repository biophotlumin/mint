import os
import pandas as pd
from pathlib import Path

input_folder = r'/media/baptiste/Windows/Users/LUMIN10/Documents/wfh/Dyna_tri Results - 20211105_091209/Dyna_tri/DYNAPYRAZOLE'
output_folder = r'/media/baptiste/Windows/Users/LUMIN10/Documents/wfh/Dyna_tri Results - 20211105_091209/deephl/CONDITION'

for path, subfolder, files in os.walk(input_folder): #Scan entire folder structure for files
        for name in files:
            if name.endswith('_rejoined.csv') == False:  #Check for correct file
                continue

            #Build output file path
            file_path = os.path.join(path, name)
            data = pd.read_csv(file_path,sep='\t')
            
            for traj in set(data.particle.unique()):
                df_traj = data[data.particle == traj]
                df_traj.drop(['mass','size','ecc','signal','raw_mass','ep','frame','particle','rejoined_particle','N','SNR','feet','n_particles'], axis='columns', inplace=True)
                df_traj['time'] = df_traj.index
                df_traj.reset_index(drop=True, inplace=True)
                df_traj.time = df_traj.time - int(df_traj.time.min())
                df_traj.time = df_traj.time * 0.05
                df_traj = df_traj.reindex(columns=['time','x','y'])
                df_traj.reset_index(drop=True, inplace=True)
                df_traj.to_csv(Path(output_folder).joinpath(str(name)+str(traj)+".csv"),sep=',',index=False)