import os
import pandas as pd
from pathlib import Path
import imageio.v3 as imageio

def longest_streak(lst):
    max_streak = 0
    current_streak = 0
    current_element = None

    for element in lst:
        if element == current_element:
            current_streak += 1
        else:
            current_streak = 1
            current_element = element

        max_streak = max(max_streak, current_streak)

    return max_streak

def get_file_list(input_folder, extension):

    path_list = []
    name_list = []

    if input_folder.endswith(f'{extension}'):
        path_list.append(Path(input_folder))
        name_list.append(Path(input_folder).name)

    else:
        for path, subfolder, files in os.walk(input_folder): # Scan entire folder structure for files
                for name in files:
                    if name.endswith(f'{extension}') == False:  # Check for correct file extension
                        continue # Skip to next file if not correct extension        

                    file_path = os.path.join(path, name) # Get file path of current file
                    path_list.append(file_path) # Append to file path list
                    name_list.append(name)
    
    return path_list, name_list


def GFP_mask2(path,name,trajectories):

    folder = Path(path).parent
    img = name[:-7]+'GFP.tif'

    try:
        img = imageio.imread(Path(folder).joinpath(img))
        # print('File found !')

    except FileNotFoundError:
        # print('File not found')
        print(Path(folder).joinpath(img))
        return trajectories

    threshold = 120
    mask = img
    mask[mask >= threshold] = 4096
    mask[mask < threshold] = 0
    
    mask_list = []

    for x,y in zip(trajectories.x, trajectories.y):
        if mask[int(y),int(x)]==4096:
            mask_list.append('in')
        else:
            mask_list.append('out')

    trajectories['GFP_mask'] = mask_list

    ratio_list = []
    streak_list = []

    df = pd.DataFrame()

    for traj in trajectories.rejoined_particle.unique():
        # print(traj)
        traj = trajectories[trajectories.rejoined_particle==traj].copy()
        ratio = len(traj[traj.GFP_mask=='in'])/len(traj)
        ratio = [ratio]*len(traj)
        ratio_list.extend(ratio)
        ratio = pd.DataFrame(ratio,columns=['ratio_df'])
        traj['ratio_df'] = ratio.values
        streak = longest_streak(list(traj.GFP_mask))
        streak = [streak]*len(traj)
        streak = pd.DataFrame(streak,columns=['streak_df'])
        traj['streak_df'] = streak.values
        streak_list.extend(streak)
        df = pd.concat((df,traj))

    trajectories['ratio'] = ratio_list
    # trajectories['streak'] = streak_list

    return df

def data_extraction(input_folder):

    path_list, name_list = get_file_list(str(input_folder), '_rejoined.csv')
    
    for (path, name) in zip(path_list,name_list): # Looping over file path list
        # print(name)
        # print(path)
        new_path = path.replace(str(input_folder),'/media/lumin/DATA')
        new_path = new_path.replace('_rejoined.csv','')
        new_path = new_path.replace((name.replace('_rejoined.csv','')),'',1)
        new_path = new_path.replace('//','/')
        new_path = new_path.replace('.nd2 ','.nd2')
        new_name = name.replace('_rejoined.csv','')
        # print(new_path)
        new_path = Path(new_path)
        processed_trajectory = pd.read_csv(path,sep='\t')
        try:
            processed_trajectory = processed_trajectory.drop(columns=['ratio','streak','GFP_mask'])
        except KeyError:
            processed_trajectory = processed_trajectory
        processed_trajectory = GFP_mask2(new_path,new_name,processed_trajectory)
        new_save = path.replace('/DATA_DEVRIM/','/DATA_DEVRIM_GFP/')
        processed_trajectory.to_csv(new_save,sep='\t')

data_extraction(r'/media/lumin/DATA/DATA_DEVRIM Results - 20231214_132940 gfp ok')