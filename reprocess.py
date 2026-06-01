from mint.utils import folder_structure_creation, csv_sniffer, print_pb
from mint.traj_calc import MSD_filtering, rejoining
from mint.output import trajectory_output, image_output
from pathlib import Path
import os
import pandas as pd
import warnings
import trackpy as tp
import numpy as np
import time
import yaml

root_input_folder = r'/media/lumin/DATA/TIRF/'
input_folder = r'/media/lumin/DATA/TIRF/RVG29C Results - 20251105_113343 reprocess MSD'

settings = {
  'MSD': True,
  'rejoining': True,
  'SNR_estimation': False,
  'individual_images': False,
  'individual_txt': False,
  'group_image': True,
  }

parameters = {
  'threshold_t': 10,
  'threshold_r': 40,
  'msd': 2
  }

path_list = []
frames = np.zeros((1, 1002, 1004))

output_folder, identifier = folder_structure_creation(input_folder)[0:2]

if output_folder.exists() is False:
    os.makedirs(output_folder)

for path, subfolder, files in os.walk(input_folder):
    for name in files:
        if name.endswith('nd2.csv') is False:
            continue

        # Build output file path
        file_path = os.path.join(path, name)
        path_list.append(file_path)

j_max = len(path_list)

start = time.time()

already = yaml.safe_load(
        open(Path(r'/media/lumin/DATA/TIRF/RVG29C/already_reprocessed.yml')))


for (path, j) in zip(path_list, [j for j in range(len(path_list))]): # enumerate
    name = Path(path).name.rstrip("_rejoined.csv")

    raw_trajectory = pd.read_csv(path, sep=csv_sniffer(path))

    output_subfolder = str(path).replace(str(root_input_folder), '')
    out_path = Path(output_folder).joinpath(output_subfolder)

    print_pb(f'\nProcessing {name}', j, j_max)
    os.makedirs(out_path)

    print_pb('\tMSD filtering', j, j_max)
    processed_trajectory = MSD_filtering(raw_trajectory, parameters['msd'])
    if len(processed_trajectory) == 0: # Check if any trajectories were found.
        # If not, the threshold might be too high.
        warnings.warn('No trajectories retained, '
                        'MSD threshold might be too high')

        n_particles = len(raw_trajectory[raw_trajectory.frame == 0])
        static = tp.filter_stubs(raw_trajectory, int(2310*0.9))
        n_static = static.particle.nunique()

        pd.DataFrame.from_dict(
            {
                'n_particles': [n_particles],
                'n_static': [n_static],
            }
        ).to_csv(Path(out_path).joinpath('static.csv'))
        continue

    if settings['rejoining']:
        print_pb('\tRejoining', j, j_max)
        processed_trajectory, n_rejoined = rejoining(processed_trajectory,
                                                        parameters['threshold_t'],
                                                        parameters['threshold_r'])

    else:
        processed_trajectory['rejoined_particle'] = processed_trajectory['particle']

    # Estimating ratio of moving particles
    first_frame = raw_trajectory[raw_trajectory.frame == 0]
    n_particles = len(first_frame)
    n_particles = [n_particles]*len(raw_trajectory)
    n_particles = pd.DataFrame(n_particles, columns=['n_particles'])
    processed_trajectory = processed_trajectory.reset_index(drop=True)
    processed_trajectory = pd.concat([processed_trajectory, n_particles],
                                        axis=1, join='inner')

    # Number of static particles
    static = MSD_filtering(raw_trajectory,
                            parameters['msd'],
                            highpass=False)
    if len(static) == 0:
        n_static = 0
    else:
        static = tp.filter_stubs(static, int(2310*0.9))
        n_static = static.particle.nunique()

    n_static = [n_static]*len(raw_trajectory)
    n_static = pd.DataFrame(n_static, columns=['n_static'])
    processed_trajectory = processed_trajectory.reset_index(drop=True)
    processed_trajectory = pd.concat([processed_trajectory, n_static],
                                        axis=1, join='inner')

    # Dumping rejoined trajectories into csv file
    trajectory_output(out_path, name, "_rejoined", processed_trajectory)

    # Per trajectory data extraction
    if (settings['individual_images'] or
        settings['individual_txt'] or
        settings['group_image']):
        print_pb('\tSaving plots and trajectories', j, j_max)
        for item in set(processed_trajectory.particle):
            sub_trajectory = processed_trajectory[processed_trajectory
                                                    .particle == item]

        # Plot all trajectories onto the first frame of the video
        if settings['group_image']:
            image_output(out_path,
                            name,
                            frames,
                            processed_trajectory,
                            False)
end = time.time()
duration = end - start
f_duration = f'{int(duration//3600)}h{int((duration%3600)/60):02d}'
print(f'Total runtime : {f_duration}')