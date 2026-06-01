from mint.utils import get_file_list
import imageio.v2 as imageio
from pathlib import Path
from mint.denoising import wavelet
import trackpy as tp
from os.path import isfile
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

gfp_folder = r'/media/lumin/DATA/Coloc/20240906_Lamelle_1/GFP/'
nv_folder = r'/media/lumin/DATA/Coloc/20240906_Lamelle_1/NV/'

gfp_paths, gfp_names = get_file_list(Path(gfp_folder), '.tif')
nv_paths, nv_names = get_file_list(Path(nv_folder), '.tif')

gfp_root = gfp_names[0].split('#')[0]

for nv_path, nv_name in zip(nv_paths, nv_names):

    print(f'Opening {nv_name}')

    nv_frame = imageio.imread(nv_path)
    nv_frame = nv_frame.astype('float64')
    nv_frame = wavelet(nv_frame)

    nv_id = nv_name.split('#')[1].rstrip('.tif')

    nv_coords = tp.locate(nv_frame,
                          minmass=15000, #3000
                          diameter=13,
                          separation=10,
                          preprocess=False,
                          engine='numba')

    filecsv = Path(
        r'/media/lumin/DATA/Coloc/Lamelle1_results').joinpath(f'NV #{nv_id}.csv')
    nv_coords.to_csv(filecsv, sep='\t')

    gfp_file_name = gfp_root + f'#{nv_id}.tif'
    gfp_path = Path(gfp_folder).joinpath(gfp_file_name)
    print(f'Looking for {gfp_file_name}')

    if isfile(gfp_path):
        print('\tFound')
        gfp_frame = imageio.imread(gfp_path)
        gfp_frame = gfp_frame.astype('float64')
        gfp_frame = wavelet(gfp_frame)

        gfp_coords = tp.locate(gfp_frame,
                               minmass=2540, # 2540
                               diameter=13,
                               separation=10,
                               preprocess=False,
                               engine='numba')
        filecsv = Path(
        r'/media/lumin/DATA/Coloc/Lamelle1_results').joinpath(f'GFP #{nv_id}.csv')
        gfp_coords.to_csv(filecsv, sep='\t')

        joined = gfp_coords[['x', 'y']
                  ].merge(nv_coords[['x', 'y']], how='cross')
        joined['r_dist'] = np.sqrt((joined.x_x - joined.x_y)**2 +
                               (joined.y_x - joined.y_y)**2)

        filtered = joined.loc[joined.r_dist < 4]
        filtered = filtered.sort_values(by=['r_dist'])
        plot_coords = filtered[['x_y', 'y_y']]
        plot_coords = plot_coords.rename(columns={'x_y': 'x', 'y_y': 'y'})

        tp.annotate(nv_coords,
                    nv_frame,
                    plot_style={'markersize': 10},
                    color='red')
        tp.annotate(gfp_coords,
                    nv_frame,
                    plot_style={'markersize': 10},
                    color='green')
        tp.annotate(plot_coords,
                    nv_frame,
                    plot_style={'markersize': 10},
                    color='blue')
        plt.savefig(Path(
        r'/media/lumin/DATA/Coloc/Lamelle1_results').joinpath(f'MERGE #{nv_id}.png'),
            dpi=400)
        plt.close()

        filecsv = Path(
        r'/media/lumin/DATA/Coloc/Lamelle1_results').joinpath(f'MERGE #{nv_id}.csv')
        filtered.to_csv(filecsv, sep='\t')

        del gfp_frame
        del gfp_coords
        del joined
        del filtered
        del plot_coords
    del nv_coords
    del nv_frame
