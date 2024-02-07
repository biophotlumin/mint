"""Functions used to output calculation results into files.
"""

import os
import sys
import scs
import yaml
import cvxpy

import numpy as np
import pandas as pd
import trackpy as tp
import matplotlib.pyplot as plt

from fpdf import FPDF
from io import StringIO
from pathlib import Path, PosixPath
from datetime import datetime

plt.switch_backend('agg')

def trajectory_output(output_file_path, name, process, output):
    """Writes trajectories into a .csv file.

    Uses pandas to_csv to write the DataFrame containing trajectories into a .csv file.

    :param output_file_path: Current output file path.
    :type output_file_path: str or Path
    :param name: Name of the file currently being processed.
    :type name: str
    :param process: Suffix appended to the file name.
    :type process: str
    :param output: DataFrame that is being dumped.
    :type output: DataFrame
    """

    # Create file path and name
    filecsv = Path(str(output_file_path)).joinpath(name+process+".csv")
    output.to_csv(filecsv, sep='\t')

def image_output(output_file_path, name, frames, trajectory, item):
    """Plots trajectories onto the first frame of a file. Inputs a dictionary,
        a string, a NumPy array, and a DataFrame.

    :param output_file_path: Current output file path.
    :type output_file_path: str or Path
    :param name: Name of the file currently being processed.
    :type name: str
    :param frames: Name of the file currently being processed.
    :type frames: str
    :param trajectory: DataFrame containing the trajectory to be plotted.
    :type trajectory: DataFrame
    :param item: Number of the trajectory to be plotted.
        Pass False to plot all trajectories.
    :type item: int, boolean
    """

    # Initialize plot
    plt.imshow(frames[0])
    # Plotting the y axis inverts it by default, so it must be inverted again
    plt.gca().invert_yaxis()
    bbox_props = dict(boxstyle="round", fc="w", ec="w")
    arrow_props = dict(arrowstyle="simple")

    trajectories = []
    if isinstance(item, bool):
        plt.title('Trajectories from '+name)
        filepng = Path(str(output_file_path)).joinpath(name+".png")
        for item in set(trajectory.rejoined_particle):
            trajectories.append(item)
    else:
        plt.title(f'Trajectory {item} from {name}')
        trajectories.append(item)
        filepng = Path(str(output_file_path)).joinpath(name+"-"+str(item)+".png")

    for item in trajectories: # Loop for each sub trajectory
        x = trajectory[trajectory.rejoined_particle == item].x
        y = trajectory[trajectory.rejoined_particle == item].y
        plt.plot(x, y, lw=0.5)
        if len(trajectories) > 1:
            x0, y0 = x.iloc[0], y.iloc[0]
            plt.annotate(str(item), xy=(x0, y0), xytext=(x0+30, y0+30), bbox=bbox_props,
                         arrowprops=arrow_props, size=4)

    plt.savefig(filepng, dpi=300)
    plt.close()

def trajectory_separation(output_file_path, name, trajectory, settings, sub):
    """Separates individual trajectories and writes them into .txt files.

    For each trajectory contained in the DataFrame, removes unnecessary columns,
        reindexes the remaining ones, and writes it into a .txt file.

    :param output_file_path: Current output file path.
    :type output_file_path: str or Path
    :param name: Name of the file currently being processed.
    :type name: string
    :param trajectory: Number of the trajectory being plotted.
    :type trajectory: int
    :param settings: Dictionary specifying wether or not
        a specific process should be applied.
    :type settings: dict
    :param sub: Subset of the DataFrame containing the trajectories.
    :type sub: DataFrame
    """

    file = Path(str(output_file_path)).joinpath(name+"-"+str(trajectory)+".txt")
    pd.options.mode.chained_assignment = None
    if settings['SNR_estimation'] is True: # Check wether SNR_estimation was used
        sub.drop(sub.columns.difference(['y', 'x', 'mass', 'feet']),
                 axis='columns', inplace=True) # If it was, drops the following columns
        sub.columns = ['Yraw', 'Xraw', 'mass', 'feet']
        sub = sub.reindex(columns=['Xraw', 'Yraw', 'mass', 'feet'])
        sub = sub.rename(columns={"mass": "Signal", "feet": "Noise"})
    else:
        # If not, drops only the columns that aren't added by SNR_estimation
        sub.drop(sub.columns.difference(['y', 'x', 'mass']),
                 axis='columns', inplace=True)
        sub.columns = ['Yraw', 'Xraw', 'mass']
        sub = sub.reindex(columns=['Xraw', 'Yraw', 'mass'])
        sub = sub.rename(columns={"mass": "Signal"})

    sub.to_csv(file, sep='\t', index=False)

def dict_dump_legacy(path, dict, file_name):
    """Writes the content of a dictionary into a .txt file

    :param path: Folder where the dictionary will be saved.
    :type path: str or Path
    :param dict: Dictionary to dump.
    :type dict: dict
    :param file_name: Name of the text file.
    :type file_name: string
    """

    with open(Path(path).joinpath(str(file_name)+".txt"), 'a') as dict_txt:
        for k, v in dict.items():
            dict_txt.write(f'{str(k)} : {str(v)}\n')

def dict_load_legacy(input_folder, dict):
    """Load a dictionary from a text file.

    :param input_folder: Folder where the dictionary is located.
    :type input_folder: str or Path
    :param dict: Name of text file to be loaded.
    :type dict: str
    :return: Loaded dictionary.
    :rtype: dict
    """

    loaded_dict = {}

    with open(input_folder.joinpath(f'{str(dict)}.txt')) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            k, v = line.split(' : ')
            loaded_dict[k] = v

    return loaded_dict

def dict_dump(path, data: dict, file_name: str, overwrite: bool=False):

    data = data.copy() # Just in case the dict has to be modified

    for k, v in data.items(): # YAML doesn't like PosixPath, convert to string
        if isinstance(v, PosixPath):
            data[k] = str(data[k])
        elif isinstance(v, np.float64):
            data[k] = float(data[k])

    file_path = Path(path).joinpath(f'{file_name}.yml')
    if file_path.is_file(): # Check for existing file
        if overwrite is True:
            with open(file_path, 'w') as f: # Overwrite if required
                yaml.dump(data, f)
        else:
            old_dict = yaml.safe_load(open(file_path))
            old_dict.update(data) # Otherwise update the preexisting dict
            with open(file_path, 'w') as f:
                yaml.dump(old_dict, f)
    else:
        with open(file_path, 'w') as f: # Otherwise create the file
            yaml.dump(data, f)

def dict_load(path, name: str):

    return yaml.safe_load(open(Path(path).joinpath(f'{name}.yml')))

# Functions generating PDF reports through FPDF
# Code is a bit messy and likely to break, but works for now
# (as long as there is no space in the name of the conditons)

def ind_page(pdf, list_var, input_folder, order, parameters, settings, n_cond):
    """Generate individual page.

    :param pdf: PDF to which the page is appended.
    :type pdf: FPDF instance
    :param list_var: Variable-specific results.
    :type list_var: List of strings
    :param input_folder: Folder where plots are located.
    :type input_folder: str or Path
    :param order: Order of experimental conditions.
    :type order: List of strings
    :param parameters: Dictionary containing parameters.
    :type parameters: dict
    :param settings: Dictionary containing settings.
    :type settings: dict
    :param n_cond: Number of experimental conditions.
    :type n_cond: int
    """

    pdf.set_font('Helvetica', 'B', size=13)
    pdf.cell(0, None, list_var[0].strip('\n').replace('_', ' '), align='C')
    pdf.set_font('Helvetica', size=10)
    pdf.ln(5)
    top_of_page = pdf.get_y()

    for path, subfolder, files in os.walk(input_folder): # Embed plots side by side
        for name in files:
            if name == 'Barplot '+str(order).replace(',', '')+'\
                  '+list_var[0].lower().strip('\n')+'.'+str(parameters['extension_out']).strip("'"):
                file_path = os.path.join(path, name)
                pdf.image(file_path, w=pdf.epw/2)
                if settings['clean_up'] is True:
                    os.remove(file_path)
            elif name == 'Boxplot '+str(order).replace(',', '')+'\
                  '+list_var[0].lower().strip('\n')+'.'+str(parameters['extension_out']).strip("'"):
                file_path = os.path.join(path, name)
                pdf.image(file_path, w=pdf.epw/2, x=(pdf.epw/2)+15, y=top_of_page)
                if settings['clean_up'] is True:
                    os.remove(file_path)
    pdf.ln(5)

    # Check wether or not a Dunn's test is included in the results
    if list_var[3] == '\n':
        dunn = 6
    else:
        dunn = 0

    line_height = pdf.font_size * 2.5
    pdf.cell(0, None, list_var[1].strip('\n'),
             new_y="LAST", new_x="LEFT")
    pdf.cell(0, None, list_var[2].strip('\n'),
             new_y="TOP", new_x="RIGHT", align='R')
    pdf.ln(5)
    pdf.cell(0, None, list_var[3+dunn].strip('\n'),
             new_y="LAST", new_x="LEFT")
    pdf.cell(0, None, list_var[4+dunn].strip('\n'),
             new_y="TOP", new_x="RIGHT", align='R')
    pdf.ln(5)
    pdf.cell(0, None, list_var[5+dunn].strip('\n'),
             new_y="LAST", new_x="LEFT")
    pdf.cell(0, None, list_var[6+dunn].strip('\n'),
             new_y="TOP", new_x="RIGHT", align='R')
    pdf.ln(10)

    top_of_ind_table = pdf.get_y()
    col_width = 20

    if dunn == 6: # Write Dunn's test table
        t_sring = str()
        for i in range(4, (4+n_cond+1), 1):
            t_sring = t_sring + (list_var[i])

        string_io = StringIO(t_sring)

        df = pd.read_csv(string_io, sep='\s+')

        data = []

        for i, v in enumerate(df.iloc[0].items()):
            data.append([x for x in df.iloc[i].items()])

        pdf.set_font('Helvetica', size=8)
        line_height = pdf.font_size * 2.5

        pdf.multi_cell(col_width, line_height, str("Dunn's test"), border=1, align='C',
                        new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size)

        for lb in df.columns:
            pdf.multi_cell(col_width, line_height, str(lb), border=1, align='C',
                        new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size)

        pdf.ln()
        for i, row in enumerate(data):
            pdf.multi_cell(col_width, line_height, str(df.index.values[i]),
                           border=1, align='C', new_x="RIGHT", new_y="TOP",
                              max_line_height=pdf.font_size)
            for datum in row:
                pdf.multi_cell(col_width, line_height, str(round(datum[1], 6)),
                               border=1, align='C', new_x="RIGHT", new_y="TOP",
                                  max_line_height=pdf.font_size)
            pdf.ln(line_height)
        # end_of_ind_table = pdf.get_y()
        dunn += 1

    string2 = str()
    for i in range(7+dunn, (7+dunn+n_cond+3), 1):
        string2 = string2 + (list_var[i])

    tt2 = StringIO(string2)

    df2 = pd.read_csv(tt2, sep='\s+')

    df2 = df2.drop([0])
    df2 = df2.drop(['deviation'], axis=1)
    df2 = df2.rename(columns={'Sample': 'Condition', 'size': 'Sample size',
                              'Standard': 'Standard deviation'})
    df2 = df2.round(decimals=6)

    pdf.set_xy(-110, top_of_ind_table)

    for lb in df2.columns: # Sample distribution parameters
        pdf.multi_cell(col_width, line_height, str(lb), border=1, align='C',
                    new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size)
    pdf.ln(line_height)
    pdf.set_xy(-110, pdf.get_y())


    for index, value in df2.iterrows():
        for v in value:
            pdf.multi_cell(col_width, line_height, str(v),
                           border=1, new_x="RIGHT", new_y="TOP",
                           max_line_height=pdf.font_size, align='C')
        pdf.ln(line_height)
        pdf.set_xy(-110, pdf.get_y())

    pdf.ln(5)

def generate_report(input_folder):
    """Generate an experiment report.

    :param input_folder: Folder containing plots and statistical test results.
    :type input_folder: str or Path
    """

    print('Generating report...\r', end='')

    input_folder = Path(input_folder)

    # Load dictionaries

    log = dict_load(input_folder, 'log')
    parameters = dict_load(input_folder, 'parameters')
    settings = dict_load(input_folder, 'settings')
    vars = dict_load(input_folder, 'vars')

    if settings['ordering'] is True:
        order = parameters['order']
    else:
        order = log['raw_order']

    # Only relevant calculation parameters and settings are kept in the report,
    # the complete dicts can still be found in parameters.txt and settings.txt

    disp_params = parameters.copy()
    disp_settings = settings.copy()

    for k in ['dpi', 'extension_out', 'order', 'extension_in']:
        disp_params.pop(k)

    if settings['stub_filtering'] is False:
        disp_params.pop('stub_filtering')

    if settings['MSD'] is False:
        disp_params.pop('msd')

    if settings['rejoining'] is False:
        for k in ['threshold_t', 'threshold_r']:
            disp_params.pop(k)

    if settings['SNR_estimation'] is False:
        disp_params.pop('base_level')

    if settings['polynomial_fit'] is False:
        for k in ['len_cutoff', 'threshold_poly3']:
            disp_params.pop(k)

    if settings['minimization'] is False:
        disp_params.pop('sigma')

    for k in ['ordering', 'clean_up']:
        disp_settings.pop(k)

    # order = ast.literal_eval(str(order))

    n_cond = len(order)

    start_time_date = str(log['identifier'])[-15:]
    start_date = f'{start_time_date[0:4]}/{start_time_date[4:6]}/{start_time_date[6:8]}'
    start_time = f'{start_time_date[9:11]}h{start_time_date[11:13]}'

    # Main page

    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('Helvetica', 'B', size=18)
    pdf.cell(0, 10, "M.I.N.T analysis report", align='C')
    pdf.ln()
    pdf.set_font('Helvetica', size=12)
    pdf.cell(0, 10, (f"Python version : {sys.version} | "
                     "Conda environment : {sys.executable.split('/')[-3]}"), align='C')
    pdf.ln()
    pdf.cell(0, 10, (f"Trackpy version {tp.__version__}"
                     f" | NumPy version {np.__version__}"
                     f" | CVXPY version {cvxpy.__version__}"
                     f" | SCS version {scs.__version__}"), align='C')
    pdf.ln()
    pdf.cell(0, 10, f"Input folder :  {log['input_folder']}", align='C')
    pdf.ln()
    pdf.cell(0, 10, (f"Started : {start_date} at {start_time}"
                     f" | Report generated : "
                     f"{str(datetime.now().strftime('%Y/%m/%d at %Hh%M'))}"
                     f" | Duration : {log['duration']}"), align='C')
    pdf.ln(15)
    pdf.set_font('Helvetica', 'B', size=12)
    pdf.cell(45, 10, 'Parameters', align='C')
    pdf.cell(60, 10, 'Settings', align='C')
    pdf.ln(10)
    pdf.set_font('Helvetica', size=8)
    line_height = pdf.font_size * 2.5
    top_of_table = pdf.get_y()

    # Building tables

    for k, v in disp_params.items():
        pdf.multi_cell(30, line_height, f'{k}', border=1,
                new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, align='C')
        pdf.multi_cell(15, line_height, f'{v}', border=1,
                new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, align='C')
        pdf.ln(line_height)

    pdf.set_xy(-150, top_of_table)
    for k, v in disp_settings.items():
        pdf.multi_cell(35, line_height, f'{k}', border=1,
                new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, align='C')
        pdf.multi_cell(15, line_height, 'Yes' if v else 'No', border=1,
                new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, align='C')
        pdf.ln(line_height)
        pdf.set_xy(-150, pdf.get_y())

    # Additional logging

    pdf.ln()
    pdf.set_font('Helvetica', size=10)
    pdf.set_xy(112, top_of_table)
    pdf.cell(None, None, f"Number of files analyzed : {log['n_files']}",
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, f"Number of files retained in tracking: {log['n_a_file']}",
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, f"Number of files retained in analysis: {log['n_t_file']}",
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, f"Number of trajectories found : {log['n_traj']}",
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, f"Number of trajectories rejoined : {str(log['n_rejoined'])}",
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, (f"Number of trajectories rejected by polynomial fit : "
                          f"{log['r_poly']}"),
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, (f"Number of trajectories rejected by abnormal speed : "
                          f"{log['r_speed']}"),
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, (f"Number of trajectories retained for analysis : "
                          f"{log['n_a_traj']}"),
             align='L', new_y="NEXT", new_x="LEFT")
    pdf.cell(None, None, f"Number of conditions : {n_cond}",
             align='L', new_y="NEXT", new_x="LEFT")

    for i in order:
        sub_x = pdf.get_x()
        pdf.ln(5)
        pdf.set_xy(sub_x, pdf.get_y())
        pdf.cell(None, None, (f"Number of animals analyzed for condition {i} : "
                              f"{log['n_a_animal_'+i]}"), align='L',
                              new_y="NEXT", new_x="LEFT")
        pdf.cell(None, None, (f"Number of files analyzed for condition {i} : "
                              f"{log['n_a_file_'+i]}"), align='L',
                              new_y="NEXT", new_x="LEFT")
        pdf.cell(None, None, (f"Number of trajectories analyzed for condition {i} : "
                              f"{log['n_a_traj_'+i]}"), align='L',
                              new_y="NEXT", new_x="LEFT")

    pdf.add_page()

    pdf.set_font('Helvetica', size=10)
    page_break = 0
    list_var = []

    stats_path = Path(input_folder).joinpath('Statistical test results.txt')

    # Looping over statistical variables for subsequent pages

    with open(stats_path) as f:
        lines = f.readlines()
        for line in lines:
            lt = (line.lower().strip('\n'))
            if lt in vars:
                if list_var:
                    if page_break >= 2: # Maximum of two variables per page
                        pdf.add_page()
                        page_break = 0
                    ind_page(pdf, list_var, input_folder, order,
                             parameters, settings, n_cond)
                    page_break += 1
                list_var = []
                list_var.append(line.strip(''))
            else:
                list_var.append(line.strip(''))
        ind_page(pdf, list_var, input_folder, order, parameters, settings, n_cond)

    pdf_path = Path(input_folder).joinpath(
        str(datetime.now().strftime('%Y%m%d_%H%M%S'))+'_Experiment_report.pdf')
    pdf.output(pdf_path)
    print('Generating report... Done\r')

if __name__ == '__main__':
    generate_report(r'')
