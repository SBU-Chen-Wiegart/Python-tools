"""
File: CMS GIWAXS AND GISAXS
Name: Cheng-Chu Chung
TODO: Auto import ref peaks from PDF card data, then put data into dictionary
----------------------------------------
Color palettes for Python: https://jiffyclub.github.io/palettable/#palette-interface
"""

from pathlib import Path
from collections import defaultdict
import pandas as pd
from pprint import pprint
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy import signal
import configparser
import palettable as pltt
import peakutils
import matplotlib
import sys
from colorama import Fore, Back, Style, init

from torch.backends.quantized import engine

print("Before, Backend used by matplotlib is: ", matplotlib.get_backend())
matplotlib.rcParams['backend'] = 'wxAgg'
print("After, Backend used by matplotlib is: ", matplotlib.get_backend())

# Step 1: Give your data directory
# GISAXS
# INPUT_PATH = r"D:\Research data\SSID\202411\20241104 CMS AE\saxs\analysis\b58-04-ScVMnSc-AE_afterAE"
# CONFIG_FILE = r"D:\Research data\SSID\202411\20241104 CMS AE\saxs\analysis\b58-04-ScVMnSc-AE_afterAE\Plot\b58-04-ScVMnSc-AE_afterAE-SAXS.ini"

# GIWAXS
INPUT_PATH = r"D:\Research data\SSID\202506\20250623 CMS CrCuNiCr AE\maxs\analysis\b57-01_x29.601_5.00s_204829_dioptas"
CONFIG_FILE = r"D:\Research data\SSID\202506\20250623 CMS CrCuNiCr AE\maxs\analysis\b57-01_x29.601_5.00s_204829_dioptas\Plot\b57-01-CrCuNiCr.ini"

OUTPUT_FOR_JADE = True
IF_SAVE = False

# OUTPUT_PATH = Path(f'{INPUT_PATH}\Output_files')
# Save with config file
OUTPUT_PATH = Path(f'{Path(CONFIG_FILE).parent}\Output_files_{Path(CONFIG_FILE).stem}')
# If you are the Mac user, please change the above line to:
# OUTPUT_PATH = Path(CONFIG_FILE).parent / f'Output_files_{Path(CONFIG_FILE).stem}'

# Step 2: Confirm your config file
CONFIG = configparser.ConfigParser()

if Path(CONFIG_FILE).is_file():
    CONFIG.read(CONFIG_FILE, encoding="utf8")
    is_ini = True
    print("=================")
    print("Use .ini input")
    print("-----------------")
else:
    is_ini = False
    print("=================")
    print("Manually input so please remove all eval commands if error occurs or file is not found")
    print("-----------------")

# FILE_TYPE = '.dat'                                                                                 # Archive - Check your file type
SCANID_INDEX = -2                                                                                    # Use scan ID to sort the files
PATTERN = eval(CONFIG['samples']['pattern'])
FILENAME_KEYWORD = '_Tc'                                                                          # Find the last index of the keyword in the filename
FILENAME_KEYWORD_OFFSET = 6                                                                          # "th"0.250
LEGEND_HEAD_KEYWORD = 'x'
LEGEND_TAIL_KEYWORD = '_th'
ANGLE_RANGE = eval(CONFIG['samples']['angle_range'])
SAMPLE_LIST = eval(CONFIG['samples']['sample_list'])
SAXS_COLUMN_NAME = ['#', 'qr', 'I']                                                                  # May be updated
WAXS_COLUMN_NAME = ['#', 'q', 'qerr', 'I(q)', 'I(q)err']                                             # 2024-3 update to data_dict['I_list'][index] = dataframe[:, 2]
# WAXS_COLUMN_NAME = ['#', 'q', 'I(q)err', 'q']                                                      # May be updated, then update: data_dict['I_list'][index] = dataframe[:, 1]
# WAXS_COLUMN_NAME = ['#', 'q', 'qerr', 'I(q)']                                                      # Before 2023
DIOPTAS_COLUMN_NAME = ['#', 'q_A^-1', 'I']
FIGURE_SIZE = eval(CONFIG['format']['figure_size'])
BATCH_NUMBER, COMPOSITION, CONDITION, INCIDENT_ANGLE = eval(CONFIG['legends']['sample_condition'])   # Whether you want to show them in the legend
PALETTE = eval(CONFIG['format']['palette'])                                                          # pld.Spectral_4_r  # _r if you want to reverse the color sequence
CMAP = PALETTE.mpl_colormap                                                                          # .mpl_colormap attribute is a continuous, interpolated map
COLOR_INCREMENT = eval(CONFIG['format']['color_increment'])                                          # Adjust your color gradient, 0 for default
OFFSET = eval(CONFIG['format']['offset'])                                                            # Value you want to add to an y offset for each curve.
SAMPLE_LABEL = eval(CONFIG['legends']['sample_label'])
GISAXS_MODE = eval(CONFIG['data_processing']['gisaxs_mode'])
FIRST_DATAPOINT = eval(CONFIG['data_processing']['first_datapoint'])
XRANGE = eval(CONFIG['format']['xrange'])
YRANGE = eval(CONFIG['format']['yrange'])
LEGEND_LOCATION = eval(CONFIG['format']['legend_location'])
TITLE = eval(CONFIG['format']['output_filename'])
# OUTPUT_FOR_JADE = eval(CONFIG['format']['output_for_jade'])
# IF_SAVE = eval(CONFIG['format']['if_save'])
SUB_DEGREE = eval(CONFIG['data_processing']['bgsub_degree'])


def main():
    files = Path(INPUT_PATH).glob(f'*{PATTERN}')  # Call Path again to grab the file

    if len(list(files)) == 0:
        print("No files found, please check your input path or pattern")
        return
    else:
        files = Path(INPUT_PATH).glob(f'*{PATTERN}')    # Call Path again to grab the file

        # Sort the filename by scan ID:
        files = sorted(files, key=lambda x: int(x.name.split('_')[SCANID_INDEX]))

    if OUTPUT_FOR_JADE:
        OUTPUT_PATH.mkdir(exist_ok=True)  # Creates the folder only if it doesn't exist

    if ANGLE_RANGE == 'wide':
        print('=============================================')
        print('Grazing-incidence wide-angle X-ray Scattering')
        print('---------------------------------------------')
        giwaxs(files)
    elif ANGLE_RANGE == 'small':
        print('=============================================')
        print('Grazing-incidence small-angle X-ray Scattering')
        print('---------------------------------------------')
        gisaxs(files)
    else:
        print("Please enter the ANGLE_RANGE ('small' or 'wide')")


def giwaxs(files):
    q_and_I_list = sorted_data(files, mode=ANGLE_RANGE)
    if SUB_DEGREE != 0:
        background_subtraction(q_and_I_list, degree=SUB_DEGREE)   # Add a new list with background subtraction
        giwaxs_plot(q_and_I_list, mode='bg_sub')
    else:
        giwaxs_plot(q_and_I_list)



def gisaxs(files):
    q_and_I_list = sorted_data(files, mode=ANGLE_RANGE)
    gisaxs_plot(q_and_I_list, mode=GISAXS_MODE)


def sorted_data(files, mode=ANGLE_RANGE):
    data_dict = {'q_list': {}, 'I_list': {},'filename_list': {}, 'background_subtraction_list': {}}   # Create a dictionary to store all information
    for index, file in enumerate(files):
        scattering_data = file.resolve()  # Make the path absolute, resolving any symlinks
        data_dict['filename_list'][index] = scattering_data.name
        print(f'{index:>2} {scattering_data.name}')
        if mode == 'wide':
            # dataframe = pd.read_table(scattering_data, sep="\s+",
            #                           usecols=WAXS_COLUMN_NAME).to_numpy() \
            #     if scattering_data.suffix == ".dat" \
            #     else pd.read_table(scattering_data, sep="\s+", usecols=DIOPTAS_COLUMN_NAME, skiprows=22).to_numpy() \
            #     # '#' is q column and 'qerr' is I(q) column
            # data_dict['q_list'][index] = dataframe[:, 0]
            # data_dict['I_list'][index] = dataframe[:, 2]          # WAXS column may be different

            # Determine which columns to use and how many rows to skip based on file extension
            if scattering_data.suffix == ".dat":
                columns_to_use = WAXS_COLUMN_NAME                   # WAXS column may be different
                skip_rows = 0
                intensity_col_index = 2                             # I(q) is in the third column
            else:
                columns_to_use = DIOPTAS_COLUMN_NAME
                skip_rows = 22
                intensity_col_index = 1                             # I(q) is in the second column

            # Read the file into a DataFrame
            df = pd.read_table(
                scattering_data,
                sep=r"\s+",
                # usecols=columns_to_use,
                # skiprows=skip_rows,
                comment="#",                # Skip comment lines starting with #
                header=None,                # Don't use the first row as column headers
                names=columns_to_use[1:],   # Use the second column name as header
                index_col=False             # Ensure the first column is not used as the index
            )

            # Convert DataFrame to NumPy array
            data_array = df.to_numpy()

            # Extract q and I(q) values
            q_values = data_array[:, 0]
            intensity_values = data_array[:, intensity_col_index]

            # Store the results in the dictionary
            data_dict['q_list'][index] = q_values
            data_dict['I_list'][index] = intensity_values

            if OUTPUT_FOR_JADE:
                    out_file(data_dict['q_list'][index], data_dict['I_list'][index], f'C_{scattering_data.name}')

        elif mode == 'small':
            dataframe = pd.read_table(scattering_data, sep="\s+",
                                      usecols=SAXS_COLUMN_NAME) \
                .to_numpy()  # '#' is q column and 'qerr' is I(q) column

            # q^2 may include the negative q range so we do a filter
            data_dict['q_list'][index] = dataframe[FIRST_DATAPOINT:, 0]
            # Find the index to have the same shape for y dataframe
            data_dict['I_list'][index] = dataframe[FIRST_DATAPOINT:, 1]

        # if index == 0:
        #     pprint(data_dict)

    return data_dict


def background_subtraction(q_and_I_list, degree=3):
    """
    Subtract background from data and save data w/ background subtracted
    :return: None
    """
    # list_dict = sorted_data()
    print('\n'+'='*50)
    print(f'Background subtraction was applied with degree = {degree}')
    print('-'*50)
    for index in q_and_I_list['q_list']:
        x = q_and_I_list['q_list'][index]
        y = q_and_I_list['I_list'][index]
        filename = q_and_I_list['filename_list'][index]
        # print(f'{index:>2} {filename}')
        # pos = filename[filename.find('pos'):filename.find('pos') + 4]
        base_line = peakutils.baseline(y, degree)  # generate baseline
        y_corrected = y - base_line  # subtract baseline from y-values
        q_and_I_list['background_subtraction_list'][index] = y_corrected
        if OUTPUT_FOR_JADE:
                out_file(x, y_corrected, f'C_bgsub_{filename}')


def giwaxs_plot(q_and_I_list, mode='raw'):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    title = TITLE
    y_max_lim = 0
    y_min_lim = 10000
    y = 0   # Initial value
    increment = 0
    print('\n=============================================')
    print('Plot:')
    print('---------------------------------------------')
    for index in SAMPLE_LIST:
        try:
            x = q_and_I_list['q_list'][index]
        except KeyError:
            print(f"Index {index} not found in the data folder. Please check {Style.BRIGHT}{Fore.RED}SAMPLE_LIST{Style.RESET_ALL}...")
            return
        if mode == 'raw':
            y = q_and_I_list['I_list'][index]
        elif mode == 'bg_sub':
            y = q_and_I_list['background_subtraction_list'][index]
        filename = q_and_I_list['filename_list'][index]
        print(f'{index:>2} {filename}')

        # Plot label
        curve_info = filename[filename.find('_') + 1:filename.find('_', 3)]
        batch_number = filename[filename.find('b'):filename.find('b') + 6] if BATCH_NUMBER else '..'
        composition = filename[filename.find('b') + 7:filename.find('-', 10)] if COMPOSITION else '..'
        condition = filename[filename.find('-', 10) + 1:filename.find('_pos')] if CONDITION else '..'
        incident_angle = f"{filename[filename.find('th') + 2:filename.find('th') + 6]} degree" if INCIDENT_ANGLE else '..'
        if len(SAMPLE_LABEL) == 0:
            plot_label = f"{batch_number}/{composition}/{condition}/{incident_angle}"
            plot_label = filename[filename.find(LEGEND_HEAD_KEYWORD):filename.find(LEGEND_TAIL_KEYWORD)] \
                if LEGEND_HEAD_KEYWORD in filename else plot_label
        else:
            plot_label = SAMPLE_LABEL[SAMPLE_LIST.index(index)] if len(SAMPLE_LABEL) == len(SAMPLE_LIST) else filename

        # Title
        if TITLE != 'Auto' and TITLE != '':
            title = TITLE
        elif TITLE == '':
            title = Path(CONFIG_FILE).stem
        elif incident_angle == '..':
            title = f"{filename[filename.find('th') + 2:filename.find('th') + 6]} degree"
        else:
            title = curve_info
        
        if mode == 'bg_sub':
            title += ' bg_sub'

        # Do the plot
        color_idx = np.linspace(0, 1, len(SAMPLE_LIST) + COLOR_INCREMENT)
        plt.plot(x, y + OFFSET * increment, linewidth=2, color=CMAP(color_idx[increment+COLOR_INCREMENT]), label=plot_label)  # full info
        increment += 1
        if y.max() > y_max_lim:
            y_max_lim = y.max()
        if y.min() < y_min_lim:
            y_min_lim = y.min()

    print(f'Total number of curves: {increment}')

    # Frame linewidth
    spineline = ['left', 'right', 'top', 'bottom']
    for direction in spineline:
        ax.spines[direction].set_linewidth(3)

    # Plotting format
    plt.xlabel('q ($\mathregular{\AA}^{-1}$)', fontsize=18)
    plt.ylabel('I(q)', fontsize=18, labelpad=10)
    plt.xticks(fontsize=14)
    # plt.yticks([])  # Disable ticks
    ax.tick_params(width=3)
    if len(XRANGE) != 0:
        plt.xlim(XRANGE[0], XRANGE[1])
    if len(YRANGE) != 0:
        plt.ylim(YRANGE[0], YRANGE[1])
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc=LEGEND_LOCATION, framealpha=1, frameon=True, fontsize=12, reverse=True)
    plt.title(title, fontsize=18, pad=15)
    # ax.set_aspect('auto')    # Aspect ratio
    plt.tight_layout()
    if IF_SAVE:
        plt.yticks([])  # Disable ticks
        plt.tight_layout()
        config_file_location = Path(CONFIG_FILE).parent
        # If you are the Mac user, please change the above line to:
        # config_file_location = Path(CONFIG_FILE).parent
        output_filename = check_filename_repetition(title, config_file_location)
        plt.savefig("{}/{}.png".format(config_file_location, output_filename), dpi=300, transparent=False)
    plt.show()


def gisaxs_plot(q_and_I_list, mode='Intensity', xrange=(0.004, 0.1), yrange=(0, 120)):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    title = TITLE
    y_max_lim = 0
    y_min_lim = 10000
    y = 0  # Initial value
    increment = 0
    print('\n=============================================')
    print('Plot:')
    print('---------------------------------------------')

    for index in SAMPLE_LIST:
        x = q_and_I_list['q_list'][index]
        y = q_and_I_list['I_list'][index]
        filename = q_and_I_list['filename_list'][index]
        print(f'{index:>2} {filename}')

        # Plot label
        curve_info = filename[filename.find('_') + 1:filename.find('_', 3)]
        batch_number = filename[filename.find('b'):filename.find('b') + 6] if BATCH_NUMBER else '..'
        composition = filename[filename.find('b') + 7:filename.find('-', 10)] if COMPOSITION else '..'
        condition = filename[filename.find('-', 10) + 1:filename.find('_pos')] if CONDITION else '..'
        incident_angle = f"{filename[filename.find('th') + 2:filename.find('th') + 6]} degree" if INCIDENT_ANGLE else '..'
        if len(SAMPLE_LABEL) == 0:
            plot_label = f"{batch_number}/{composition}/{condition}/{incident_angle}"
            plot_label = filename[filename.find(LEGEND_HEAD_KEYWORD):filename.find(LEGEND_TAIL_KEYWORD)] \
                if LEGEND_HEAD_KEYWORD in filename else plot_label
        else:
            plot_label = SAMPLE_LABEL[SAMPLE_LIST.index(index)]

        # Title
        if TITLE != 'Auto' and TITLE != '':
            title = TITLE
        elif TITLE == '':
            title = Path(CONFIG_FILE).stem
        elif incident_angle == '..':
            title = f"{filename[filename.find('th') + 2:filename.find('th') + 6]} degree"
        else:
            title = curve_info
        title += ' ' + mode

        # Do the plot

        color_idx = np.linspace(0, 1, len(SAMPLE_LIST) + COLOR_INCREMENT)
        if mode == 'Guinier':
            plt.plot(x**2, np.log(y) + OFFSET * increment, 'o', linewidth=3, color=CMAP(color_idx[increment + COLOR_INCREMENT]),
                     label=plot_label)  # full info
            guinier_region_fit(x, y, color=CMAP(color_idx[increment + COLOR_INCREMENT]))
        
        elif mode == 'Guinier Peak':
            plt.plot(x**2, np.log(x*y) + OFFSET * increment, linewidth=3, color=CMAP(color_idx[increment + COLOR_INCREMENT]),
                     label=plot_label)  # full info
            new_x_for_our_of_range = x[x < np.sqrt(XRANGE[1])]
            new_y_for_our_of_range = y[:len(new_x_for_our_of_range)]
            q_max = x[np.where(x * y == np.max(new_x_for_our_of_range * new_y_for_our_of_range))][0]
            q_max = np.round(q_max, 4)
            print(f'q_max = {q_max}')
            Rg = np.round(1.5**0.5/q_max, 2)
            print(f'Rg = {Rg}')
            print(f'Check guinier region q x Rg: {q_max*Rg}')
            text = f' Rg = {Rg}'
            ax.text(q_max**2, np.log(np.max(x*y)), text, fontsize=14, color='b')
        
        elif mode == 'Peak Enhanced':
            plt.plot(x, x*y + OFFSET * increment, linewidth=3, color=CMAP(color_idx[increment + COLOR_INCREMENT]),
                     label=plot_label)  # full info
            # fitting_start = 670-FIRST_DATAPOINT
            # guassian_peak_fit(x[fitting_start:fitting_start+50], x[fitting_start:fitting_start+50]*y[fitting_start:fitting_start+50])
        
        else:
            plt.plot(x, y + OFFSET * increment, linewidth=3, color=CMAP(color_idx[increment + COLOR_INCREMENT]),
                 label=plot_label)  # full info

        increment += 1
        if y.max() > y_max_lim:
            y_max_lim = y.max()
        if y.min() < y_min_lim:
            y_min_lim = y.min()

    print(f'Total number of curves: {increment}')

    # Frame linewidth
    spineline = ['left', 'right', 'top', 'bottom']
    for direction in spineline:
        ax.spines[direction].set_linewidth(3)

    if mode == 'Guinier':
        plt.xlabel('$\mathregular{q_r^2 (\AA^{-1})}$', fontsize=20)
        plt.ylabel('Log[I(q)]', fontsize=20, labelpad=10)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    elif mode == 'Guinier Peak':
        plt.xlabel('$\mathregular{q_r^2 (\AA^{-1})}$', fontsize=20)
        plt.ylabel('Log[qI(q)]', fontsize=20, labelpad=10)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        plt.xlim(xrange)
        plt.ylim(yrange)
    elif mode == 'Paper':
        plt.xlabel('$\mathregular{q_r}$ ($\mathregular{\AA}^{-1}$)', fontsize=20)
        plt.ylabel('I(q)', fontsize=20, labelpad=10)
        plt.xscale('log')
        plt.yscale('log')
    elif mode =='Peak Enhanced':
        plt.xlabel('$\mathregular{q_r}$ ($\mathregular{\AA}^{-1}$)', fontsize=22)
        plt.ylabel('qI(q)', fontsize=22, labelpad=10)
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.xlabel('$\mathregular{q_r}$ ($\mathregular{\AA}^{-1}$)', fontsize=20)
        plt.ylabel('I(q)', fontsize=20, labelpad=10)
        plt.yscale('log')

    # Plotting format
    if len(XRANGE) != 0:
        # xrange == XRANGE
        plt.xlim(XRANGE)
    if len(YRANGE) != 0:
        # yrange == YRANGE
        plt.ylim(YRANGE)

    plt.xticks(fontsize=20)
    # ax.get_xaxis().set_tick_params(which='minor', width=3)
    plt.yticks(fontsize=20)
    # plt.yticks([])  # Disable ticks
    ax.tick_params(which='major', length=8, width=3)
    ax.tick_params(which='minor', length=5, width=3)
    plt.legend(loc=LEGEND_LOCATION, framealpha=1, frameon=False, fontsize=18, reverse=True)
    plt.title(title, fontsize=18, pad=15)
    plt.tight_layout()
    if IF_SAVE:
        config_file_location = Path(CONFIG_FILE).parent
        # If you are the Mac user, please change the above line to:
        # config_file_location = Path(CONFIG_FILE).parent
        output_filename = check_filename_repetition(title, config_file_location)
        plt.savefig("{}/{}.png".format(config_file_location, output_filename), dpi=300, transparent=False)
    plt.show()


def guinier_region_fit(x, y, color='k'):
    new_x_for_our_of_range = x[x < np.sqrt(XRANGE[1])]
    new_y_for_our_of_range = y[:len(new_x_for_our_of_range)]
    q_max_from_GPA = x[np.where(x * y == np.max(new_x_for_our_of_range * new_y_for_our_of_range))][0]
    print(f'q_max_from_GPA = {np.round(q_max_from_GPA, 4)}')
    Rg = np.round(1.5 ** 0.5 / q_max_from_GPA, 2)
    print(f'Rg = {Rg:>> 10}')
    print(f'Check guinier region qRg: {q_max_from_GPA * Rg}')

    q_max_index = np.where(new_x_for_our_of_range == q_max_from_GPA)[0][0]
    # slope, intercept, rvalue, pvalue, stderr = stats.linregress(
    #     new_x_for_our_of_range[:q_max_index] ** 2, np.log(new_y_for_our_of_range[:q_max_index]))
    q_min_index, q_max_index, slope, intercept, rvalue, pvalue, stderr = guinier_region_check(new_x_for_our_of_range,
                                                                                 new_y_for_our_of_range,
                                                                                 q_max_index, q_max_from_GPA)
    print(f'Slope: {slope}')
    print(f'Intercept: {intercept}')
    # print(f'Pearson correlation coefficient: {rvalue}')
    # print(f'p-value: {pvalue}')
    # print(f'Standard error: {stderr}')
    # print(f'Slope: {intercept_stderr}')
    RG = np.sqrt(slope*-3)
    print(f'RG = {RG:>> 50}')
    print(f'Check qRG maximum: {x[q_max_index] * RG}')
    q_min = new_x_for_our_of_range[q_min_index]
    print(f'Check qRG minimum: {q_min * RG}')
    # plt.plot(x**2, np.log(y), 'o', label='original data')
    plt.plot(x[q_min_index:q_max_index]**2, intercept + slope * x[q_min_index:q_max_index]**2, '--k')   # , label='fitted line'
    plt.plot(x[q_min_index:q_max_index] ** 2, intercept + slope * x[q_min_index:q_max_index] ** 2 - np.log(y[q_min_index:q_max_index]), '--', color=color) #, label='residual line'
    print('---------------------------------------------')


def guinier_region_check(x, y, q_max_index, q_max_from_GPA):
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(
        x[:q_max_index] ** 2, np.log(y[:q_max_index]))
    RG = np.sqrt(slope * -3)
    print(f'RG = {RG:>> 30}')
    print(f'Check qRG maximum from GPA: {q_max_from_GPA * RG}')
    q_max = q_max_from_GPA
    q_min_index = 0
    # q max approach-
    while q_max * RG > 1.3:
        q_max_index -= 1
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(
            x[q_min_index:q_max_index] ** 2, np.log(y[q_min_index:q_max_index]))
        RG = np.sqrt(slope * -3)
        q_max = x[q_max_index]
    # q max approach+
    while q_max * RG < 1.3:
        q_max_index += 1
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(
            x[q_min_index:q_max_index] ** 2, np.log(y[q_min_index:q_max_index]))
        RG = np.sqrt(slope * -3)
        q_max = x[q_max_index]
    # q min approach+
    while x[q_min_index] * RG < 0.9:
        q_min_index += 1
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(
            x[q_min_index:q_max_index] ** 2, np.log(y[q_min_index:q_max_index]))
        RG = np.sqrt(slope * -3)

    # while q_max * RG > 1.3 or x[q_min_index] * RG < 0.9:
    #     if q_max * RG > 1.3:
    #         q_max_index -= 1
    #     if x[q_min_index] * RG < 0.9:
    #         q_min_index += 1
    #     slope, intercept, rvalue, pvalue, stderr = stats.linregress(
    #         x[q_min_index:q_max_index] ** 2, np.log(y[q_min_index:q_max_index]))
    #     RG = np.sqrt(slope * -3)
    return q_min_index, q_max_index-1, slope, intercept, rvalue, pvalue, stderr


def guassian_peak_fit(xdata, ydata):
    def Gauss(x, A, B):
        y = A*np.exp(-1*B*x**2)
        return y
    
    parameters, covariance = curve_fit(Gauss, xdata, ydata, method='lm')
    print(parameters)
    print('X data')
    print(xdata)
    print('Y data')
    print(ydata)
    fit_A = parameters[0]
    fit_B = parameters[1]

    fit_y = Gauss(xdata, fit_A, fit_B)
    plt.plot(xdata, ydata, 'o', label='data')
    plt.plot(xdata, fit_y, '-', label='fit')
    



def check_filename_repetition(output_filename, directory):
    """
    :param output_filename: string, output filename
    :return: string, new output filename
    """
    print("\n==============================")
    print('Check filename repetition')
    print("------------------------------")
    files = Path(directory).glob(f'*.png')
    png_list = []
    for index, file in enumerate(files):
        png_list.append(file.name[:-4])

    print(output_filename)
    while output_filename in png_list:
        output_filename = output_filename + '_1'
        print(output_filename)
    return output_filename


def out_file(q, intensity, filename):
    """
    :param q: Array, an array stores 2theta
    :param intensity: Array, an array stores intensity
    :param filename: List, a list stores filenames
    :return: None
    """
    print('=================================================================================')
    # Shorten the filename for q output
    parts = filename.split('_')
    base_name = filename[:filename.find(FILENAME_KEYWORD)]
    short_filename = f"{base_name}_{parts[-6]}_{parts[-5]}_{parts[-4]}_{parts[-2]}_q.xy"

    print(f'Converting CMS GIWAXS data to --> {short_filename}')
    output_filename = OUTPUT_PATH / short_filename
    try:
        with open(output_filename, 'w') as out:
            out.write('q I(q)\n')
            for i in range(len(q)):
                out.write(str('{:.5f}'.format(q[i]))+' '+str('{:.5f}'.format(intensity[i]))+'\n')
    except Exception as e:
        print(f"{Style.BRIGHT}{Fore.RED}Shorten your ini filename or update FILENAME_KEYWORD to avoid the error:{Style.RESET_ALL} {e}")

    # Save as image
    if PATTERN=='.xy':
        # Use the same base to build png name, and consistently save to OUTPUT_PATH
        short_filename_png = f"{base_name}_{parts[-6]}_{parts[-5]}_{parts[-4]}_{parts[-2]}.png"
        output_file_png = Path(INPUT_PATH) / short_filename_png

        plt.figure()
        plt.plot(q, intensity)
        plt.xlabel('q ($\mathregular{\AA}^{-1}$)')
        plt.ylabel('I(q)')
        plt.title(f'{short_filename_png[:-4]}')  # Title for context
        plt.savefig(output_file_png, bbox_inches='tight')
        print(f"Plot saved successfully at: {output_file_png}")

    # Shorten the filename for 2theta output
    parts = filename.split('_')
    base_name = filename[:filename.find(FILENAME_KEYWORD)]
    short_filename = f"{base_name}_{parts[-6]}_{parts[-5]}_{parts[-4]}_{parts[-2]}_tth.xy"

    print(f'Converting CMS GIWAXS data to --> {short_filename}')
    output_filename = OUTPUT_PATH / short_filename
    with open(output_filename, 'w') as out:
        out.write('tth I(tth)\n')
        for i in range(len(q)):
            tth = 2 * np.arcsin(q[i] * 1.5406 / 4 / np.pi) * 180 / np.pi    # 0.9184, 1.5406 q = 4*pi*sin(theta)/lambda
            out.write(str('{:.5f}'.format(tth))+' '+str('{:.5f}'.format(intensity[i]))+'\n')
    print('=================================================================================')


if __name__ == '__main__':
    main()