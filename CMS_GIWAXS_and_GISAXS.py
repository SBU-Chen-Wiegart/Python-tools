"""
File: CMS GIWAXS AND GISAXS
Name: Cheng-Chu Chung
TODO: Auto import ref peaks from PDF card data, then put data into dictionary
----------------------------------------
"""

from pathlib import Path
from collections import defaultdict
import pandas as pd
from pprint import pprint
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy import signal
import os
import configparser
import palettable as pltt
import peakutils

INPUT_PATH = r'D:\Research data\SSID\202210\20221003 CMS b32\waxs\analysis\circular_average'
CONFIG_FILE = r"D:\Research data\SSID\202210\20221003 CMS b32\waxs\CMS_plot_config.ini"
CONFIG = configparser.ConfigParser()

if Path(CONFIG_FILE).is_file():
    CONFIG.read(CONFIG_FILE, encoding="utf8")
    print("-----------------")
    print("Use .ini input")
    print("-----------------")
else:
    print("-----------------")
    print("Manually input so please remove all eval commands if error occurs or file is not found")
    print("-----------------")

FILE_TYPE = '.dat'
SAMPLE_LIST = eval(CONFIG['samples']['sample_list'])
BATCH_NUMBER, COMPOSITION, CONDITION, INCIDENT_ANGLE = eval(CONFIG['legends']['sample_condition'])   # Whether you want to show them in the legend
PALETTE = eval(CONFIG['format']['palette'])                                                          # pld.Spectral_4_r  # _r if you want to reverse the color sequence
CMAP = PALETTE.mpl_colormap                                                                          # .mpl_colormap attribute is a continuous, interpolated map
COLOR_INCREMENT = eval(CONFIG['format']['color_increment'])                                          # Adjust your color gradient, 0 for default
OFFSET = eval(CONFIG['format']['offset'])                                                            # Value you want to add to an y offset for each curve.
SAMPLE_LABEL = eval(CONFIG['legends']['sample_label'])
TITLE = eval(CONFIG['format']['output_filename'])
OUTPUT_FOR_JADE = eval(CONFIG['format']['output_for_jade'])
IF_SAVE = eval(CONFIG['format']['if_save'])
SUB_DEGREE = eval(CONFIG['data_processing']['bgsub_degree'])

def main():
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')
    for index, dat_file in enumerate(files):
        print(index, dat_file.name)
    print('-----------------------------------------')
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')  # Call Path again to grab the file
    giwaxs(files)


def sorted_data(files):
    data_dict = {'q_list': {}, 'I_list': {},'filename_list': {}, 'background_subtraction_list': {}}   # Create a dictionary to store all information
    for index, file in enumerate(files):
        waxs_data = file.resolve()  # Make the path absolute, resolving any symlinks
        data_dict['filename_list'][index] = waxs_data.name
        # print(index, waxs_data.name)
        dataframe = pd.read_table(waxs_data, sep="\s+",
                                  usecols=['#', 'q', 'qerr', 'I(q)'])\
            .to_numpy()  # '#' is q column and 'qerr' is I(q) column
        data_dict['q_list'][index] = dataframe[:, 0]
        data_dict['I_list'][index] = dataframe[:, 2]
        # if index == 0:
        #     pprint(data_dict)
        if OUTPUT_FOR_JADE:
                out_file(dataframe[:, 0], dataframe[:, 2], f'Converted_{waxs_data.name}')
    return data_dict


def giwaxs(files):
    q_and_I_list = sorted_data(files)
    background_subtraction(q_and_I_list, degree=5)   # Add a new list with background subtraction
    intensity_plot(q_and_I_list)
    intensity_plot(q_and_I_list, mode='bg_sub')


def background_subtraction(q_and_I_list, degree=3):
    """
    Subtract background from data and save data w/ background subtracted
    :return: None
    """
    # list_dict = sorted_data()
    for index in q_and_I_list['q_list']:
        x = q_and_I_list['q_list'][index]
        y = q_and_I_list['I_list'][index]
        filename = q_and_I_list['filename_list'][index]
        # pos = filename[filename.find('pos'):filename.find('pos') + 4]
        base_line = peakutils.baseline(y, degree)  # generate baseline
        y_corrected = y - base_line  # subtract baseline from y-values
        q_and_I_list['background_subtraction_list'][index] = y_corrected
        if OUTPUT_FOR_JADE:
                out_file(x, y_corrected, f'Converted_bgsub_{filename}')


def intensity_plot(q_and_I_list, mode='raw'):
    fig, ax = plt.subplots()
    title = TITLE
    y_max_lim = 0
    y_min_lim = 10000
    increment = 0
    for index in SAMPLE_LIST:
        x = q_and_I_list['q_list'][index]
        if mode == 'raw':
            y = q_and_I_list['I_list'][index]
        elif mode == 'bg_sub':
            y = q_and_I_list['background_subtraction_list'][index]
        filename = q_and_I_list['filename_list'][index]
        
        # Plot label
        curve_info = filename[filename.find('_') + 1:filename.find('_', 3)]
        batch_number = filename[filename.find('b'):filename.find('b') + 6] if BATCH_NUMBER else '..'
        composition = filename[filename.find('b') + 7:filename.find('-', 10)] if COMPOSITION else '..'
        condition = filename[filename.find('-', 10) + 1:filename.find('_pos')] if CONDITION else '..'
        incident_angle = f"{filename[filename.find('th') + 2:filename.find('th') + 6]} degree" if INCIDENT_ANGLE else '..'
        if len(SAMPLE_LABEL) == 0:
            plot_label = f"{batch_number}/{composition}/{condition}/{incident_angle}"
        else:
            plot_label = SAMPLE_LABEL[SAMPLE_LIST.index(index)]
        # Title
        if TITLE != 'Auto':
            title = TITLE
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

    # Frame linewidth
    spineline = ['left', 'right', 'top', 'bottom']
    for direction in spineline:
        ax.spines[direction].set_linewidth(2)

    # Plotting format
    plt.xlabel('q ($\mathregular{\AA}^{-1}$)', fontsize=18)
    plt.ylabel('I(q)', fontsize=18, labelpad=10)
    plt.xticks(fontsize=14)
    plt.yticks([])  # Disable ticks
    ax.tick_params(width=2)
    plt.xlim(q_and_I_list['q_list'][SAMPLE_LIST[0]].min(), q_and_I_list['q_list'][SAMPLE_LIST[0]].max())
    plt.ylim(y_min_lim-50, y_max_lim+200)
    plt.legend(loc='upper left', framealpha=1, frameon=False, fontsize=12)
    plt.title(title, fontsize=18)
    plt.tight_layout()
    if IF_SAVE:
        output_filename = check_filename_repetition(title)
        plt.savefig("{}/{}.png".format(Path(INPUT_PATH), output_filename), dpi=300, transparent=False)
    plt.show()


def check_filename_repetition(output_filename):
    """
    :param output_filename: string, output filename
    :return: string, new output filename
    """
    print("\n==============================")
    print('Check filename repetition')
    print("------------------------------")
    files = Path(INPUT_PATH).glob(f'*.png')
    png_list = []
    for index, file in enumerate(files):
        png_list.append(file.name[:-4])

    print(output_filename)
    while output_filename in png_list:
        output_filename = output_filename + '_1'
        print(output_filename)
    return output_filename


def out_file(tth, intensity, filename):
    """
    :param tth: Array, an array stores 2theta
    :param intensity: Array, an array stores intensity
    :param filename: List, a list stores filenames
    :return: None
    """
    print('=================================================================================')
    short_filename = filename[:filename.find('_pos1')] + '-' + filename[filename.find('th'):filename.find('th')+6]
    print(f'Converting CMS GIWAXS data to --> {short_filename}')
    filename = os.path.join(INPUT_PATH+'/', short_filename + '.xy')
    with open(filename, 'w') as out:
        out.write('tth intensity\n')
        for i in range(len(tth)):
            out.write(str('{:.2f}'.format(tth[i]))+' '+str('{:.5f}'.format(intensity[i]))+'\n')
    print('=================================================================================')
    print(' ')


if __name__ == '__main__':
    main()