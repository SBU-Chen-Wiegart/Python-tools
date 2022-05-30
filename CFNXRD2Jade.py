"""
File: CFNXRD2Jade
Name: Cheng-Chu Chung
----------------------------------------
TODO: Convert CFN XRD data to the data readable in Jade (10 --> 10.00 for x_axis value)
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

# Step 1: Paste your directory
INPUT_PATH = r'D:\Research data\SSID\202205\20220526 XRD b31 NbAl'   # <-- Enter the folder position you want to explore

# Step 2: Set up your plotting parameters
# CONSTANT
FILE_TYPE = '.xy'
PLOT_LIST = [5, 4, 3, 2]    # [] for default or [1, 7, 5, 3] index list for the index sequence you desire
OUTPUT = False   # "True" if you want to save the converted file
PLOT_OFFSET = 500    # Value you want to add to an offset for each curve.
PLOT_FIGURE = True  # "True" if you want to show the plots
IF_LEGEND = True    # "True" if you want to show the legend
# Good luck for your data conversion!


def main():
    # Path function converts \ to / and glob method returns .xy files in a generator ---> Very powerful!
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')
    dictionary_of_I_and_q = convert_format(files)
    intensity_plot(dictionary_of_I_and_q)


# @pysnooper.snoop()
def convert_format(files):
    """
    :param files: generator, list of files in the directory
    :return: dict, a dictionary stores 2theta, intensity, and filename
    """
    list_dict = {'q_list': {}, 'I_list': {}, 'filename_list': {}}   # Create a dictionary to store all information
    # Import all the data
    print('Index Filename')
    for index, file_directory in enumerate(files):
        file = file_directory.resolve()  # Make the path absolute, resolving any symlinks
        filename = file.name
        if '.xy' in filename and 'Converted' not in filename:
            list_dict['filename_list'][index] = filename    # Append index and filename for the outline
            print(index, filename)  # Print all the files
            df = pd.read_table(
                INPUT_PATH + '/' + filename, header=None
                )
            x = np.array(df[0].tolist())   # q
            y = np.array(df[1].tolist())   # I(q)
            list_dict['q_list'][index] = x
            list_dict['I_list'][index] = y
            # print(df)
            if OUTPUT:
                out_file(x, y, f'Converted_{filename}')
    return list_dict


def out_file(tth, intensity, filename):
    """
    :param tth: Array, an array stores 2theta
    :param intensity: Array, an array stores intensity
    :param filename: List, a list stores filenames
    :return: None
    """
    print('=================================================================================')
    print(f'Converting CFN XRD data to --> {filename}')
    filename = os.path.join(INPUT_PATH+'/', filename)
    with open(filename, 'w') as out:
        out.write('tth intensity\n')
        for i in range(len(tth)):
            out.write(str('{:.2f}'.format(tth[i]))+' '+str('{:.5f}'.format(intensity[i]))+'\n')
    print('=================================================================================')
    print(' ')


def intensity_plot(dictionary_of_I_and_q):
    """
    :param dictionary_of_I_and_q: dict, a dictionary contains an intensity, q, and filename list
    :return: None
    """
    # Import the data
    if len(PLOT_LIST) == 0:
        index = dictionary_of_I_and_q['filename_list']     # Select the index from the list_dict['filename_list']
    else:
        index = PLOT_LIST
    plot_sequence = 0
    print('Plot:')
    fig, ax = plt.subplots()
    for i in index:
        x = dictionary_of_I_and_q['q_list'][i]
        y = dictionary_of_I_and_q['I_list'][i] + plot_sequence*PLOT_OFFSET
        filename = dictionary_of_I_and_q['filename_list'][i]
        print(i, filename)
        plt.plot(x, y, label=f'{filename}')
        plot_sequence += 1
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    plt.xlabel('$\mathregular{2\\theta \\ (degree)}$', fontsize=12)
    plt.ylabel('Intensity (arb. units)', fontsize=12)
    plt.xlim(10, 80)
    if IF_LEGEND:
        plt.legend(title=r'Sample name')
    # plt.tight_layout()
    plt.title('CFN XRD')
    if PLOT_FIGURE:
        plt.show()


if __name__ == '__main__':
    main()
