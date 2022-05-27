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

# Step 1: Paste your directory
INPUT_PATH = 'D:/Research data/SSID/202205/20220526 XRD b31 NbAl'   # <-- Enter the folder position you want to explore

# Step 2: Set up your plotting parameters
# CONSTANT
PLOT_LIST = [11, 9, 7, 5]    # [] for default or [1, 7, 5, 3] index list for the index sequence you desire
OUTPUT = False   # "True" if you want to save the converted file
PLOT_OFFSET = 500    # Number you want to add to an offset for each curve.
PLOT_FIGURE = True  # "True" if you want to show the plots
IF_LEGEND = True    # "True" if you want to show the legend
# Good luck for your data conversion!


def main():
    intensity_plot()


# @pysnooper.snoop()
def convert_format():
    """
    :return: dict, a dictionary stores 2theta, intensity, and filename
    """
    list_dict = {'q_list': {}, 'I_list': {}, 'filename_list': {}}   # Create a dictionary to store all information
    # Import all the data
    for home, dirs, files in os.walk(INPUT_PATH):
        print('This is home: ' + home)
        print(f'This is dirs: {dirs}')
        print('Index and filename')
        for index, filename in enumerate(files):
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


def intensity_plot():
    """
    :return: None
    """
    list_dict = convert_format()   # Import the data
    if len(PLOT_LIST) == 0:
        index = list_dict['filename_list']     # Select the index from the list_dict['filename_list']
    else:
        index = PLOT_LIST
    plot_sequence = 0
    print('Plot:')
    fig, ax = plt.subplots()
    for i in index:
        x = list_dict['q_list'][i]
        y = list_dict['I_list'][i] + plot_sequence*PLOT_OFFSET
        filename = list_dict['filename_list'][i]
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
