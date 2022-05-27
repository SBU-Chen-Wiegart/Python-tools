"""
File: Plot_XPD_in-situ_heating
Name: Cheng-Chu Chung
----------------------------------------
TODO: Plot XPD in-situ heating
"""
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import palettable.colorbrewer.sequential as pbs
import palettable.cartocolors.diverging as pcd
import palettable.cartocolors.sequential as pcs
import palettable.scientific.sequential as pss
import time as t
import pysnooper

# Constant
INPUT_PATH = r'D:\Research data\SSID\202204\20220406 XPD ex-situ check\LK_b30-14_Nb40Al60Sc_SiO2Si_pristine_heating'
TIMESTAMP_FILE = 'sample_LK_heating_20220408_172022.xlsx'   # '' will return index as the color bar
FILE_TYPE = '.xy'
HEADERS = ['Q or 2theta', 'intensity']
SKIPROWS = 23
PLOT_LIST = list(np.arange(0, 191, 5))   # [] for default or [1, 7, 5, 3] index list for the index sequence you desire
COLORBAR_TICKS = 10
PALETTE = pcd.Earth_7
CMAP = PALETTE.mpl_colormap
PLOT_OFFSET = 0.5    # Number you want to add to an offset for each curve.
PLOT_FIGURE = True  # "True" if you want to show the plots
SAVE_IMG = False   # "True" if you want to save the converted file


# @pysnooper.snoop()
def main():
    data_info_list = file_preprocessing()
    data_number = len(data_info_list['filename_list'])
    timestamp_info = read_timestamp(data_number)
    plot_data(data_info_list, timestamp_info)


def file_preprocessing():
    list_dict = {'x_list': {}, 'y_list': {}, 'filename_list': {}}
    # Path function converts \ to / and glob method returns .xy files in a generator ---> Very powerful!
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')
    number_of_files = len(list(files))   # Path function ends
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')      # Call Path again to execute the for-loop
    print('Index Filename')
    for index, file_directory in enumerate(files):
        file = file_directory.resolve()  # Make the path absolute, resolving any symlinks
        filename = file.name
        print('{0:<5} {1}'.format(index, filename))
        list_dict['filename_list'][index] = filename
        data = pd.read_table(file,
                           delimiter='\s+',
                           engine='python',
                           skiprows=SKIPROWS,
                           names=HEADERS)
        x = np.array(data[HEADERS[0]].tolist())  # q
        y = np.array(data[HEADERS[1]].tolist())  # I(q)
        list_dict['x_list'][index] = x
        list_dict['y_list'][index] = y
        if index == number_of_files-1:
            print('=================================')
            print(data)
    return list_dict


def plot_data(data_info_list, timestamp_info):
    if len(PLOT_LIST) == 0:
        index = data_info_list['filename_list']  # Select the index from the list_dict['filename_list']
    else:
        index = PLOT_LIST
    plot_sequence = 0
    colors = [CMAP(i) for i in np.linspace(0, 1, len(index))]   # Create a color bar scale
    print('=================================')
    print('Plot:')

    # Create a color bar based on the number of data you want to plot
    if TIMESTAMP_FILE == '':
        time_interval = 60
    else:
        time = timestamp_info['time'].tolist()
        time.reverse()
        time_interval = (time[1] - time[0]).seconds
    color_x = range(len(index))
    color_y = range(len(index))
    color_scale = np.array(index)*time_interval/60
    color_bar = plt.scatter(color_x, color_y, c=color_scale, cmap=CMAP)
    plt.close()

    # Start to plot
    fig, ax = plt.subplots(linewidth=50)
    for i in index:
        x = data_info_list['x_list'][i]
        y = data_info_list['y_list'][i] + plot_sequence * PLOT_OFFSET
        filename = data_info_list['filename_list'][i]
        print(i, filename)
        plt.plot(x, y, color=colors[index.index(i)], linewidth=2)    # , label=f'{filename}'
        plot_sequence += 1

    # Plot format
    for axis in ['top', 'bottom', 'left', 'right']:     # Change all spines
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params(width=1.5)                           # Increase tick width
    plt.xlabel('$\mathregular{q \ (\\AA^{-1})}$', fontsize=18)
    plt.ylabel('Intensity (arb. units)', fontsize=18)
    plt.xlim(1, 10)
    # plt.ylim(1, 40)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(title=f'           Sample name\n'
                     f'b30-14_Nb40Al60Sc_SiO2Si', title_fontsize=12)
    plt.title(f'XPD in-situ heating experiment - {len(index)}/{len(data_info_list["filename_list"])}', fontsize=14)

    # Color bar setting
    ticks_setting = np.linspace(color_scale.min(), color_scale.max(), COLORBAR_TICKS, endpoint=True)
    cbar = fig.colorbar(color_bar, ticks=ticks_setting, pad=0.05)
    if TIMESTAMP_FILE == '':
        cbar.set_label('Index', fontsize=18)
    else:
        cbar.set_label('Time (mins)', fontsize=18, labelpad=10)
    cbar.ax.tick_params(labelsize=12, width=1.5)

    # Export the figure
    plt.tight_layout(pad=2)
    if SAVE_IMG:
        fig.savefig("{}/XPD_{}.png".format(Path(INPUT_PATH), t.time()))  # <--- added to save figure
    if PLOT_FIGURE:
        plt.show()


def read_timestamp(number_of_rows):
    if TIMESTAMP_FILE == '':
        pass
    else:
        file = (Path(INPUT_PATH) / f'{TIMESTAMP_FILE}').resolve()   # Make the path absolute, resolving any symlinks
        data = pd.read_excel(file)[:number_of_rows]
        return data


if __name__ == '__main__':
    main()
