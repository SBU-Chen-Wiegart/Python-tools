"""
File: Larch_XAS
Name: Cheng-Chu Chung
----------------------------------------
TODO: XAS data processing
Github source: https://github.com/xraypy/xraylarch
Color palettes for Python: https://jiffyclub.github.io/palettable/#palette-interface
"""


import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from larch.xafs import pre_edge
from larch.io import read_ascii, write_ascii, merge_groups, read_athena, create_athena
import athena_project
from numba import njit, jit
import palettable.colorbrewer.diverging as pld
import time as t
from random import randint

# Constant
FILE_TYPE = '.prj'  # ".prj" if you want to merge the scans; ".txt" if you want to plot the scans
INPUT_PATH = r'D:\Research data\Conversion coating\202203\20211029 BMM\merge_test'

# Merge Constant
IF_NOR = False   # Do normalization
SHOW_DATA_INFORMATION = False   # List athena parameters, such as atomic symbol, edge, label, etc.

# Plot Constant
"""
You could set FILE_INDEX = 0, SAMPLE_LIST = [], STANDARD_LIST = [], 
SAMPLE_LABEL = [], ENERGY_RANGE = () as a default for your first try.
"""
FILE_INDEX = 7  # Which file in file list you want to plot
SAMPLE_LIST = [7, 8]     # [] for default or [1, 7, 5, 3] for a index list you want to plot
STANDARD_LIST = []      # [] if no standards in the SAMPLE_LIST or [5, 3] in the SAMPLE_LIST become dash lines
FIGURE_SIZE = (6.4, 4.8)  # Cheng-hung uses (6, 7.5), but the default is (6.4, 4.8)
SAMPLE_LABEL = ['ISS twin Cu foil', 'ISS twin Cu10-PAA50']  # [] for default or adding a specific name list.
PALETTE = pld.Spectral_4_r  # _r if you want to reverse the color sequence
CMAP = PALETTE.mpl_colormap     # .mpl_colormap attribute is a continuous, interpolated map
OFFSET = 0  # Value you want to add to an offset for each curve.
ENERGY_RANGE = (17900, 18301)   # () for default
ENERGY_INTERVAL = 100   # This parameter work when you set a ENERGY_RANGE
IF_SAVE = True
OUTPUT_FILENAME = 'test_plot'


def main():
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')

    if FILE_TYPE == '.prj':
        # new_merge_project = create_athena(f'default.prj')     # Call the athena_project.py in larch.io
        new_merge_project = athena_project.create_athena(f'default.prj')    # Call the athena_project.py in current folder
        filename = ''
        for index, file_prj in enumerate(files):
            if 'Created' not in file_prj.name:
                merge_scan(file_prj, new_merge_project)
                filename = file_prj.name
        new_merge_project.save(f'{Path(INPUT_PATH)}/Created_{filename[:filename.find("b") + 3]}.prj')
        print('\n=================================================================================')
        print(f'Save merge project into ---> Created_{filename[:filename.find("b") + 3]}.prj')
        print('=================================================================================')

    elif FILE_TYPE == '.txt':
        plot_xas(files)


def plot_xas(files):
    f_list = []
    print("==============================")
    print('Files')
    print("------------------------------")
    for index, file in enumerate(files):
        f_list.append(file)
        print(index, file)

    print("==============================")
    print(f'Data column in file number {FILE_INDEX}')
    print("------------------------------")
    file = read_ascii(f_list[FILE_INDEX])
    file_keys = file.__dir__()  # Make a list
    for index, key in enumerate(file_keys):
        print(index, key)

    energy = getattr(file, file_keys[6])

    # Do the plotting
    f1, ax1 = plt.subplots(1, 1, figsize=FIGURE_SIZE)
    print("==============================")
    print('Index   Filename')
    print("------------------------------")
    increment = 0   # Increment for offset
    if len(SAMPLE_LIST) == 0:
        color_idx = np.linspace(0, 1, len(file_keys)-7)   # All plots have their own color
        for i in range(len(file_keys)-7):   # Start from the first sample name because file key 0-6 are data information
            sample_index = i + 7
            sample_name = file_keys[sample_index]
            if len(SAMPLE_LABEL) > i:
                sample_label = SAMPLE_LABEL[i]
            else:
                sample_label = sample_name
            mu = getattr(file, file_keys[sample_index])
            ax1.plot(energy, mu + OFFSET * increment, color=CMAP(color_idx[increment]), label=sample_label)
            increment += 1
            print('{:>3}     {}'.format(sample_index, sample_name))
    else:
        color_idx = np.linspace(0, 1, len(SAMPLE_LIST))     # Only the plots you want have their own color
        for sample_index in SAMPLE_LIST:
            sample_name = file_keys[sample_index]
            if len(SAMPLE_LABEL) > SAMPLE_LIST.index(sample_index):
                sample_label = SAMPLE_LABEL[SAMPLE_LIST.index(sample_index)]
            else:
                sample_label = sample_name
            mu = getattr(file, file_keys[sample_index])
            if sample_index in STANDARD_LIST:
                ax1.plot(energy, mu + OFFSET * increment, '--', color=CMAP(color_idx[increment]), label=sample_label)
            else:
                ax1.plot(energy, mu + OFFSET * increment, color=CMAP(color_idx[increment]), label=sample_label)
            increment += 1
            print('{:>3}     {}'.format(sample_index, sample_name))

    # Plot format
    if ENERGY_RANGE == ():
        ax1.set_xlim(energy.min() // 1 + 1, energy.max() // 1 - 1)
        plt.xticks(fontsize=14)
    else:
        ax1.set_xlim(ENERGY_RANGE)
        plt.xticks(np.arange(ENERGY_RANGE[0], ENERGY_RANGE[1], step=ENERGY_INTERVAL), fontsize=14)
    plt.title(OUTPUT_FILENAME, fontsize=20)
    x_label = r'$\mathregular{Energy\ (eV)}$'
    y_label = r'$\mathregular{Normalized\ \mu(E)}$'
    plt.yticks(fontsize=14)
    ax1.set_xlabel(x_label, fontsize=18)
    ax1.set_ylabel(y_label, fontsize=18)
    plt.legend(loc='lower right', framealpha=1, frameon=False)
    plt.tight_layout()
    if IF_SAVE:
        plt.savefig("{}/{}.png".format(Path(INPUT_PATH), OUTPUT_FILENAME), dpi=300, transparent=False)
    plt.show()


def merge_scan(file_prj, new_merge_project):
    """
    :param file_prj: prj file, a prj file from the current folder
    :param new_merge_project: prj file, an empty prj file to store merged data
    :return: None
    """
    filename = file_prj.name
    print(f'\nWhether the file "{filename}" exists?', file_prj.exists())

    # Read Athena project
    scans = read_athena(f'{file_prj}')
    # print(scans.groups)
    scans_namelist = []
    scans_grouplist = []    # Each scan is a group

    # Print each scan name in the prj file
    print("\n==============================")
    print('Scan name')
    print("------------------------------")
    for name, group in scans._athena_groups.items():
        scans_namelist.append(name)
        scans_grouplist.append(group)
        print(name, group)

    # Print scan information
    first_scan_information = scans_grouplist[0]
    if SHOW_DATA_INFORMATION:
        print("\n==============================")
        print(f'Athena parameters in {scans_namelist[0]}')
        print("------------------------------")
        for scan_attribute in dir(first_scan_information):
            print(scan_attribute, type(getattr(first_scan_information, scan_attribute)))

    # Plot scans
    for index, scan in enumerate(scans_grouplist):
        if IF_NOR:  # Do normalization
            pre_edge(scan.energy, scan.mu, group=scan)
            plt.plot(scan.energy, scan.norm, label=scan.label)
        else:
            plt.plot(scan.energy, scan.mu, label=scan.label)

    # Do the merge and plot
    merges = merge_groups(scans_grouplist)
    if IF_NOR:  # Do normalization
        pre_edge(merges.energy, merges.mu, group=merges)
        plt.plot(merges.energy, merges.norm, label=f'{first_scan_information.label[:-4]}_merged')
    else:
        plt.plot(merges.energy, merges.mu, label=f'{first_scan_information.label[:-4]}_merged')
    if SHOW_DATA_INFORMATION:
        print("\n==============================")
        print(f'Merged parameters in {first_scan_information.label[:-4]}_merged')
        print("------------------------------")
        for scan_attribute in dir(merges):
            print(scan_attribute, type(getattr(merges, scan_attribute)))

    # Replace '-' with '_' because '-' will cause error and add the merge into the new prj
    scan_name = f'{first_scan_information.label[:-4]}_merged'.replace('-', '_')
    new_merge_project.add_group(merges, scan_name)

    # Figure information
    plt.xlabel('$\mathregular{Energy\ (eV)}$', fontsize=12)
    plt.ylabel('$\mathregular{Normalized\ \mu(E)}$', fontsize=12)
    plt.title(f'{first_scan_information.atsym} {first_scan_information.edge}-edge')
    plt.legend(loc='lower right', framealpha=1, frameon=False)
    if IF_SAVE:
        plt.savefig("{}/{}.png".format(Path(INPUT_PATH), first_scan_information.label[:-4]), dpi=300, transparent=False)
        print('\n=================================================================================')
        print(f'Save merge image into ---> {first_scan_information.label[:-4]}.png')
        print('=================================================================================')
    plt.close()


if __name__ == '__main__':
    main()
