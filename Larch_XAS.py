"""
File: Larch_XAS
Name: Cheng-Chu Chung
----------------------------------------
TODO: XAS data processing
https://github.com/xraypy/xraylarch
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
FILE_TYPE = '.prj'
INPUT_PATH = r'D:\Research data\Conversion coating\202203\20211029 BMM\merge_test'
IF_NOR = False   # Do normalization
SHOW_DATA_INFORMATION = False   # List athena parameters, such as atomic symbol, edge, label, etc.

FILE_INDEX = 6
OFFSET = 0
SAMPLE_LIST = [7, 8, 9, 10, 11, 12, 13, 14, 15]     # [] for default or [1, 7, 5, 3] for index list you want
STANDARD_LIST = [9, 11, 13, 15]
IF_SAVE = True
OUTPUT_FILENAME = 'merge_mu_01.png'
PALETTE = pld.Spectral_4_r
CMAP = PALETTE.mpl_colormap


def main():
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')
    # plot_xas(files)
    # new_merge_project = create_athena(f'default.prj')     # Call the athena_project.py in larch.io
    new_merge_project = athena_project.create_athena(f'default.prj')    # Call the athena_project.py in current folder
    filename = ''
    for index, file_prj in enumerate(files):
        if 'Created' not in file_prj.name:
            merge_scan(file_prj, new_merge_project)
            filename = file_prj.name
    new_merge_project.save(f'{Path(INPUT_PATH)}/Created_{filename[:filename.find("b") + 3]}.prj')


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
    f0 = read_ascii(f_list[FILE_INDEX])
    f0_keys = f0.__dir__()  # Make a list
    for index, item in enumerate(f0_keys):
        print(index, item)

    energy = getattr(f0, f0_keys[6])

    f1, ax1 = plt.subplots(1, 1, figsize=(6, 7.5))
    print("==============================")
    print('Index   Filename')
    print("------------------------------")
    increment = 0
    if len(SAMPLE_LIST) == 0:
        color_idx = np.linspace(0, 1, len(f0_keys)-7)   # All plots have their own color
        for i in range(len(f0_keys)-7):
            sample_index = i + 7
            sample_name = f0_keys[sample_index]
            mu = getattr(f0, f0_keys[sample_index])
            ax1.plot(energy, mu + OFFSET * increment, color=CMAP(color_idx[increment]), label=sample_name)
            increment += 1
            print('{:>3}     {}'.format(sample_index, sample_name))
    else:
        color_idx = np.linspace(0, 1, len(SAMPLE_LIST))     # Only the plots you want have their own color
        for sample_index in SAMPLE_LIST:
            sample_name = f0_keys[sample_index]
            mu = getattr(f0, f0_keys[sample_index])
            if sample_index in STANDARD_LIST:
                ax1.plot(energy, mu + OFFSET * increment, '--', color=CMAP(color_idx[increment]), label=sample_name)
            else:
                ax1.plot(energy, mu + OFFSET * increment, color=CMAP(color_idx[increment]), label=sample_name)
            increment += 1
            print('{:>3}     {}'.format(sample_index, sample_name))

    ax1.set_xlim(energy.min()//1+1, energy.max()//1-1)
    x_label = r'$Energy\ (eV)$'  # r'$Capacity\ (mAh g^{-1})$'
    y_label = r'$Normalized\ \mu(E)$'  # r'$Voltage\ (V)$'
    ax1.set_xlabel(x_label, fontsize=18)  # fontweight=fontweight
    ax1.set_ylabel(y_label, fontsize=18)  # fontweight=fontweight

    plt.legend(loc='lower right')
    if IF_SAVE:
        plt.savefig(OUTPUT_FILENAME, dpi=300, transparent=True)
    plt.show()


def merge_scan(file_prj, new_merge_project):
    """
    :param file_prj: a prj file
    :return: None
    """
    filename = file_prj.name
    print(f'\nWhether the file "{filename}" exists?', file_prj.exists())

    # Read Athena project
    scans = read_athena(f'{file_prj}')
    # print(scans.groups)
    scans_namelist = []
    scans_grouplist = []

    print("\n==============================")
    print('Scan name')
    print("------------------------------")
    for name, group in scans._athena_groups.items():
        scans_namelist.append(name)
        scans_grouplist.append(group)
        print(name, group)

    # Create Athena project
    first_scan_information = scans_grouplist[0]

    if SHOW_DATA_INFORMATION:
        print("\n==============================")
        print(f'Athena parameters in {scans_namelist[0]}')
        print("------------------------------")
        for scan_attribute in dir(first_scan_information):
            print(scan_attribute, type(getattr(first_scan_information, scan_attribute)))

    for index, scan in enumerate(scans_grouplist):
        if IF_NOR:  # Do normalization
            pre_edge(scan.energy, scan.mu, group=scan)
            plt.plot(scan.energy, scan.norm, label=scan.label)
        else:
            plt.plot(scan.energy, scan.mu, label=scan.label)

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

    # Replace '-' with '_' because '-' will cause error
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
