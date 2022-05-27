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
from larch.io import read_ascii, write_ascii, merge_groups, read_athena
from numba import njit, jit
import palettable.colorbrewer.diverging as pld

FILE_TYPE = '.txt'
INPUT_PATH = r'D:\Research data\Conversion coating\202205\20220511 ISS samples'
FILE_INDEX = 6
OFFSET = 0
SAMPLE_LIST = [7, 8, 9, 10, 11, 12, 13, 14, 15]     # [] for default or [1, 7, 5, 3] for index list you want
STANDARD_LIST = [9, 11, 13, 15]
IF_SAVE = False
OUTPUT_FILENAME = 'merge_mu_01.png'


def main():
    # plt.close('all')
    palette = pld.Spectral_4_r
    cmap = palette.mpl_colormap

    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')
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
            ax1.plot(energy, mu + OFFSET * increment, color=cmap(color_idx[increment]), label=sample_name)
            increment += 1
            print('{:>3}     {}'.format(sample_index, sample_name))
    else:
        color_idx = np.linspace(0, 1, len(SAMPLE_LIST))     # Only the plots you want have their own color
        for sample_index in SAMPLE_LIST:
            sample_name = f0_keys[sample_index]
            mu = getattr(f0, f0_keys[sample_index])
            if sample_index in STANDARD_LIST:
                ax1.plot(energy, mu + OFFSET * increment, '--', color=cmap(color_idx[increment]), label=sample_name)
            else:
                ax1.plot(energy, mu + OFFSET * increment, color=cmap(color_idx[increment]), label=sample_name)
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
    plt.draw()
    plt.show()


def merge_scan():
    """
    TODO: import list
    """
    file_path = r'D:\Research data\Conversion coating\202203\20211029 BMM\Cu-b26_coating_on_Al'
    file_folder = Path(file_path)
    file_name = 'Cu-b26-4_Al_Cu20_PAMAM.prj'
    file_prj = file_folder/file_name

    # files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')
    # f_list = []
    # print("==============================")
    # print('Files')
    # print("------------------------------")
    # for index, file in enumerate(files):
    #     f_list.append(file)
    #     print(index, file)

    print(f'Whether the file ({file_name}) exists?', file_prj.exists())
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

    print("\n==============================")
    print('Athena parameters')
    print("------------------------------")
    for scan_attribute in dir(scans_grouplist[0]):
        print(scan_attribute, type(getattr(scans_grouplist[0], scan_attribute)))

    for index, scan in enumerate(scans_grouplist):
        pre_edge(scan.energy, scan.mu, group=scan)   # Do normalization
        plt.plot(scan.energy, scan.norm, label=scan.label)

    merges = merge_groups(scans_grouplist)
    pre_edge(merges.energy, merges.mu, group=merges)  # Do normalization
    plt.plot(merges.energy, merges.norm, label=f'{scans_grouplist[0].label[:-4]}_merged')

    plt.xlabel('Energy')
    plt.ylabel('mu')
    plt.title(f'{scans_grouplist[0].atsym} {scans_grouplist[0].edge}-edge')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
    # merge_scan()