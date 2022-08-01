"""
File: Bar_Chart.py
Name: Cheng-Chu Chung
----------------------------------------
TODO: Bar chart in publication quiality
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import palettable.colorbrewer.diverging as pld
import matplotlib.lines as lines

PALETTE = pld.Spectral_4_r  # _r if you want to reverse the color sequence
CMAP = PALETTE.mpl_colormap     # .mpl_colormap attribute is a continuous, interpolated map


def main():
    fig, ax = plt.subplots()
    colors = ["#69b3a2", "#4374B3"]
    # sns.set_palette(sns.color_palette(colors))

    samples = ['NbAl\npristine', 'Sc/NbAl\npristine', 'Sc/NbAl\n900C60M', 'Sc/NbAl\n1100C60M']
    element_A_at = [53.00, 51.23, 50.95, 50.52]
    elememt_B_at = [47.00, 48.77, 49.05, 49.48]
    x = np.arange(len(samples))
    width = 0.4
    color_idx = np.linspace(0, 1, len(samples))

    bar_A = plt.bar(x, element_A_at, width, color=CMAP(color_idx[0]), label='Nb', edgecolor='white')
    bar_B = plt.bar(x + width, elememt_B_at, width, color=CMAP(color_idx[1]), label='Al', edgecolor='white')

    for index, bar in enumerate(bar_A):
        height_A = bar.get_height()
        height_B = 100 - height_A

        ax.text(bar.get_x() + bar.get_width() / 2, 1.05 * height_A,
                f'{element_A_at[index]}',
                ha='center', va='bottom', rotation=0, fontsize=12)
        ax.text(bar.get_x() + bar.get_width() / 2 + width, 1.05 * height_B,
                f'{elememt_B_at[index]}',
                ha='center', va='bottom', rotation=0, fontsize=12)

    spineline = ['left', 'right', 'top', 'bottom']
    for direction in spineline:
        ax.spines[direction].set_linewidth('2')
    plt.yticks(fontsize=14)
    ax.tick_params(width=2)
    plt.xticks(x + width / 2, samples, fontsize=14)
    plt.ylabel('Atomic percentage', fontsize=14)
    plt.ylim(0, 80)
    plt.title('')
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()