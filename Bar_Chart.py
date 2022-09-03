"""
File: Bar_Chart.py
Name: Cheng-Chu Chung
----------------------------------------
TODO: Bar chart in publication quality
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import palettable.colorbrewer.diverging as pld
import matplotlib.lines as lines
import configparser
from collections import defaultdict
from pprint import pprint

PALETTE = pld.Spectral_4_r  # _r if you want to reverse the color sequence
CMAP = PALETTE.mpl_colormap     # .mpl_colormap attribute is a continuous, interpolated map
CONFIG_FILE = r'D:\Research data\SSID\Advanced Computer Python\Python-tools\bar_chart_config.ini'


def main():
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE, encoding="utf8")
    print('\n====================')
    print(f'Sections')
    print('====================')
    pprint(config.sections())

    fig, ax = plt.subplots()
    colors = ["#69b3a2", "#4374B3"]
    # sns.set_palette(sns.color_palette(colors))

    print('\n====================')
    print(f'Sample name')
    print('====================')
    sample_list = []
    for item in config['samples']['sample_name'].split(', '):
        print(item)
        sample_list.append(item.replace(' ', '\n'))
    samples = sample_list

    print('\n====================')
    print(f'Legends')
    print('====================')
    element_dictionary = defaultdict(list)
    for element in config['quantity']:
        element_dictionary[element] = list(np.array(config['quantity'][element].split(', ')).astype(float))
        print(f"{element}: {config['legends'][element]}")

    print('\n====================')
    print(f'Quantity')
    print('====================')
    pprint(element_dictionary)

    x = np.arange(len(samples))
    bar_width = config.getfloat('format', 'bar_width')
    bar_width = 0.8/len(config['legends'])   # A relation in serendipity
    color_idx = np.linspace(0, 1, len(samples))

    # Plotting
    first_bar = None
    for index, bar in enumerate(element_dictionary):
        if index == 0:
            first_bar = plt.bar(x + bar_width * index, element_dictionary[bar], bar_width, color=CMAP(color_idx[index]),
                label=config['legends'][bar], edgecolor=config['format']['bar_edgecolor'])
        else:
            plt.bar(x + bar_width * index, element_dictionary[bar], bar_width, color=CMAP(color_idx[index]),
                    label=config['legends'][bar], edgecolor=config['format']['bar_edgecolor'])

    # Texts
    if config.getint('text_on_bar', 'text_on'):
        for bar_index, bar in enumerate(first_bar):
            for element_index, element in enumerate(element_dictionary):
                height_bar = element_dictionary[element][bar_index]
                ax.text(bar.get_x() + bar.get_width() / 2 + bar_width * element_index,
                        config.getfloat('text_on_bar', 'text_height_factor') * height_bar, f'{height_bar}',
                        ha=config['text_on_bar']['horizontalalignment'], va=config['text_on_bar']['verticalalignment'],
                        rotation=config.getint('text_on_bar', 'rotation'),
                        fontsize=config.getint('text_on_bar', 'fontsize'))

    # Frame linewidth
    spineline = ['left', 'right', 'top', 'bottom']
    for direction in spineline:
        ax.spines[direction].set_linewidth(config['format']['spinelinewidth'])

    # Formatting
    plt.yticks(fontsize=config.getint('format', 'fontsize'))
    ax.tick_params(width=config.getint('format', 'tick_width'))
    plt.xticks(x + bar_width / 2 * (len(config['legends'])-1), samples, fontsize=config.getint('format', 'fontsize'))
    plt.ylabel(config['format']['ylabel'], fontsize=config.getint('format', 'fontsize'))
    plt.ylim(tuple(np.array(config['format']['ylim'].split(', ')).astype(int)))
    plt.title(config['format']['title'], fontsize=config.getint('format', 'title_fontsize'),
              pad=config.getint('format', 'title_pad'))
    plt.legend(fontsize=config.getint('format', 'fontsize'), ncol=config.getint('format', 'legend_column'),
               frameon=config.getint('format', 'frameon'))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()