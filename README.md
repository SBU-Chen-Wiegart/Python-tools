# Python-tools
Data Analysis Tools for X-ray Spectroscopy

## Table of Contents

- [Larch XAS.py](#larch-XASpy)
- [CFNXRD2Jade.py](#cFNXRD2Jadepy)
- [plot xrd SSID insitu.py](#plot-xrd-SSID-insitupy)

## Larch XAS.py
> Citation: M. Newville, Larch: An Analysis Package For XAFS And Related Spectroscopies. Journal of Physics: Conference Series, 430:012007 (2013).
- Documentation: http://xraypy.github.io/xraylarch
- Code: http://github.com/xraypy/xraylarch

### [Install](https://xraypy.github.io/xraylarch/installation.html)

---
* Download [Anaconda](https://www.anaconda.com/)
* Download **athena_project.py** and put it at the same folder with Larch_XAS.py

Within a shell:

1. Activate your conda environment (called base by default) and update it:
```
conda activate
conda update -y conda python pip
```
2. (optional/expert) Create a dedicated environment for Larch and activate it:
```
conda create -y --name xraylarch python==3.9
conda activate xraylarch
conda update --all
```
3. Install main dependencies:
```
conda install -y "numpy=>1.20" "scipy=>1.6" "matplotlib=>3.0" scikit-learn pandas
conda install -y -c conda-forge wxpython pymatgen tomopy pycifrw
```
4. Install Larch (latest release):
```
pip install xraylarch
```
5. Set Conda Environment by changing Python Interpreters as Anaconda\envs\xraylarch\python.exe 

### Usage

---
#### Data processing

FILE_TYPE =
* '.prj': merge fluorescence scans (_e.g._ BNL NSLS-II BMM)
* '': merge transmission scans (_e.g._ BNL NSLS-II BMM)
* '.dat': merge transmission scans (_e.g._ BNL NSLS-II ISS)
* '.txt': Plot scans file exported from Athena

INPUT_PATH = 
* Copy and paste data location

#### Data display (optional)

SKIP_SCANS = 
* ['MnO2_45_16C_Charge_Mn_001']: add scan name you want to exclude

IF_NOR = 
* True / False: do normalization for fluorescence scans

ADD_DEV =
* True / False: add plus and minus standard deviation lines for fluorescence scans

SHOW_DATA_INFORMATION = 
* True / False: display Athena parameters, such as atomic symbol, edge, label, etc.

#### Plot spectrum (optional)

* Download **larch_plot_config.ini**

CONFIG_FILE = 
* Copy and paste the absolute directory of the config file for plotting .txt files

## CFNXRD2Jade.py
### Usage

---
#### Step 1: Paste your data directory

#### Step 2: Set up your plotting parameters

## plot xrd SSID insitu.py
<p align="center">
  <img src="./XPD_time.struct_time(tm_year=2022, tm_mon=4, tm_mday=21, tm_hour=23, tm_min=9, tm_sec=11, tm_wday=3, tm_yday=111, tm_isdst=1).png" alt="Size Limit CLI" width="738">
</p>

### Usage

---
INPUT_PATH = 
* Copy and paste data location

TIMESTAMP_FILE = 
* Absolute directory or ' ' will return index as the color bar

PLOT_LIST = 
* [  ] for default or [1, 7, 5, 3] index list for the index sequence you desire

PLOT_OFFSET = 
* Number you want to add to an offset for each curve.

PLOT_FIGURE = 
* "True" if you want to show the plots

SAVE_IMG = 
* "True" if you want to save the converted file

## License