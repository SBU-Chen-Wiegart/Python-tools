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

If merging fluorescence scans (_e.g._ BNL NSLS-II BMM)
```
FILE_TYPE = '.prj'
```
If merging transmission scans (_e.g._ BNL NSLS-II BMM)
```
FILE_TYPE = ''
```
If merging transmission scans (_e.g._ BNL NSLS-II ISS)
```
FILE_TYPE = '.dat'
```
If plotting scans file exported from Athena
```
FILE_TYPE = '.txt'
```

Copy and paste the directory of the data folder
```
INPUT_PATH = r'D:\Research data\SSID\202209\20220330 ISS NbOx'
```

#### Data display (optional)

Add scan name you want to exclude
```
SKIP_SCANS = ['MnO2_45_16C_Charge_Mn_001']
```

Do normalization for fluorescence scans
```
IF_NOR = False
```

Add plus and minus standard deviation lines for fluorescence scans
```
ADD_DEV = False
```

Display Athena parameters, such as atomic symbol, edge, label, etc.
```
SHOW_DATA_INFORMATION = False
```

#### Plot spectrum (optional)

* Download **larch_plot_config.ini**
* Copy and paste the absolute directory of the config file for plotting .txt files
```
CONFIG_FILE = r"D:\Research data\SSID\202205\20220509 20210221 BMM\b28_Sc_pure_config.ini"
```
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
Copy and paste data location
```
INPUT_PATH = r'D:\Research data\SSID\202204\20220406 XPD ex-situ check\LK_b30-14_Nb40Al60Sc_SiO2Si_pristine_heating'
```
Absolute directory or ' ' will return index as the color bar
```
TIMESTAMP_FILE = 'sample_LK_heating_20220408_172022.xlsx'
```
[  ] for default or [1, 7, 5, 3] index list for the index sequence you desire
```
PLOT_LIST = list(np.arange(0, 191, 5)) 
```
Number you want to add to an offset for each curve
```
PLOT_OFFSET = 0.5
```

"True" if you want to show the plots
```
PLOT_FIGURE = True
```
"True" if you want to save the converted file
```
SAVE_IMG = False
```
## License