# Python-tools
Data Analysis Tools for X-ray Spectroscopy

## Table of Contents

- [Larch XAS.py](#larch-XASpy)
- [CFNXRD2Jade.py](#CFNXRD2Jadepy)
- [plot xrd SSID insitu.py](#plot-xrd-SSID-insitupy)
- [CMS GIWAXS and GISAXS.py](#CMS-GIWAXS-and-GISAXSpy)
- [CMS SAXS data export](#CMS-SAXS-data-export)

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
* Paste the absolute directory of **larch_plot_config.ini** in **Larch_XAS.py** file 
```
CONFIG_FILE = r"D:\Research data\SSID\202205\20220509 20210221 BMM\b28_Sc_pure_config.ini"
```
* Update the config file. 

The key parameters:
```
sample_list = [0, 1]
standard_list = [1]
sample_label = ['Pure Sc']
energy_range = (4425, 4625)
output_filename = 'Sc-b33-NbAlSc-SP'
```
## CFNXRD2Jade.py
### Usage

---
* Paste your data directory
```
INPUT_PATH = r"D:\Research data\SSID\202205\20220526 XRD b31 NbAl"
```
* Set up your plotting parameters
```
FILE_TYPE = '.xy'
PLOT_LIST = [5, 4, 2]    
SAMPLE_LABEL = ['Substrate', 'Pristine', '900C60M']  
OUTPUT = False   
Y_RANGE = (-100, 500)   # Increment of ylim. ex: (ymin-100, ymax+500)
PLOT_OFFSET = 500    
PLOT_FIGURE = True  
IF_LEGEND = True    
LEGEND_LOCATION = 'upper left'
PALETTE = pld.Spectral_4_r  
CMAP = PALETTE.mpl_colormap     
OUTPUT_FILENAME = 'b31-NbAl-SiO2Si'
```
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

## CMS GIWAXS and GISAXS.py
### Usage
* Download **CMS_plot_config.ini**
* Paste your data directory in the **CMS_GIWAXS_and_GISAXS.py**
```
INPUT_PATH = r"D:\Research data\SSID\202302\20230228 CMS b33 SP\saxs\analysis\qz=0.07_dq=0.02_b33"
```
* Paste your **CMS_plot_config.ini** absolute directory
```
CONFIG_FILE = r"D:\Research data\SSID\202302\20230228 CMS b33 SP\saxs\b33-NbAlSc-SP-th0.2_CMS_plot_config.ini"
```
* Update the config file for your data plot.

The key paramteres:
```
sample_list = [0]
angle_range = 'wide' or 'small'
sample_label = ['Pristine']
output_filename = 'b33-NbAl and Sc-SP-th0.2'
output_for_jade = False or True    # Converted file for Jade reading
if_save = True or False
```

## CMS SAXS data export
* Download **SciAnalysis** folder
* Open saxs\analysis\runXS.py
* Paste the directory of the **SciAnalysis** folder
```
SciAnalysis_PATH=r"D:\Research data\SSID\Advanced Computer Python\Python-tools\SciAnalysis"
```
* Update the protocols command to reduce 2D scattering data
```
protocols = [Protocols.linecut_qr(name='qz=0.07_dq=0.02_b33', qz=0.07, dq=0.02, xlog=False, ylog=True, show_region=True, gridlines=True, plot_range=[0.004, None, 1, None])]
```
name: Exported folder name

qz: Integration center along qz direction

dq: Integration range

plot_range: [x1, y1, x2, y2]
## License