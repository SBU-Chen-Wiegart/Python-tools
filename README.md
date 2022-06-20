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
* Download [Anaconda](https://www.anaconda.com/)
* Download athena_project.py and put it at the same folder with Larch_XAS.py

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

## CFNXRD2Jade.py
Coming soon

## plot xrd SSID insitu.py
Coming soon

## License

