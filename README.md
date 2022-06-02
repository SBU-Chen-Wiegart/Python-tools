# Python-tools
Data Analysis Tools for X-ray Spectroscopy

## Table of Contents

- [Larch_XAS.py](#larch-XASpy)
- [CFNXRD2Jade.py](#cFNXRD2Jadepy)
- [plot_xrd_SSID_insitu.py](#plot-xrd-SSID-insitupy)

## Larch_XAS.py
- Documentation: http://xraypy.github.io/xraylarch
- Code: http://github.com/xraypy/xraylarch
### [Install](https://xraypy.github.io/xraylarch/installation.html)
* Download [Anaconda](https://www.anaconda.com/)
* Download athena_project.py

Within a shell:

1. activate your conda environment (called base by default) and update it:
```
conda activate
conda update -y conda python pip
```
2. (optional/expert) create a dedicated environment for Larch and activate it:
```
conda create -y --name xraylarch python==3.9
conda activate xraylarch
conda update --all
```
3. install main dependencies:
```
conda install -y "numpy=>1.20" "scipy=>1.6" "matplotlib=>3.0" scikit-learn pandas
conda install -y -c conda-forge wxpython pymatgen tomopy pycifrw
```
4. install Larch (latest release):
```
pip install xraylarch
```
5. Set Python Interpreters as Anaconda\envs\xraylarch\python.exe

## CFNXRD2Jade.py

## plot_xrd_SSID_insitu.py

## License

MIT © Richard McRichface
