# CDRMEx
Carbon Dioxide Removal (CDR) Modeling Experiments

##### CC-BY-4.0, 2019 Shannon A. Fiume

This project models highly speculative Carbon Dioxide Removal to understand
its effects and speculate how much carbon may need to be removed to return to a
carbon dioxide concentration of 280 ppm. The experiments are performed in MAGICC6.8
and have been run on pymagicc. The repo contains the scenario input files for MAGICC
and a notebook that outlines the experiments and results.

The experiments are shown in [ONCtests.ipynb](ONCtests.ipynb) which is 
a jupyter notebook that runs pymagicc, and requires windows or 
wine when run on a non-windows platform. To run these experiments, download 
[wine](https://sourceforge.net/projects/wine/files/latest/download),
[python](https://www.python.org/downloads/), pip, 
[pymagicc](https://github.com/openscm/pymagicc), this repo, and open the 
notebook in jupyter.

#### Install and run the workbook
Download/install [wine](https://sourceforge.net/projects/wine/files/latest/download)

Next open a terminal, and add wine to the path.

Then run:
```
pip install -r requirements.txt
jupyter-notebook ONCtests.ipynb
```

#### Install for development 
Open a terminal and do something like the following:

```
which wine
git clone https://github.com/hsbay/cdrmex
git clone https://github.com/openscm/pymagicc
cd pymagicc
make venv
./venv/bin/pip install --editable .
./venv/bin/pip install ipywidgets appmode
./venv/bin/pip install -r requirements.txt
jupyter nbextension enable --py --sys-prefix widgetsnbextension
jupyter nbextension     enable --py --sys-prefix appmode
jupyter serverextension enable --py --sys-prefix appmode
./venv/bin/jupyter-notebook ../cdrmex/ONCtests.ipynb
```

After the notebook is up, run all the cells, if they haven't already been populated.

This workbook uses [pymagicc](https://pymagicc.readthedocs.io/en/latest/) by R. Gieseke, S. N. Willner and M. Mengel, (2018). 
Pymagicc: A Python wrapper for the simple climate model MAGICC. 
   Journal of Open Source Software, 3(22), 516, 
   https://doi.org/10.21105/joss.00516

[MAGICC](http://magicc.org/) is by:
    M. Meinshausen, S. C. B. Raper and T. M. L. Wigley (2011). 
    “Emulating coupled atmosphere-ocean and carbon cycle models with a simpler model, MAGICC6: Part I “Model Description and Calibration.” 
    Atmospheric Chemistry and Physics 11: 1417-1456. 
    https://doi.org/10.5194/acp-11-1417-2011

This software is CC-BY-4.0 and carries no warranty towards any liability, use at your own risk.
See license.txt for more information.
