# FuelLib
FuelLib (SWR-25-26) utilizes the tables and functions of the Group Contribution Method (GCM) as proposed by [Constantinou and Gani (1994)](https://doi.org/10.1002/aic.690401011) and [Constantinou, Gani and O'Connel (1995)](https://doi.org/10.1016/0378-3812(94)02593-P), with additional physical properties discussed in [Govindaraju & Ihme (2016)](https://doi.org/10.1016/j.ijheatmasstransfer.2016.06.079).  The code is based on Pavan B. Govindaraju's [Matlab implementation](https://github.com/gpavanb-old/GroupContribution) of the GCM, and has been expanded to include additional thermodynamic properties and mixture properties.  The fuel library contains gas chromatography (GC x GC) data for a variety of fuels ranging from simple single component fuels to complex jet fuels.  The GC x GC data for POSF jet fuels comes from [Edwards (2020)](https://apps.dtic.mil/sti/pdfs/AD1093317.pdf).  

## Python Environment
The following conda environment is required to run this code:
~~~
conda create --name fuellib-env matplotlib pandas black 
~~~

## Running the code
This repository includes an example `ex-mixtureProperties.py` that calculates a given mixture's density, viscosity and vapor pressure from GC x GC data of the mixture.  The results are plotted against data from NIST and [Edwards (2020)](https://apps.dtic.mil/sti/pdfs/AD1093317.pdf).

# Contributing
New contributions are always welcome.  If you have an idea for a new feature follow these steps:
1. Fork the main repository
2. Create a `newFeature` branch that contains your changes
3. Update the sphinx documentation in `newFeature`
4. Format the source code files using the [Black code formatter](https://github.com/psf/black) by running the following command:
   ~~~
   find . -name "*.py" -print0 | xargs -0 black
   ~~~
5. Open a Pull Request (PR) from `newFeature` on your fork to branch `main` FuelLib repository.

## Sphinx Documentation
This repository uses [Sphinx](https://www.sphinx-doc.org/en/master/usage/quickstart.html) to generate documentation.  This requires the following Conda environment:
~~~
conda create --name sphinx-env
conda activate sphinx-env
conda install python
conda install anaconda::sphinx
conda install conda-forge::sphinx_rtd_theme 
conda install conda-forge::sphinxcontrib-bibtex
~~~

To view the documentation locally, build the html using the following: 
~~~
cd FuelLib/docs/
sphinx-build -M html . _build/
~~~
You should now be able to view the html by opening `FuelLib/docs/_build/html/index.html` in a web browser. 
