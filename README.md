# P1414/HL-P24050 ROI
See https://hydrologicbv.github.io/P1414/intro.html for documentation.

# Creation of an anaconda environment for ROI
1. conda env create -n %NAME% --file %filepath_to_environment.yml%
2. conda activate %NAME%
3. pip install -e %filepath_to_HYDROLIB_adapted%

# Code structure
The main scripts to build a DHYDRO model using the GIS database are located under Code -> Scripts (e.g. combined_model_WBD.py)
In HYDROLIB_adapted, all the necessary functionalities are present to actually build the model. This folder is based on the Hydrolib D-HyDAMO scripting, but has some necessary adjustments needed for the ROI.
