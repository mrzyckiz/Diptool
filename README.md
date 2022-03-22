# Diptool
 Diptool is a screening tool for a rapid determination of the Gemini agent affinity to various types of homogenous lipid membranes delivering particle trajectory visualization and free energy analysis. It's significantly faster than classical methods, reaching a one million-fold compared to the MD approach. 
 
 
 

## Requirements
 -Python 3.7 or higher
 
 -matplotlib, math, pyplot, numpy, seaborn, scipy libraries 
 
 
## General info
Diptool is provided with a dedicated engine written in C++ and visualization package implemented in Python (Python Software Foundation, Wilmington, DA, USA.) and requires version 3.7 or higher. Three-dimensional trajectories and energy profile plots are generated based on delivered input from the Diptool engine. 

Diptool should be executed from an included python script (Diptool_run.py), in which the tool engine is embedded and automatically launched, and after finished calculations a trajectory and energy plots are generated. In the meantime, three files with plotted membrane position, agent trajectory, and energy profiles are produced, respectively. In membrane.txt file the dipoles arrangement in X, Y, Z direction are stored, data.txt include the agent trajectory in the X, Y, Z axes, and finally energy.txt contain the energy profile with respect to the bilayer normal—Z direction. The python visualization code was commented including several hints to simplify the usage.


# The methodology and general description you can find in the article: 

 ### M. Rzycki, S. Kraszewski, M. Gładysiewicz-Kudrawiec, *Diptool – a novel numerical tool for membrane interactions analysis, applying to antimicrobial detergents and drug delivery aids.* Materials 2021, 14, 6455. 

