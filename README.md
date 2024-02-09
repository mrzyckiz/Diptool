# Diptool

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10641713.svg)](https://doi.org/10.5281/zenodo.10641713)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10641713.svg)](https://doi.org/10.3390/ma14216455)

<div align="right">
  <img src="pics/diptool_logo.png" alt="logo" style="float:right; width:100px;">
</div>


 Diptool is a screening tool for a rapid determination of the antimicrobial agent affinity to various types of homogenous lipid membranes delivering particle trajectory visualization and free energy analysis. It's significantly faster than classical methods, reaching a one million-fold compared to the MD approach. 
 It provides a graphical user interface (GUI) with features for visualizing energy plots, calculating statistics, and running simulations with customizable parameters.
 
 
## Features

- **Visualize Energy Plots**: Diptool allows you to generate energy plots based on membrane and agent features. You can analyze changes in energy with respect to molecular displacement.
  
- **Customizable Parameters**: You can set various parameters for the membrane and agent properties, such as size, dipole, mass, and more.

- **Histogram Options**: Diptool enables you to adjust histogram settings or apply plot filtering, for more accurate energy analysis.

- **Simulation Runs**: The application supports running multiple simulation runs with user-defined parameters to gather comprehensive data.

- **User-Friendly GUI**: The GUI is designed with the user in mind, providing an intuitive interface for parameter input and result visualization.

<div align="center">
  <img src="pics/diptool_gui.png" alt="GUI">
</div>

## General info

Diptool combines a dedicated C++ engine and a Python visualization package (Python Software Foundation, Wilmington, DA, USA.) to offer a versatile solution for bilayer-agent analysis.
Three-dimensional trajectories and energy profile plots are generated based on delivered input from the Diptool engine. 

To initiate Diptool, simply execute the provided Python script, Diptool_run.py. This script encapsulates the tool's engine and automates its launch. Once calculations conclude, Diptool generates trajectory and energy plots.
Throughout this process, five distinct files are created:

- **input_params.txt**: Stores input membrane and agent parameters
- **membrane.txt**: Contains the arrangement of dipoles in the X, Y, and Z directions.
- **data.txt**: Encompasses the agent's trajectory along the X, Y, and Z axes.
- **energy.txt**: Presents the energy profile concerning the bilayer's normal in the Z direction.
- **energy_plot.csv**: A histogram of energy profiles further for visualization purpose.

Exemplary results can be found in the [`tests`](./tests/) folder. 


## Getting Started

1. Clone this repository to your local machine.

2. Install the required dependencies by running:

```bash
pip install -r requirements.txt
```

3. Launch Diptool by running:

```bash
python diptool_gui.py
```


### Prerequisites
 - Python 3.7
 -  Required Python packages: **matplotlib**, **numpy**, **pandas**, **pillow**, **scipy**, **tkinter**, **ttkthemes**

 
 
 ## Executable Version

An executable version of Diptool is available, which includes all the necessary libraries (all in one). 

To run the executable, download Zenodo repository content, unpack and run **Diptool.exe**. This eliminates the need to install Python, dependencies, or libraries separately, making it a convenient option for some users.

***Please note that the executable versions may take a moment to launch as they need to unpack all necessary components.***


## Usage

1. Choose the lipid type and customize membrane and agent parameters in the GUI.

2. Adjust histogram settings and the number of simulation runs.

3. Click the **Calculate** button to start simulations and energy analysis.

4. View the energy plot and analysis results generated by Diptool.


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


# The methodology and general description you can find in the article: 

 ### M. Rzycki, S. Kraszewski, M. Gladysiewicz-Kudrawiec, *Diptool a novel numerical tool for membrane interactions analysis, applying to antimicrobial detergents and drug delivery aids.* Materials 2021, 14, 6455. 
 
 ---

