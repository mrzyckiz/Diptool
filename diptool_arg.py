import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import style
from scipy.stats import binned_statistic
from scipy.ndimage.filters import gaussian_filter1d
import subprocess
import argparse

# Create ArgumentParser object
parser = argparse.ArgumentParser(description='Perform Diptool analysis.')

# Add argument for the number of runs
parser.add_argument('-runs', type=int, help='Number of runs')

# Parse the command-line arguments
args = parser.parse_args()

# Access the value of the argument
num_runs = args.runs

# Use the number of runs 
subprocess.call(('core/Diptool_engine_v2.14.exe')) 

#definition of number of runs
for i in range(num_runs):
    subprocess.call(('core/Diptool_engine_v2.14i.exe'))



energy_file = open("energy.txt", encoding="utf8")
df = pd.read_csv(energy_file, delimiter="\s+", engine='python', skiprows=1)
df['Energy']=df['[A]']
df = df.drop(labels=['[A]', '[kcal/mol]'], axis=1)
energy_file.close()

#dataframe to array conversion
df = df.to_numpy() 

#working energy list
arr=[] 
for (z, e) in df:

    # initialization of energy plot
    if (z == 40.0):  
        arr.append([40.0, 0.0])
        e0 = e
        z0 = z
    else:
        if (z0 > z):

            # Subtract e0 for normalization in far water point
            e = e - e0  
        else:

            # Update energy for molecule reverse displacement
            e = -e + e0  
        arr.append([z, e])
        z0 = z

#list to array conversion
arra = np.array(arr) 

# number of bins to smoothen the plot
binsregulate = 40  

xarr = []
for (x, y) in arr:
    xarr.append(x)

yarr = []
for (x, y) in arr:
    yarr.append(y)

mean_stat = binned_statistic(xarr, yarr,
                             statistic='mean',
                             bins=binsregulate,
                             range=(0., 40.))

plotarray = np.insert(mean_stat.statistic, binsregulate, 0)
ysmoothed = gaussian_filter1d(plotarray, sigma=2)
pd.DataFrame(mean_stat.statistic).to_csv("energy_plot.csv")

fig1 = plt.figure()
plt.plot(mean_stat.bin_edges, ysmoothed, color='skyblue', label='E.total', linewidth=3)

# Set labels and title
plt.xlabel('Z [A]', fontsize=12)
plt.ylabel('ΔG [kcal/mol]', fontsize=12)
plt.title('Diptool energy', fontsize=14)

# Add legend and grid
plt.grid(color='black', linestyle='--', linewidth=0.5)

# Add some color to the plot background
plt.gca().set_facecolor('white')

membrane_system = pd.read_csv('lipid_coordinates.csv', delimiter=',',skiprows=1,names=['dipoles_X','dipoles_Y','dipoles_Z','total_dipole', 'X','Y','Z','lipid_type'])

# your agent trajectory from Diptool
trajectory = pd.read_csv('data.txt', delim_whitespace=True, usecols= (0,1,2), skiprows=2, names=['X','Y','Z']) 

# Data breakdown by lipid type
lipid_A = membrane_system[membrane_system['lipid_type'] == 'A']
lipid_B = membrane_system[membrane_system['lipid_type'] == 'B']
lipid_C = membrane_system[membrane_system['lipid_type'] == 'C']

# Scatter plot for Lipids A, B, and C
fig2 = plt.figure(figsize=(10, 7))
ax = fig2.add_subplot(111, projection='3d')

# Plot each lipid type if data is present
if not lipid_A.empty:
    ax.scatter(lipid_A['X'], lipid_A['Y'], lipid_A['Z'], c='salmon', label='Lipid A', marker='o', s=100, alpha=0.9, edgecolors='k')

if not lipid_B.empty:
    ax.scatter(lipid_B['X'], lipid_B['Y'], lipid_B['Z'], c='limegreen', label='Lipid B', marker='o', s=100, alpha=0.9, edgecolors='k')

if not lipid_C.empty:
    ax.scatter(lipid_C['X'], lipid_C['Y'], lipid_C['Z'], c='skyblue', label='Lipid C', marker='o', s=100, alpha=0.9, edgecolors='k')


# Plot trajectory
ax.plot(trajectory['X'], trajectory['Y'], trajectory['Z'], c='black', linestyle='--', alpha=0.5, linewidth=2,label='Agent trajectory')


# Plot surface to fill between Z coordinates
ax.set_xlabel('X [nm]')
ax.set_ylabel('Y [nm]')
ax.set_zlabel('Z [nm]')
ax.view_init(elev=10, azim=125)
ax.legend(fontsize=12)

plt.style.use('ggplot')
plt.title('Diptool trajectory')
plt.legend()
plt.show()