import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import style
from scipy.stats import binned_statistic
from scipy.ndimage.filters import gaussian_filter1d
from pandas import read_csv
import subprocess

subprocess.call(('Diptool_engine_extra_err.exe'))

membrane_system = np.loadtxt('membrane.txt', dtype= np.float64,delimiter='  ',usecols=(0,1,2), unpack=True, skiprows=2) # Membrane dipoles distribution
trajectory = np.loadtxt('data.txt', dtype= np.float64, delimiter='  ', usecols= (0,1,2), unpack=True, skiprows =2) # your agent trajectory from Diptool



plt.rcParams.update({'font.size': 20})
#print(plt.style.available)
style.use('seaborn-poster')
sns.set_palette(sns.color_palette("husl", 5))

#----------------------------------------------------------------------------------- Figure 1 "Particle trajectory" -----------------------------------------------------------------------------------
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')

ax.plot(membrane_system[0]*10, membrane_system[1]*10, membrane_system[2]*10, 'o', markersize='15', label ='Membrane Dipoles', color='teal') # visualization details of membrane
ax.plot(trajectory[0]*10, trajectory[1]*10, trajectory[2]*10, linewidth=3, label='Agent trajectory', linestyle='dashed') # visualization details of agent trajectory



ax.axis('auto')
ax.set_xlabel('X [A]')
ax.set_ylabel('Y [A]')
ax.set_zlabel('Z [A]')
ax.yaxis.labelpad=30
ax.zaxis.labelpad=20
ax.xaxis.labelpad=30
plt.title('Diptool trajectory')
plt.legend()


#----------------------------------------------------------------------------------- Figure 2 "Particle energy vs time" -----------------------------------------------------------------------------------
#fig1, (ax1) = plt.subplots(1, sharex=True)

energy_file = open("energy.txt", encoding="utf8") # your energy file from Diptool
df = read_csv(energy_file, delimiter="\s+", engine='python', skiprows=1)
df['Energy']=df['[A]']
df = df.drop(labels=['[A]', '[kcal/mol]'], axis=1)
energy_file.close()
df = df.to_numpy() #dataframe to array conversion

arr=[] #working energy list
for (z, e) in df:
    if(z==45.0): #initialization of energy plot
        arr.append([45.0, 0.0])
        e0=e
        z0=z
    else:
        if(z0>z): e=e-e0 #e0 subsrtacion for normalization in far woater point
        else: e=-e+e0 #energy sign change in case of molecule reverse displacement
        arr.append([z, e])
        z0=z

arra = np.array(arr) #list to array conversion


binsregulate = 45  # number of bins to smoothen the plot

xarr = []
for (x, y) in arr:
    xarr.append(x)

yarr = []
for (x, y) in arr:
    yarr.append(y)

mean_stat = binned_statistic(xarr, yarr,
                             statistic='mean',
                             bins=binsregulate,
                             range=(0., 45.))

plotarray = np.insert(mean_stat.statistic, binsregulate, 0)
ysmoothed = gaussian_filter1d(plotarray, sigma=2)

fig2 = plt.figure(2)
plt.plot(mean_stat.bin_edges, ysmoothed, label= 'E.total', linewidth=3)
plt.ylabel('\u0394G [kcal/mol]')
plt.title('Diptool energy')
plt.xlabel('Z [A]')
plt.legend()
plt.grid(color='black', linestyle= ':', linewidth=0.5)
plt.show()






# a.create_animation()

