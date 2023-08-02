import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import style
from scipy.stats import binned_statistic
from scipy.ndimage.filters import gaussian_filter1d
import pandas as pd
from pandas import read_csv
import subprocess
subprocess.call(('Diptool_engine_extra_err.exe'))  #definition of number of runs
for i in range(1000):
    subprocess.call(('Diptool_engine_extra_err_append.exe'))

membrane_system = np.loadtxt('membrane.txt', dtype= np.float64,delimiter='  ',usecols=(0,1,2), unpack=True, skiprows=2) # Membrane dipoles distribution
trajectory = np.loadtxt('data.txt', dtype= np.float64, delimiter='  ', usecols= (0,1,2), unpack=True, skiprows =2) # your agent trajectory from Diptool
#trajectory_ref = np.loadtxt('dane_chx.txt', dtype= np.float64, delimiter='  ', usecols= (0,1,2), unpack=True) # if you have additional agent trajectory to visualize please adjust that line
#energy_abf = np.loadtxt('chx_pc.txt', dtype=np.float64, delimiter=' ', usecols= (0,2), unpack=True) # your reference trajectory from MD
#energy_total= np.loadtxt('energy.txt', dtype=np.float64, delimiter='  ', usecols= (0,1), unpack=True, skiprows=2) # your energy file from Diptool


plt.rcParams.update({'font.size': 20})
#print(plt.style.available)
style.use('seaborn-poster')
sns.set_palette(sns.color_palette("husl", 5))

#----------------------------------------------------------------------------------- Figure 1 "Particle trajectory" -----------------------------------------------------------------------------------
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')

ax.plot(membrane_system[0]*10, membrane_system[1]*10, membrane_system[2]*10, 'o', markersize='15', label ='Membrane Dipoles', color='teal') # visualization details of membrane
ax.plot(trajectory[0]*10, trajectory[1]*10, trajectory[2]*10, linewidth=3, label='Agent trajectory', linestyle='dashed') # visualization details of agent trajectory
#ax.plot(trajectory_ref[0]*10, trajectory_ref[1]*10, trajectory_ref[2]*10, linestyle='dotted', color ='black', linewidth=3, label=' trajectory') # if you have additional agent to visualize please adjust that line


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



energy_file = open("energy.txt", encoding="utf8")
df = read_csv(energy_file, delimiter="\s+", engine='python', skiprows=1)
df['Energy']=df['[A]']
df = df.drop(labels=['[A]', '[kcal/mol]'], axis=1)
energy_file.close()
df = df.to_numpy() #dataframe to array conversion

arr=[] #working energy list
for (z, e) in df:
    if(z==35.0): #initialization of energy plot
        arr.append([35.0, 0.0])
        e0=e
        z0=z
    else:
        if(z0>z): e=e-e0 #e0 subsrtacion for normalization in far water point
        else: e=-e+e0 #energy update in case of molecule reverse displace
        # ment
        arr.append([z, e])
        z0=z

arra = np.array(arr) #list to array conversion


binsregulate = 35  # number of bins to smoothen the plot

xarr = []
for (x, y) in arr:
    xarr.append(x)

yarr = []
for (x, y) in arr:
    yarr.append(y)

mean_stat = binned_statistic(xarr, yarr,
                             statistic='mean',
                             bins=binsregulate,
                             range=(0., 35.))

plotarray = np.insert(mean_stat.statistic, binsregulate, 0)
ysmoothed = gaussian_filter1d(plotarray, sigma=2)
pd.DataFrame(mean_stat.statistic).to_csv("energy_plot.csv")

fig2 = plt.figure(2)
plt.plot(mean_stat.bin_edges, ysmoothed, label= 'E.total', linewidth=3)
plt.ylabel('\u0394G [kcal/mol]')
plt.title('Diptool energy')
plt.xlabel('Z [A]')
plt.legend()
plt.grid(color='black', linestyle= ':', linewidth=0.5)
plt.show()




#ax1.plot(energy_total[0], energy_total[1], linewidth=3,label= 'E.total' ) # Diptool free energy plot
#ax1.set_ylabel('\u0394G [kcal/mol]')
#ax1.set_title('Diptool energy')
#plt.rcParams.update({'font.size': 22})
#ax1.legend()
#ax1.set_xlabel('z [A]')

#ax2.plot(energy_abf[0], energy_abf[1], color='steelblue',linewidth=3,label= 'ABF' ) # MD free energy plot
#ax2.set_ylabel('\u0394G [kcal/mol]')
#ax2.set_title('MD energy')
#ax2.set_xlabel('Z-coordinate [A]')
#plt.rcParams.update({'font.size': 20})



#plt.show()



