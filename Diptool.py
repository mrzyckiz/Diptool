import tkinter as tk
from tkinter import ttk
from ttkthemes import ThemedStyle
import subprocess
from tkinter import *
from PIL import ImageTk, Image
from scipy.stats import binned_statistic
from scipy.ndimage.filters import gaussian_filter1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Global variable to indicate whether the calculation should be stopped
stop_requested = False


def is_valid_number(text):
    # Convert the input text to a float, replacing commas with dots
    try:
        float(text.replace(',', '.'))
        return True
    except ValueError:
        return False

def stop_calculation():
    # Stop button function
    global stop_requested
    stop_requested = True

def calculate():
    global stop_requested

    try:
        # Placeholder function for the calculation
        num_bins = int(bins_entry.get())
        num_runs = int(runs_entry.get())  # Get the number of runs from the entry
        num_sigma = int(gaussian_entry.get()) # Get the smoothing level
        size_z_value = membrane_param_entries[2].get()


        # Check the validity of the input values
        for entry in membrane_param_entries or agent_param_entries:
            if not is_valid_number(entry.get()):
                status_label.config(text="Invalid input format.\nPlease use dot-decimal notation.", anchor="w")
                return

        if int(runs_entry.get()) < 50:
            status_label.config(text="Insufficient sampling!\n"
                                     "At least 50 runs are recommended.", anchor="w")
            root.update()  # Update the GUI to show the new label
            return

    except ValueError:
        status_label.config(text="Error. Please check your input.")
        return


    with open("input_params.txt", "w") as file:

        # Write selected membrane type
        if membrane_var.get() == "POPC":
            default_values = default_popc_values
        else:
            default_values = default_popg_values
        file.write("Membrane params: " +" ".join(map(str, default_values))+ " ")

        # Write membrane parameters
        file.write(" ".join(entry.get() for entry in membrane_param_entries)+ " "+"\n")

        # Write agent parameters
        file.write("Agent Params: " + " ".join(entry.get() for entry in agent_param_entries))

    plt.close()

    # Update the GUI to show the new label
    calculate_button.config(state='disabled')
    stop_button.config(state='normal')
    status_label.config(text="Diptool is running, please wait ...")
    root.update()


    # Run Diptool_engine with the specified number of runs
    subprocess.call(('core/Diptool_engine_1.exe'))

    for i in range(num_runs):
        #if stop_requested:
            #break
        subprocess.call(('core/Diptool_engine_2.exe'))

    '''
    if not stop_requested:
        status_label.config(text="Generating your energy files ...")
        root.update()
    '''
    # Update the GUI to show the new label
    status_label.config(text="Generating you energy files ...")
    root.update()

    # Read the energy data
    energy_file = open("energy.txt", encoding="utf8")
    df = pd.read_csv(energy_file, delimiter="\s+", engine='python', skiprows=1)
    df['Energy'] = df['[A]']
    df = df.drop(labels=['[A]', '[kcal/mol]'], axis=1)
    energy_file.close()
    df = df.to_numpy()

    # Initialize a working energy list
    arr = []
    for (z, e) in df:
        if (z == 40.0):  # initialization of energy plot
            arr.append([40.0, 0.0])
            e0 = e
            z0 = z
        else:
            if (z0 > z):
                e = e - e0  # Subtract e0 for normalization in far water point
            else:
                e = -e + e0  # Update energy for molecule reverse displacement
                # ment
            arr.append([z, e])
            z0 = z

    # Histogram preparation
    arra = np.array(arr)

    if num_bins < 20:
        num_bins = 35
    else:
        num_bins = int(num_bins)
    binsregulate = num_bins

    xarr = []
    for (x, y) in arr:
        xarr.append(x)

    yarr = []
    for (x, y) in arr:
        yarr.append(y)

    # Calculate binned statistics using scipy's binned_statistic
    mean_stat = binned_statistic(xarr, yarr,
                                    statistic='mean',
                                    bins=binsregulate,
                                    range=(0., 40.0))

    if num_sigma<2: num_sigma=2

    # Prepare data for plotting after Gaussian smoothing
    plotarray = np.insert(mean_stat.statistic, binsregulate, 0)
    ysmoothed = gaussian_filter1d(plotarray, sigma=num_sigma)
    pd.DataFrame(mean_stat.statistic).to_csv("energy_plot.csv")

    # Update the GUI to show the new label
    stop_requested = False
    calculate_button.config(state='normal')
    stop_button.config(state='disabled')
    status_label.config(text="Calculation complete.")
    root.update()

    # Plot the energy data with Matplotlib
    plt.figure(num="Diptool Energy Plot")
    plt.plot(mean_stat.bin_edges, ysmoothed, label='E_tot', linewidth=3)
    plt.ylabel('\u0394G [kcal/mol]')
    plt.title('Diptool energy')
    plt.xlabel('Z [A]')
    plt.legend()
    plt.grid(color='black', linestyle=':', linewidth=0.5)
    plt.show()



root = tk.Tk()
root.title("Diptool")
#root.geometry("800x600")  # Set the initial size of the window

# Create small_logo
icon = PhotoImage(file='pics/logo.png')
root.iconphoto(True,icon)

# Apply the ttkthemes style with "arc" theme
style = ThemedStyle(root)
style.set_theme("arc") 

# Create a label for the membrane type
membrane_label = ttk.LabelFrame(root, text="Lipid type")
membrane_label.grid(row=0, column=0, padx=20, pady=20, sticky="w")

# Create a drop-down menu for membrane type
membrane_var = tk.StringVar(value="POPC")
membrane_choices = ["POPC", "POPG"]
membrane_menu = ttk.Combobox(membrane_label, textvariable=membrane_var, values=membrane_choices)
membrane_menu.grid(row=0, column=0, padx=20, pady=20, sticky="w")

# Set default values for POPC and POPG parameters
default_popc_values = [0.34, 11.29, 0.37, 11.39, 1.65, 8.44, 18.27]
default_popg_values = [0.14, 10.17, 0.29, 10.09, 35.08, 8.29, 39.19]

# Load Diptool graphics
diptool_image = Image.open("pics/diptool_logo.png")
diptool_image = diptool_image.resize((240, 120), Image.LANCZOS)  # Resize the image
diptool_image = ImageTk.PhotoImage(diptool_image)

# Create a label for Diptool graphics
diptool_label = tk.Label(root, image=diptool_image)
diptool_label.grid(row=0, column=1, padx=20, pady=20, sticky="e")

# Create a frame for membrane parameters
membrane_frame = ttk.LabelFrame(root, text="Membrane Parameters")
membrane_frame.grid(row=1, column=0, padx=20, pady=10, sticky="nsew")

# Create input areas for membrane parameters
membrane_param_labels = ["Size X [A]", "Size Y [A]", "Size Z [A]", "Lipid APL [A²]"]
membrane_param_entries = []
for i, label_text in enumerate(membrane_param_labels):
    label = ttk.Label(membrane_frame, text=label_text)
    entry = ttk.Entry(membrane_frame)
    label.grid(row=i, column=0, padx=10, pady=5, sticky="w")
    entry.grid(row=i, column=1, padx=10, pady=5, sticky="e")
    membrane_param_entries.append(entry)

# Create a frame for agent parameters
agent_frame = ttk.LabelFrame(root, text="Agent Parameters")
agent_frame.grid(row=1, column=1, padx=20, pady=10, sticky="nsew")

# Create input areas for agent parameters
agent_param_labels = ["Mass [g/mol]", "Dipole X [D]", "dX", "Dipole Y [D]",
                      "dY", "Dipole Z [D]", "dZ", "LogP"]
agent_param_entries = []
for i, label_text in enumerate(agent_param_labels):
    label = ttk.Label(agent_frame, text=label_text)
    entry = ttk.Entry(agent_frame)
    label.grid(row=i, column=0, padx=10, pady=5, sticky="w")
    entry.grid(row=i, column=1, padx=10, pady=5, sticky="e")
    agent_param_entries.append(entry)

# Create a frame for histogram options
histogram_frame = ttk.LabelFrame(root, text="Histogram Options")
histogram_frame.grid(row=2, column=0, padx=20, pady=10, sticky="nsew")

# Create input areas for histogram options
bins_label = ttk.Label(histogram_frame, text="Bins")
bins_entry = ttk.Entry(histogram_frame)
gaussian_label = ttk.Label(histogram_frame, text="Gaussian Filter")
gaussian_entry = ttk.Entry(histogram_frame)

bins_label.grid(row=0, column=0, padx=10, pady=5, sticky="w")
bins_entry.grid(row=0, column=1, padx=10, pady=5, sticky="e")
gaussian_label.grid(row=1, column=0, padx=10, pady=5, sticky="w")
gaussian_entry.grid(row=1, column=1, padx=10, pady=5, sticky="e")

# Create a frame for Diptool runs
runs_frame = ttk.LabelFrame(root, text="Diptool Runs")
runs_frame.grid(row=2, column=1, padx=20, pady=10, sticky="nsew")

# Create an input area for number of runs
runs_label = ttk.Label(runs_frame, text="Number of Runs")
runs_entry = ttk.Entry(runs_frame)
runs_label.grid(row=0, column=0, padx=10, pady=5, sticky="w")
runs_entry.grid(row=0, column=1, padx=10, pady=5, sticky="e")


# Create a Calculate button
calculate_button = ttk.Button(root, text="Calculate", command=calculate)
calculate_button.grid(row=3, column=0, padx=20, pady=20, sticky="e", columnspan=2)

# Create a Stop button
stop_button = ttk.Button(root, text="Stop", command=stop_calculation, state='disabled')
stop_button.grid(row=3, column=1, padx=110, pady=20, sticky="e")


#Create a label to display status
status_label = tk.Label(root, text="", font=("Helvetica", 12),justify="left")
status_label.grid(row=3, column=0, padx=10, pady=10, sticky='w')

# Adjust column and row weights for resizing
root.columnconfigure(1, weight=1)
root.rowconfigure(4, weight=1)

root.mainloop()