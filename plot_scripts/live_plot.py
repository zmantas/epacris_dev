import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd

# Get step count and output directory from command line arguments
total_step_count = int(sys.argv[1])
live_plot_dir = sys.argv[2]

# Path to species mapping file (adjust as needed)
species_file = 'Condition/SpeciesList/species_helios_comp3.dat'
species_names = {}
if os.path.exists(species_file):
    with open(species_file, 'r') as f:
        f.readline()  # skip header
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    std_num = int(parts[2])
                    name = parts[0]
                    species_names[std_num] = name
                except Exception:
                    continue

# Define colors for cycling through - using 10 different colors
COLORS = ['k', 'r', 'g', 'b', 'm', 'c', 'y', 'orange', 'purple', 'brown']

# Load data
data = np.loadtxt('plot_scripts/temp_data.txt', skiprows=1)  # Skip header row
t_new = data[:, 0]
pressure = data[:, 1]  # Pressure in bar

# Get species numbers from header
with open('plot_scripts/temp_data.txt', 'r') as f:
    header = f.readline().strip('# \n').split()

# Separate gas species and cloud species
gas_indices = []
cloud_indices = []
gas_species_numbers = []
cloud_species_numbers = []

for i, col_name in enumerate(header[2:], 2):  # Skip Temperature and Pressure
    if col_name.startswith('c'):
        # This is a cloud species (with 'c' prefix)
        cloud_indices.append(i)
        cloud_species_numbers.append(int(col_name[1:]))  # Remove 'c' prefix to get species number
    else:
        # This is a gas species
        gas_indices.append(i)
        gas_species_numbers.append(int(col_name))

# Extract gas and cloud data
gas_data = data[:, gas_indices]
cloud_data = data[:, cloud_indices]

# Calculate average VMR for each gas species (excluding zeros and NaNs)
gas_avg_vmr = {}
for i, col_idx in enumerate(gas_indices):
    vmr = data[:, col_idx]
    vmr = vmr[~np.isnan(vmr)]  # Remove NaNs
    vmr = vmr[vmr > 0]         # Remove zeros
    if len(vmr) > 0:
        gas_avg_vmr[i] = np.mean(vmr)

# Get top 10 gas species by average VMR
top_10_gas = sorted(gas_avg_vmr.items(), key=lambda x: x[1], reverse=True)[:10]

# Create figure with two subplots
fig = plt.figure(figsize=(20, 10))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])  # Equal width panels
ax1 = fig.add_subplot(gs[0])  # Temperature plot
ax2 = fig.add_subplot(gs[1])  # Abundance plot

# Temperature plot
ax1.plot(t_new, pressure, 'k-', linewidth=3, label='Temperature')
ax1.set_yscale('log')
ax1.set_ylim(max(pressure), min(pressure))
ax1.set_xlim(0,1500)
ax1.set_xlabel('Temperature (K)', fontsize=12)
ax1.set_ylabel('Pressure (Bar)', fontsize=12)
ax1.set_title(f'Temperature Profile (Step {total_step_count})', fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Create a dictionary to map species to colors for consistent coloring
species_colors = {}

# Abundance plot - top 10 gas species
for i, (idx, avg) in enumerate(top_10_gas):
    color = COLORS[i % len(COLORS)]
    vmr = gas_data[:, idx]
    std_num = gas_species_numbers[idx]
    name = species_names.get(std_num, f"Species {std_num}")
    label = f"{name} ({avg:.2e})"
    ax2.plot(vmr, pressure, color=color, linewidth=2, label=label)
    
    # Store color for this species
    species_colors[std_num] = color

# Check for and plot non-zero cloud data
has_clouds = False
for i, std_num in enumerate(cloud_species_numbers):
    cloud_vmr = cloud_data[:, i]
    
    # Check if cloud has any significant values (above threshold)
    if np.any(cloud_vmr > 1e-10):
        has_clouds = True
        # Use same color as gas species if available
        color = species_colors.get(std_num, COLORS[i % len(COLORS)])
        name = species_names.get(std_num, f"Species {std_num}")
        
        # Get non-zero values for proper averaging
        significant_values = cloud_vmr[cloud_vmr > 1e-10]
        avg_vmr = np.mean(significant_values)
        
        # Plot clouds with dashed lines and different marker
        ax2.plot(
            cloud_vmr, 
            pressure, 
            color=color, 
            linestyle='--',
            marker='o', 
            markersize=4,
            linewidth=2,
            label=f"{name} (cloud, {avg_vmr:.2e})"
        )
        print(f"Plotting cloud data for species {std_num} ({name})")

if not has_clouds:
    print("WARNING: No significant cloud abundances found to plot!")

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(max(pressure), min(pressure))
ax2.set_xlim(1e-14, 1)
ax2.set_xlabel('Volume Mixing Ratio', fontsize=12)
ax2.set_ylabel('Pressure (Bar)', fontsize=12)
ax2.set_title('Species Abundance', fontsize=14)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{live_plot_dir}TP_profile_step_{total_step_count}.png', bbox_inches='tight')
plt.close()
