import matplotlib.pyplot as plt
import numpy as np
import os
import sys

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
species_data = data[:, 2:]  # All species mixing ratios

# Get species numbers from header
with open('plot_scripts/temp_data.txt', 'r') as f:
    header = f.readline().strip('# \n').split()
species_numbers = [int(x) for x in header[2:]]  # Skip Temperature and Pressure

# Calculate average VMR for each species (excluding zeros and NaNs)
avg_vmr = {}
valid_species_indices = []
for i in range(species_data.shape[1]):
    vmr = species_data[:, i]
    vmr = vmr[~np.isnan(vmr)]  # Remove NaNs
    vmr = vmr[vmr > 0]        # Remove zeros
    if len(vmr) > 0:
        avg_vmr[i] = np.mean(vmr)
        valid_species_indices.append(i)

# Get top 10 species by average VMR (by column index)
top_10_species = sorted(avg_vmr.items(), key=lambda x: x[1], reverse=True)[:10]

# Create figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

# Temperature plot
ax1.plot(t_new, pressure, 'k-', linewidth=3, label='Temperature')
ax1.set_yscale('log')
ax1.set_ylim(max(pressure), min(pressure))
ax1.set_xlim(0, 1000)
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Pressure (Bar)')
ax1.set_title(f'Temperature-Pressure Profile (Step {total_step_count})')
ax1.legend()

# Abundance plot - only top 10 species
for i, (col_idx, avg) in enumerate(top_10_species):
    color = COLORS[i % len(COLORS)]
    vmr = species_data[:, col_idx]
    std_num = species_numbers[col_idx]
    name = species_names.get(std_num, f"Species {std_num}")
    label = f"{name} ({avg:.2e})"
    ax2.plot(vmr, pressure, color=color, label=label)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(max(pressure), min(pressure))
ax2.set_xlim(1e-14, 1)
ax2.set_xlabel('Volume Mixing Ratio')
ax2.set_ylabel('Pressure (Bar)')
ax2.set_title('Top 10 Species by VMR')
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig(f'{live_plot_dir}TP_profile_step_{total_step_count}.png', bbox_inches='tight')
plt.close()
