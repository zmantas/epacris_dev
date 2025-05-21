import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Get step count and output directory from command line arguments
step_count = int(sys.argv[1])
output_dir = sys.argv[2]

# Load data
data = np.loadtxt('live_plot/temp_data.txt')
t_new = data[:, 0]
t_old = data[:, 1]
pressure = data[:, 2] / 1e5  # Convert to bar

# Create plot
plt.figure(figsize=(8, 6))
plt.semilogy(t_new, pressure, 'k-', linewidth=3, label='T$_{new}$')
plt.semilogy(t_old, pressure, 'r--', linewidth=3, label='T$_{old}$')
plt.ylim(200, 1e-8)
plt.xlim(800, 3200)
plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Bar)')
plt.title(f'Temperature-Pressure Profile (Step {step_count})')
plt.legend()

# Save figure
plt.savefig(f'{output_dir}TP_profile_step_{step_count}.png')
plt.close()
