import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def read_vmr_file(file_path):
    """
    Read a vertical mixing ratio file and return the data
    """
    try:
        # Read the file with numpy's genfromtxt
        vmr_data = np.genfromtxt(file_path, names=True, dtype=None, skip_header=0)
        
        # Get the column names (species)
        species_names = vmr_data.dtype.names
        
        # Get the pressure column (assuming first column is pressure)
        pressure_name = species_names[0]
        pressure = vmr_data[pressure_name]
        
        print(f"Successfully read VMR file: {file_path}")
        print(f"Found species: {', '.join(species_names[1:])}")
        print(f"Pressure column: {pressure_name}")
        
        return vmr_data, pressure, species_names
    except Exception as e:
        print(f"Error reading VMR file: {e}")
        return None, None, None

# Paths to the VMR files
vmr_file_paths = [
    "./input/couplingChem/vmr_mix_1.txt",
    "./input/couplingChem/vmr_mix_2.txt",
    "./input/couplingChem/vmr_mix_3.txt"
]

# Create a figure for plotting
plt.figure(figsize=(10, 8))

# Colors for different profiles
colors = ['b', 'r', 'g']
markers = ['o', 's', '^']

# Loop through each VMR file
for i, vmr_file_path in enumerate(vmr_file_paths):
    # Read the VMR file
    vmr_data, pressure, species_names = read_vmr_file(vmr_file_path)
    
    if vmr_data is not None:
        # Check if H2O is in the species list
        if 'H2O' in species_names:
            h2o_abundance = vmr_data['H2O']
            
            print(f"\nH2O abundance profile for {vmr_file_path}:")
            print("Pressure [bar]\tH2O VMR")
            print("-" * 30)
            
            # Print pressure and H2O abundance
            for p, h2o in zip(pressure, h2o_abundance):
                print(f"{p:.6e}\t{h2o:.6e}")
            
            # Plot H2O abundance vs pressure
            label = f'H2O (Profile {i+1})'
            plt.semilogy(h2o_abundance, pressure, color=colors[i], linestyle='-', 
                         marker=markers[i], markersize=4, markevery=5, label=label)
        else:
            print(f"H2O not found in the VMR file {vmr_file_path}.")
            print(f"Available species: {', '.join(species_names[1:])}")

# Finalize the plot
plt.gca().invert_yaxis()  # Invert y-axis to have lower pressure at top
plt.xlabel('Volume Mixing Ratio')
plt.ylabel('Pressure [bar]')
plt.title('H2O Abundance Profiles Comparison')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.xscale('log')
plt.tight_layout()
plt.savefig('h2o_abundance_profiles_comparison.png')
plt.show()
