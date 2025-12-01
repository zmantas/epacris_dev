/**
 * cloud_opticsH.h
 * 
 * Header file for cloud optical property lookup tables from LX-Mie output
 * Reads Mie scattering tables from Clouds/LXMieOuput/ directory
 * Based on HELIOS clouds.py implementation
 * 
 * Tables structure:
 * - Multiple files: one per particle size (r0.010000.dat, r0.012589.dat, etc.)
 * - Each file: wavelength, extinction, scattering, absorption, albedo, asymmetry parameter g
 * - Files contain particle radius in filename (rXXX.XXXXXX.dat)
 */

#ifndef CLOUD_OPTICSH_H
#define CLOUD_OPTICSH_H

#include "constant.h"
// Forward declarations - actual definitions in global_temp.h
extern int CONDENSIBLES[];

// Number of particle size bins in LX-Mie lookup tables (will be determined dynamically)
// Typical range: ~0.01 to ~1000 microns, log-spaced
#define MAX_PARTICLE_SIZES 100

// Number of wavelengths in LX-Mie output files (will be determined from first file)
#define MAX_MIE_WAVELENGTHS 1000

// Maximum number of cloud species supported
#define MAX_CLOUD_SPECIES 10

// Cloud optical property data structures (LX-Mie format)
// Format: [particle_size_index][wavelength_index]
typedef struct {
    int n_particle_sizes;                                      // Actual number of particle sizes found
    int n_wavelengths;                                         // Actual number of wavelengths in files
    double radius[MAX_PARTICLE_SIZES];                        // Particle radii in micrometers (from filename)
    double wavelength[MAX_MIE_WAVELENGTHS];                   // Wavelengths in micrometers
    double extinction_cross[MAX_PARTICLE_SIZES][MAX_MIE_WAVELENGTHS];  // Extinction cross section (cm²) - total = scattering + absorption
    double albedo[MAX_PARTICLE_SIZES][MAX_MIE_WAVELENGTHS];           // Single scattering albedo
    double asymmetry_g[MAX_PARTICLE_SIZES][MAX_MIE_WAVELENGTHS];       // Asymmetry parameter g
} CloudOpticalTableMie;

// Global optical property tables (defined in cloud_optics.c)
// Array of optical tables, one per cloud species
#define MAX_CLOUD_SPECIES 10  // Maximum number of cloud species supported
extern CloudOpticalTableMie cloud_mie_optics[MAX_CLOUD_SPECIES];  // Array of optical tables for each cloud species
extern int cloud_species_ids[MAX_CLOUD_SPECIES];  // Array of species IDs corresponding to cloud_mie_optics
extern int n_cloud_species_loaded;  // Actual number of cloud species loaded

// Function declarations
void read_cloud_optical_tables_mie(void);
void cleanup_cloud_optical_tables_mie(void);

// Calculate cloud opacity arrays for radiative transfer
// Must be called after particle sizes are calculated (after cloud physics)
// Uses global arrays: clouds, particle_r2, particle_r0, particle_VP, particle_mass, wavelength
// Writes to global arrays: cH2O, aH2O, gH2O
void calculate_cloud_opacity_arrays(void);

#endif // CLOUD_OPTICSH_H

