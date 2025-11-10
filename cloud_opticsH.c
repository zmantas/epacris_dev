/**
 * cloud_opticsH.c
 * 
 * Functions to read and store cloud optical property lookup tables from LX-Mie output
 * Based on HELIOS clouds.py implementation
 * 
 * Data files from Clouds/LXMieOuput/H2O_s/:
 * - r0.010000.dat, r0.012589.dat, ... (particle radius in filename, in microns)
 * - Each file contains: wavelength, extinction, scattering, absorption, albedo, asymmetry g
 * 
 * File format:
 * # Header line
 * wavelength(μm)  size_param  extinction(cm²)  scattering(cm²)  absorption(cm²)  albedo  asymmetry_g
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include "cloud_opticsH.h"
#include "global_temp.h"  // For zbin, NSP, NLAMBDA, NCONDENSIBLES, CONDENSIBLES

// Global optical property table
CloudOpticalTableMie H2O_mie_optics;

/**
 * Extract particle radius from filename like "r0.010000.dat"
 * Returns radius in micrometers, or -1.0 on error
 */
static double extract_radius_from_filename(const char *filename) {
    // Skip "r" prefix
    if (filename[0] != 'r') {
        return -1.0;
    }
    
    // Find ".dat" suffix
    char *dot_dat = strstr(filename, ".dat");
    if (dot_dat == NULL) {
        return -1.0;
    }
    
    // Extract number between "r" and ".dat"
    char num_str[64];
    int len = dot_dat - (filename + 1);
    if (len <= 0 || len >= 64) {
        return -1.0;
    }
    
    strncpy(num_str, filename + 1, len);
    num_str[len] = '\0';
    
    return atof(num_str);
}

/**
 * Compare function for sorting filenames by particle radius
 */
static int compare_particle_files(const void *a, const void *b) {
    const char *file_a = *(const char **)a;
    const char *file_b = *(const char **)b;
    
    double radius_a = extract_radius_from_filename(file_a);
    double radius_b = extract_radius_from_filename(file_b);
    
    if (radius_a < radius_b) return -1;
    if (radius_a > radius_b) return 1;
    return 0;
}

/**
 * Read a single LX-Mie output file
 * 
 * File format:
 * # Header line (skip)
 * wavelength(μm)  size_param  extinction(cm²)  scattering(cm²)  absorption(cm²)  albedo  asymmetry_g
 */
static int read_mie_file(const char *filepath, double *wavelengths, 
                         double *extinction, double *albedo, double *asymmetry,
                         int *n_wavelengths) {
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) {
        printf("ERROR: Cannot open LX-Mie file: %s\n", filepath);
        return -1;
    }
    
    // Skip header line (starts with #)
    char header[1024];
    if (fgets(header, sizeof(header), fp) == NULL) {
        printf("ERROR: Failed to read header from %s\n", filepath);
        fclose(fp);
        return -1;
    }
    
    // Read data lines
    // Format: wavelength(μm)  size_param  extinction(cm²)  scattering(cm²)  absorption(cm²)  albedo  asymmetry_g
    int count = 0;
    double wl, size_param, ext, scat, abs, alb, asym_val;
    
    while (fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", 
                 &wl, &size_param, &ext, &scat, &abs, &alb, &asym_val) == 7) {
        if (count >= MAX_MIE_WAVELENGTHS) {
            printf("WARNING: Too many wavelengths in %s, truncating at %d\n", 
                   filepath, MAX_MIE_WAVELENGTHS);
            break;
        }
        
        wavelengths[count] = wl;           // Column 0: wavelength in microns
        extinction[count] = ext;           // Column 2: extinction cross section (cm²) - total = scattering + absorption
        albedo[count] = alb;               // Column 5: single scattering albedo
        asymmetry[count] = asym_val;       // Column 6: asymmetry parameter g
        
        count++;
    }
    
    fclose(fp);
    
    if (count == 0) {
        printf("ERROR: No data read from %s\n", filepath);
        return -1;
    }
    
    *n_wavelengths = count;
    return 0;
}

/**
 * Read all LX-Mie output files for H2O
 * Scans directory, finds all r*.dat files, sorts by particle radius, reads them
 */
void read_cloud_optical_tables_mie(void) {
    printf("\n");
    printf("=======================================================\n");
    printf("READING CLOUD OPTICAL PROPERTY TABLES (LX-MIE FORMAT)\n");
    printf("=======================================================\n");
    
    const char *directory = "Clouds/LXMieOuput/H2O_s";
    
    // Initialize structure
    H2O_mie_optics.n_particle_sizes = 0;
    H2O_mie_optics.n_wavelengths = 0;
    
    // Open directory and collect all r*.dat filenames
    DIR *dir = opendir(directory);
    if (dir == NULL) {
        printf("FATAL ERROR: Cannot open directory: %s\n", directory);
        exit(1);
    }
    
    // Collect all r*.dat filenames
    char **filenames = (char **)malloc(MAX_PARTICLE_SIZES * sizeof(char *));
    int file_count = 0;
    
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL && file_count < MAX_PARTICLE_SIZES) {
        const char *name = entry->d_name;
        
        // Check if filename matches pattern r*.dat
        if (name[0] == 'r' && strstr(name, ".dat") != NULL) {
            // Check if radius extraction works
            double radius = extract_radius_from_filename(name);
            if (radius > 0.0) {
                filenames[file_count] = (char *)malloc(strlen(name) + 1);
                strcpy(filenames[file_count], name);
                file_count++;
            }
        }
    }
    closedir(dir);
    
    if (file_count == 0) {
        printf("FATAL ERROR: No valid r*.dat files found in %s\n", directory);
        exit(1);
    }
    
    // Sort filenames by particle radius (ascending)
    qsort(filenames, file_count, sizeof(char *), compare_particle_files);
    
    printf("Found %d particle size files\n", file_count);
    
    // Read first file to determine number of wavelengths
    char first_filepath[512];
    sprintf(first_filepath, "%s/%s", directory, filenames[0]);
    
    //Store temp var for reading out
    double temp_wavelengths[MAX_MIE_WAVELENGTHS];
    double temp_extinction[MAX_MIE_WAVELENGTHS];
    double temp_albedo[MAX_MIE_WAVELENGTHS];
    double temp_asymmetry[MAX_MIE_WAVELENGTHS];
    
    if (read_mie_file(first_filepath, temp_wavelengths, temp_extinction, 
                     temp_albedo, temp_asymmetry, &H2O_mie_optics.n_wavelengths) != 0) {
        printf("FATAL ERROR: Failed to read first file: %s\n", first_filepath);
        exit(1);
    }
    
    // Copy wavelengths from first file (should be same for all particle sizes)
    for (int j = 0; j < H2O_mie_optics.n_wavelengths; j++) {
        H2O_mie_optics.wavelength[j] = temp_wavelengths[j];
    }
    
    printf("Number of wavelengths per file: %d\n", H2O_mie_optics.n_wavelengths);
    
    // Read all particle size files
    for (int i = 0; i < file_count; i++) {
        char filepath[512];
        sprintf(filepath, "%s/%s", directory, filenames[i]);
        
        double radius = extract_radius_from_filename(filenames[i]);
        H2O_mie_optics.radius[i] = radius;
        
        if (read_mie_file(filepath, temp_wavelengths, temp_extinction,
                         temp_albedo, temp_asymmetry, &H2O_mie_optics.n_wavelengths) != 0) {
            printf("WARNING: Failed to read %s, skipping\n", filepath);
            continue;
        }
        
        // Verify wavelength grid matches
        int wavelengths_match = 1;
        for (int j = 0; j < H2O_mie_optics.n_wavelengths; j++) {
            if (fabs(temp_wavelengths[j] - H2O_mie_optics.wavelength[j]) > 1e-10) {
                wavelengths_match = 0;
                break;
            }
        }
        
        if (!wavelengths_match) {
            printf("WARNING: Wavelength grid mismatch in %s, skipping\n", filepath);
            continue;
        }
        
        // Copy optical properties, i is the particle size index, j is the wavelength index
        for (int j = 0; j < H2O_mie_optics.n_wavelengths; j++) {
            H2O_mie_optics.extinction_cross[i][j] = temp_extinction[j];
            H2O_mie_optics.albedo[i][j] = temp_albedo[j];
            H2O_mie_optics.asymmetry_g[i][j] = temp_asymmetry[j];
        }
        
        H2O_mie_optics.n_particle_sizes++;
    }
    
    // Print debug in scientific notation
    printf("n_particle_sizes = %d\n", H2O_mie_optics.n_particle_sizes);
    printf("n_wavelengths = %d\n", H2O_mie_optics.n_wavelengths);
    printf("radius[0] = %e\n", H2O_mie_optics.radius[0]);
    printf("wavelength[0] = %e\n", H2O_mie_optics.wavelength[0]);
    printf("extinction_cross[0][0] = %e\n", H2O_mie_optics.extinction_cross[0][0]);
    printf("albedo[0][0] = %e\n", H2O_mie_optics.albedo[0][0]);
    printf("asymmetry_g[0][0] = %e\n", H2O_mie_optics.asymmetry_g[0][0]);
    printf("extinction_cross[0][1] = %e\n", H2O_mie_optics.extinction_cross[0][1]);
    printf("albedo[0][1] = %e\n", H2O_mie_optics.albedo[0][1]);
    printf("asymmetry_g[0][1] = %e\n", H2O_mie_optics.asymmetry_g[0][1]);
    printf("extinction_cross[0][2] = %e\n", H2O_mie_optics.extinction_cross[0][2]);



    // Free filename memory
    for (int i = 0; i < file_count; i++) {
        free(filenames[i]);
    }
    free(filenames);
    
    printf("\n");
    printf("=======================================================\n");
    printf("CLOUD OPTICAL TABLES (LX-MIE) SUCCESSFULLY LOADED\n");
    printf("=======================================================\n");
    printf("Particle size range: %.6f - %.6f micrometers\n", 
           H2O_mie_optics.radius[0], 
           H2O_mie_optics.radius[H2O_mie_optics.n_particle_sizes-1]);
    printf("Number of particle sizes: %d\n", H2O_mie_optics.n_particle_sizes);
    printf("Number of wavelengths: %d\n", H2O_mie_optics.n_wavelengths);
    printf("Wavelength range: %.6f - %.6f micrometers\n",
           H2O_mie_optics.wavelength[0],
           H2O_mie_optics.wavelength[H2O_mie_optics.n_wavelengths-1]);
    printf("Species loaded: H2O (solid/ice)\n");
    printf("=======================================================\n\n");
}

/**
 * Cleanup function (if needed for future memory management)
 */
void cleanup_cloud_optical_tables_mie(void) {
    // No dynamic memory to free currently
    // This function is provided for future extension
}

/**
 * Linear interpolation in log-space (for particle size interpolation)
 * Interpolates y (linear) vs log10(x) - used for albedo and asymmetry
 * Similar to MATLAB interp1(log10(x), y, log10(xi))
 */
static double interp1_log(double *x, double *y, int n, double xi) {
    if (n < 2) return y[0];
    
    double log_xi = log10(fmax(x[0], fmin(xi, x[n-1])));
    
    // Handle boundaries
    if (log_xi <= log10(x[0])) return y[0];
    if (log_xi >= log10(x[n-1])) return y[n-1];
    
    // Find interpolation interval
    for (int i = 0; i < n-1; i++) {
        double log_x0 = log10(x[i]);
        double log_x1 = log10(x[i+1]);
        
        if (log_xi >= log_x0 && log_xi <= log_x1) {
            // Linear interpolation in log space
            double t = (log_xi - log_x0) / (log_x1 - log_x0);
            return y[i] + t * (y[i+1] - y[i]);
        }
    }
    
    return y[n-1]; // fallback
}

/**
 * Log-log interpolation (for cross-section interpolation)
 * Interpolates log10(y) vs log10(x), then converts back
 * Similar to MATLAB pow(10, interp1(log10(x), log10(y), log10(xi)))
 * Matches cloud_interp.c line 167: pow(10, interp1(log10(H2OI_r), log10(H2OI_c), log10(...)))
 */
static double interp1_loglog(double *x, double *y, int n, double xi) {
    if (n < 2) return y[0]; // check if fewer than 2 points, no interpolation needed
    
    // convert r0 to log10 space while maintining boundaries
    double log_xi = log10(fmax(x[0], fmin(xi, x[n-1])));
    
    // If boundary conditions are met, return the value at the boundary
    if (log_xi <= log10(x[0])) return y[0];
    if (log_xi >= log10(x[n-1])) return y[n-1];
    
    // Find interpolation interval
    // Loop through all particle sizes
    for (int i = 0; i < n-1; i++) {
        double log_x0 = log10(x[i]);
        double log_x1 = log10(x[i+1]);
        
        // If between boundary conditions, perform linear interpolation in log-log space
        if (log_xi >= log_x0 && log_xi <= log_x1) {
            // Linear interpolation in log-log space
            double log_y0 = log10(fmax(1e-30, y[i]));      // Protect against zero/negative
            double log_y1 = log10(fmax(1e-30, y[i+1]));
            double t = (log_xi - log_x0) / (log_x1 - log_x0); // calculate the weight factor
            double log_y_interp = log_y0 + t * (log_y1 - log_y0);
            return pow(10.0, log_y_interp);
        }
    }
    
    return y[n-1]; // fallback
}

/**
 * Linear interpolation in wavelength (for wavelength grid mapping)
 */
static double interp1_wavelength(double *wavelength_mie, double *value_mie, int n_mie,
                                double wavelength_rt) {
    if (n_mie < 2) return value_mie[0];
    
    // Handle boundaries
    if (wavelength_rt <= wavelength_mie[0]) return value_mie[0];
    if (wavelength_rt >= wavelength_mie[n_mie-1]) return value_mie[n_mie-1];
    
    // Find interpolation interval
    for (int i = 0; i < n_mie-1; i++) {
        if (wavelength_rt >= wavelength_mie[i] && wavelength_rt <= wavelength_mie[i+1]) {
            double t = (wavelength_rt - wavelength_mie[i]) / (wavelength_mie[i+1] - wavelength_mie[i]);
            return value_mie[i] + t * (value_mie[i+1] - value_mie[i]);
        }
    }
    
    return value_mie[n_mie-1]; // fallback
}

/**
 * Calculate cloud optical properties for all layers and wavelengths
 * 
 * This function:
 * 1. Interpolates Mie scattering data in particle size dimension (log-space)
 * 2. Interpolates in wavelength dimension to match EPACRIS RT grid
 * 3. Converts cloud mass density to opacity (cm^-1)
 * 
 * Inputs:
 * - particle_r2: particle radii from cloud physics (micrometers, volume-weighted) [zbin+1][MAX_CONDENSIBLES]
 * - clouds: cloud mass densities (molecules/cm³) [zbin+1][NSP+1]
 * - wavelength: EPACRIS wavelength grid (nm, from wavelength[] array)
 * - cH2O, aH2O, gH2O: output arrays [zbin+1][NLAMBDA]
 * 
 * Based on cloud_interp.c logic with HELIOS Mie data format
 */
void calculate_cloud_opacity_arrays(void) {
    // All arrays accessed as globals:
    // Inputs: clouds, particle_r2, particle_r0, particle_VP, particle_mass, wavelength
    // Outputs: cH2O, aH2O, gH2O
    
    const double sig = 2.0;  // Lognormal distribution width parameter
    
    const int H2O_species_id = 7;  // H2O species ID
    
    // Temporary arrays for interpolation
    double *ext_cross_interp = (double *)malloc(H2O_mie_optics.n_wavelengths * sizeof(double));
    double *albedo_interp = (double *)malloc(H2O_mie_optics.n_wavelengths * sizeof(double));
    double *asym_interp = (double *)malloc(H2O_mie_optics.n_wavelengths * sizeof(double));
    
    // Find H2O species index in CONDENSIBLES array
    int h2o_idx = -1;
    for (int k = 0; k < NCONDENSIBLES; k++) {
        if (CONDENSIBLES[k] == H2O_species_id) {
            h2o_idx = k;
            break;
        }
    }
    
    // Check if H2O is in condensibles list
    // Arrays are already initialized to zero in epacris_test_my_v2.c, so no need to zero here
    if (h2o_idx < 0) {
        // H2O not in condensibles list - exit early (arrays already initialized to zero)
        free(ext_cross_interp);
        free(albedo_interp);
        free(asym_interp);
        return;
    }
    
    // Process each layer as in the original code
    for (int j = 1; j <= zbin; j++) {
        // double r2 = particle_r2[j][h2o_idx];  // Volume-weighted radius (micrometers)
        // double cloud_mass_density = clouds[j][H2O_species_id];  // molecules/cm³
        
        // if cloud density very small, set to zero opacity and continue to next layer
        if (clouds[j][H2O_species_id] < 1e-12) {
            // No clouds or very small particles - set to zero opacity
            for (int i = 0; i < NLAMBDA; i++) {
                cH2O[j][i] = 0.0;
                aH2O[j][i] = 1.0;
                gH2O[j][i] = 0.0;
            }
            continue;
        }
        
        // Get mode radius (r0) for this layer and condensible species
        // r0 is the nucleation/monomer radius in micrometers [μm]
        // Note: Interpolation functions handle out-of-range values automatically
        // by clamping to the actual lookup table range
        double r0 = particle_r0[j][h2o_idx];
        
        // Particle number density: stored in particles/m³, convert to particles/cm³
        double particle_number_density_cm3 = particle_number_density[j][h2o_idx] * 1.0e-6;  // particles/m³ -> particles/cm³
        
        // Debug: Print interpolation parameters only for layers with non-zero particle size
        // Check both r0 and particle number density to ensure we only print for actual clouds
        if (r0 > 1e-10 && particle_number_density[j][h2o_idx] > 1e-10) {
            printf("\n=== CLOUD OPTICS INTERPOLATION DEBUG (Layer %d) ===\n", j);
            printf("particle_r0[%d][%d] = %.6e μm\n", j, h2o_idx, r0);
            printf("Lookup table radius range: %.6e - %.6e μm\n", 
                   H2O_mie_optics.radius[0], 
                   H2O_mie_optics.radius[H2O_mie_optics.n_particle_sizes-1]);
            printf("Number of particle sizes in lookup table: %d\n", H2O_mie_optics.n_particle_sizes);
            printf("particle_number_density[%d][%d] = %.6e particles/m³\n", j, h2o_idx, particle_number_density[j][h2o_idx]);
            printf("particle_number_density_cm3 = %.6e particles/cm³\n", particle_number_density_cm3);
        }


        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Interpolate optical properties in particle size dimension
        // Following cloud_interp.c rules (matching MATLAB cloud_interp.m):
        // - Extinction cross-section: log-log interpolation (log10 vs log10) - matches H2OI_c interpolation
        // - Albedo: linear in log-radius space - matches H2OI_a interpolation  
        // - Asymmetry: linear in log-radius space - matches H2OI_g interpolation
        
        // Interpolate for each Mie wavelength
        for (int w = 0; w < H2O_mie_optics.n_wavelengths; w++) {
            
            // Extract values for this wavelength across all particle sizes
            double *ext_per_r = (double *)malloc(H2O_mie_optics.n_particle_sizes * sizeof(double));
            double *alb_per_r = (double *)malloc(H2O_mie_optics.n_particle_sizes * sizeof(double));
            double *asym_per_r = (double *)malloc(H2O_mie_optics.n_particle_sizes * sizeof(double));
            
            for (int r = 0; r < H2O_mie_optics.n_particle_sizes; r++) {
                ext_per_r[r] = H2O_mie_optics.extinction_cross[r][w];
                alb_per_r[r] = H2O_mie_optics.albedo[r][w];
                asym_per_r[r] = H2O_mie_optics.asymmetry_g[r][w];
            }
            
            // Extinction cross-section: log-log interpolation (matches cloud_interp.c line 167: H2OI_c)
            ext_cross_interp[w] = interp1_loglog(H2O_mie_optics.radius, ext_per_r, 
                                                   H2O_mie_optics.n_particle_sizes, r0);
            
            // Albedo: linear in log-radius space (matches cloud_interp.c line 168: H2OI_a)
            albedo_interp[w] = interp1_log(H2O_mie_optics.radius, alb_per_r,
                                           H2O_mie_optics.n_particle_sizes, r0);
            
            // Asymmetry: linear in log-radius space (matches cloud_interp.c line 169: H2OI_g)
            asym_interp[w] = interp1_log(H2O_mie_optics.radius, asym_per_r,
                                        H2O_mie_optics.n_particle_sizes, r0);
            
            // Debug: Print interpolation results for layers with non-zero particle size (first 3 wavelengths)
            if (r0 > 1e-10 && particle_number_density[j][h2o_idx] > 1e-10 && w < 3) {
                // Find bracketing indices for r0
                int idx_low = -1, idx_high = -1;
                for (int r = 0; r < H2O_mie_optics.n_particle_sizes - 1; r++) {
                    if (r0 >= H2O_mie_optics.radius[r] && r0 <= H2O_mie_optics.radius[r+1]) {
                        idx_low = r;
                        idx_high = r + 1;
                        break;
                    }
                }
                
                printf("\n  Wavelength %d (%.6e μm):\n", w, H2O_mie_optics.wavelength[w]);
                if (idx_low >= 0) {
                    printf("    Lookup table radii: %.6e - %.6e μm (bracketing r0=%.6e)\n",
                           H2O_mie_optics.radius[idx_low], H2O_mie_optics.radius[idx_high], r0);
                    printf("    Extinction at radii: %.6e, %.6e cm²\n",
                           ext_per_r[idx_low], ext_per_r[idx_high]);
                    printf("    Albedo at radii: %.6f, %.6f\n",
                           alb_per_r[idx_low], alb_per_r[idx_high]);
                    printf("    Asymmetry at radii: %.6f, %.6f\n",
                           asym_per_r[idx_low], asym_per_r[idx_high]);
                } else if (r0 <= H2O_mie_optics.radius[0]) {
                    printf("    r0 below lookup range - using first table value\n");
                    printf("    Extinction: %.6e cm² (from radius %.6e μm)\n",
                           ext_per_r[0], H2O_mie_optics.radius[0]);
                } else {
                    printf("    r0 above lookup range - using last table value\n");
                    printf("    Extinction: %.6e cm² (from radius %.6e μm)\n",
                           ext_per_r[H2O_mie_optics.n_particle_sizes-1],
                           H2O_mie_optics.radius[H2O_mie_optics.n_particle_sizes-1]);
                }
                printf("    Interpolated extinction: %.6e cm²\n", ext_cross_interp[w]);
                printf("    Interpolated albedo: %.6f\n", albedo_interp[w]);
                printf("    Interpolated asymmetry: %.6f\n", asym_interp[w]);
            }
            
            free(ext_per_r);
            free(alb_per_r);
            free(asym_per_r);
        }
        
        // Debug: Print summary for layers with non-zero particle size
        if (r0 > 1e-10 && particle_number_density[j][h2o_idx] > 1e-10) {
            printf("\n  Interpolation complete for layer %d\n", j);
            printf("  Sample interpolated values (first 3 wavelengths):\n");
            for (int w = 0; w < 3 && w < H2O_mie_optics.n_wavelengths; w++) {
                printf("    λ=%.6e μm: ext=%.6e cm², alb=%.6f, asym=%.6f\n",
                       H2O_mie_optics.wavelength[w],
                       ext_cross_interp[w], albedo_interp[w], asym_interp[w]);
            }
            printf("================================================\n");
        }
        
        
        // Now interpolate in wavelength dimension to EPACRIS RT grid
        for (int i = 0; i < NLAMBDA; i++) {
            double wavelength_rt_micron = wavelength[i] * 1.0e-3;  // Convert nm to microns
            
            // Interpolate optical properties to RT wavelength grid
            double ext_cross = interp1_wavelength(H2O_mie_optics.wavelength, ext_cross_interp,
                                                  H2O_mie_optics.n_wavelengths, wavelength_rt_micron);
            double alb_val = interp1_wavelength(H2O_mie_optics.wavelength, albedo_interp,
                                               H2O_mie_optics.n_wavelengths, wavelength_rt_micron);
            double asym_val = interp1_wavelength(H2O_mie_optics.wavelength, asym_interp,
                                                H2O_mie_optics.n_wavelengths, wavelength_rt_micron);
            
            // Debug: Print wavelength interpolation results for layers with non-zero particle size (first 5 wavelengths)
            if (r0 > 1e-10 && particle_number_density[j][h2o_idx] > 1e-10 && i < 5) {
                // Find bracketing indices for wavelength
                int idx_low = -1, idx_high = -1;
                for (int w = 0; w < H2O_mie_optics.n_wavelengths - 1; w++) {
                    if (wavelength_rt_micron >= H2O_mie_optics.wavelength[w] && 
                        wavelength_rt_micron <= H2O_mie_optics.wavelength[w+1]) {
                        idx_low = w;
                        idx_high = w + 1;
                        break;
                    }
                }
                
                printf("\n  RT Wavelength %d (%.6e μm = %.6e nm):\n", i, wavelength_rt_micron, wavelength[i]);
                if (idx_low >= 0) {
                    printf("    MIE lookup wavelengths: %.6e - %.6e μm (bracketing λ_rt=%.6e)\n",
                           H2O_mie_optics.wavelength[idx_low], H2O_mie_optics.wavelength[idx_high], wavelength_rt_micron);
                    printf("    Extinction at wavelengths: %.6e, %.6e cm²\n",
                           ext_cross_interp[idx_low], ext_cross_interp[idx_high]);
                    printf("    Albedo at wavelengths: %.6f, %.6f\n",
                           albedo_interp[idx_low], albedo_interp[idx_high]);
                    printf("    Asymmetry at wavelengths: %.6f, %.6f\n",
                           asym_interp[idx_low], asym_interp[idx_high]);
                } else if (wavelength_rt_micron <= H2O_mie_optics.wavelength[0]) {
                    printf("    λ_rt below lookup range - using first table value\n");
                    printf("    Extinction: %.6e cm² (from wavelength %.6e μm)\n",
                           ext_cross_interp[0], H2O_mie_optics.wavelength[0]);
                } else {
                    printf("    λ_rt above lookup range - using last table value\n");
                    printf("    Extinction: %.6e cm² (from wavelength %.6e μm)\n",
                           ext_cross_interp[H2O_mie_optics.n_wavelengths-1],
                           H2O_mie_optics.wavelength[H2O_mie_optics.n_wavelengths-1]);
                }
                printf("    Interpolated extinction: %.6e cm²\n", ext_cross);
                printf("    Interpolated albedo: %.6f\n", alb_val);
                printf("    Interpolated asymmetry: %.6f\n", asym_val);
                printf("    Final opacity cH2O[%d][%d] = %.6e cm⁻¹ (n_density=%.6e * ext_cross=%.6e)\n",
                       j, i, particle_number_density_cm3 * ext_cross, particle_number_density_cm3, ext_cross);
            }
            
            // Calculate opacity (cm^-1) - matches MATLAB: crow = cloudden/VP * 1.0E-3 * pow(10, interp1(...))
            // cH2O = particle_number_density * extinction_cross_section
            cH2O[j][i] = particle_number_density_cm3 * ext_cross;  // cm^-1
            
            // Single scattering albedo (from MIE file, already interpolated)
            aH2O[j][i] = alb_val;
            
            // Asymmetry parameter (from MIE file, already interpolated)
            gH2O[j][i] = asym_val;
        }
        
        // Debug: Print summary for wavelength interpolation
        if (r0 > 1e-10 && particle_number_density[j][h2o_idx] > 1e-10) {
            printf("\n  Wavelength interpolation complete for layer %d\n", j);
            printf("  Sample final RT values (first 5 wavelengths):\n");
            for (int i = 0; i < 5 && i < NLAMBDA; i++) {
                printf("    λ=%.6e nm: cH2O=%.6e cm⁻¹, aH2O=%.6f, gH2O=%.6f\n",
                       wavelength[i], cH2O[j][i], aH2O[j][i], gH2O[j][i]);
            }
            printf("================================================\n");
        }
    }
    
    free(ext_cross_interp);
    free(albedo_interp);
    free(asym_interp);
}

