/**
 * EPACRIS Cloud Physics Integration Example
 * 
 * This file demonstrates how to integrate the new cloud distribution physics
 * into your existing EPACRIS convection routine (ms_conv_funcs.c)
 * 
 * INTEGRATION STEPS:
 * 1. Add #include "cloud_physics.h" to your main files
 * 2. Replace simple ALPHA_RAINOUT with calculate_cloud_distribution()
 * 3. Update radiative transfer to use enhanced cloud optical properties
 * 4. Add configuration constants to AlphaAb.h
 */

#include <stdio.h>
#include <math.h>
#include "constant.h"
#include "cloud_physics.h"

//=========================================================
//=== Integration into ms_conv_funcs.c ===================
//=========================================================

void enhanced_convection_step(int layer) {
    /**
     * Enhanced convection step that replaces the simple rainout physics
     * with proper cloud distribution calculations
     * 
     * REPLACE THIS SECTION IN ms_conv_funcs.c:
     * 
     * OLD CODE (around line 500-550):
     * if (PRESSURE_CONSERVATION == 1) {
     *     for (int i = 0; i < NCONDENSIBLES; i++) {
     *         int species_id = CONDENSIBLES[i];
     *         clouds[layer][species_id] *= ALPHA_RAINOUT;
     *     }
     * }
     * 
     * NEW CODE:
     */
    
    if (PRESSURE_CONSERVATION == 1) {
        // Calculate physically realistic cloud distribution
        calculate_cloud_distribution(layer);
        
        // Optional: Apply sedimentation if you want time-dependent transport
        // apply_sedimentation(layer, dt_timestep);
        
        // Write diagnostics periodically
        static int diagnostic_counter = 0;
        if (diagnostic_counter % 100 == 0) {
            char filename[256];
            sprintf(filename, "cloud_diagnostics_%d.dat", diagnostic_counter);
            write_cloud_diagnostics(diagnostic_counter, filename);
        }
        diagnostic_counter++;
    }
}

//=========================================================
//=== Integration into ms_radtrans_test.c ================
//=========================================================

void enhanced_radiative_transfer(int layer, int wavelength_index) {
    /**
     * Enhanced radiative transfer that uses the new cloud optical properties
     * 
     * REPLACE THIS SECTION IN ms_radtrans_test.c (around lines 250-270):
     * 
     * OLD CODE:
     * // wa[j] += cH2O[j1][i]*(1.0-aH2O[j1][i])/MM[j1];  // COMMENTED OUT
     * // wa[j] += cNH3[j1][i]*(1.0-aNH3[j1][i])/MM[j1];  // COMMENTED OUT
     * 
     * NEW CODE:
     */
    
    double total_cloud_absorption = 0.0;
    double total_cloud_scattering = 0.0;
    
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        if (particle_radius[layer][i] > 0.0) {
            // Get cloud optical properties from the cloud physics calculation
            double optical_depth = cloud_optical_depth[layer][i];
            double single_scatter_albedo = cloud_single_scattering_albedo[layer][i];
            double asymmetry_param = cloud_asymmetry_parameter[layer][i];
            
            // Calculate absorption and scattering contributions
            double absorption_opacity = optical_depth * (1.0 - single_scatter_albedo);
            double scattering_opacity = optical_depth * single_scatter_albedo;
            
            total_cloud_absorption += absorption_opacity / MM[layer];
            total_cloud_scattering += scattering_opacity / MM[layer];
            
            // Debug output for water clouds
            if (species_id == 7 && layer == 60) { // H2O at mid-atmosphere
                printf("Layer %d H2O cloud: r=%.2f um, tau=%.6f, ssa=%.3f\n",
                       layer, particle_radius[layer][i]*1e6, optical_depth, single_scatter_albedo);
            }
        }
    }
    
    // Add to total absorption coefficient
    // wa[wavelength_index] += total_cloud_absorption;
    
    // Add to total scattering coefficient (if you have scattering implemented)
    // ws[wavelength_index] += total_cloud_scattering;
}

//=========================================================
//=== Configuration Constants for AlphaAb.h ==============
//=========================================================

/**
 * ADD THESE CONSTANTS TO YOUR AlphaAb.h FILE:
 * 
 * // Cloud Physics Configuration
 * #define CLOUD_DISTRIBUTION_MODE 1        // 0=simple, 1=Ackerman-Marley
 * #define EDDY_DIFFUSION_COEFF 1.0e6      // m²/s (typical for exoplanets)
 * #define MIN_PARTICLE_RADIUS 1.0e-8      // m (10 nm minimum)
 * #define MAX_PARTICLE_RADIUS 1.0e-4      // m (100 um maximum)
 * #define MIN_SEDIMENTATION_VELOCITY 1.0e-6 // m/s
 * #define MAX_SEDIMENTATION_VELOCITY 1.0e1  // m/s
 * #define MASS_DIFFUSION_COEFF 1.0e-5     // m²/s (molecular diffusion)
 * #define ACCOMMODATION_COEFF 0.1          // Accommodation coefficient
 * 
 * // Physical Constants (if not already defined)
 * #define KB 1.380649e-23                 // Boltzmann constant (J/K)
 * #define AMU 1.66053906660e-27            // Atomic mass unit (kg)
 * #define GRAVITY 9.81                     // Surface gravity (m/s²) - UPDATE FOR YOUR PLANET
 */

//=========================================================
//=== Example Usage in Main Program ======================
//=========================================================

void example_main_integration(void) {
    /**
     * Example of how to integrate into your main EPACRIS loop
     */
    
    printf("=== EPACRIS Enhanced Cloud Physics ===\n");
    printf("Initializing cloud distribution physics...\n");
    
    // Initialize cloud physics arrays (done automatically when arrays are declared)
    
    // Main convection loop (this replaces your existing loop in epacris_test_my_v2.c)
    for (int iteration = 0; iteration < 1000; iteration++) {
        
        // Your existing temperature/pressure iteration
        // ... existing EPACRIS code ...
        
        // Enhanced convection with cloud physics
        for (int layer = 1; layer <= zbin; layer++) {
            
            // Your existing convection calculations
            // ... existing ms_conv_funcs calls ...
            
            // NEW: Enhanced cloud distribution
            if (NCONDENSIBLES > 0) {
                enhanced_convection_step(layer);
            }
        }
        
        // Your existing radiative transfer
        for (int layer = 1; layer <= zbin; layer++) {
            for (int wl = 0; wl < NWAVE; wl++) {
                
                // Your existing gas opacity calculations
                // ... existing opacity calculations ...
                
                // NEW: Enhanced cloud opacity
                enhanced_radiative_transfer(layer, wl);
            }
        }
        
        // Convergence check and output
        if (iteration % 50 == 0) {
            printf("Iteration %d: Enhanced cloud physics active\n", iteration);
            
            // Output cloud diagnostics
            char filename[256];
            sprintf(filename, "clouds_iter_%d.dat", iteration);
            write_cloud_diagnostics(iteration, filename);
        }
    }
    
    printf("=== Cloud Physics Integration Complete ===\n");
}

//=========================================================
//=== Compilation Instructions ===========================
//=========================================================

/**
 * COMPILATION:
 * 
 * 1. Add to your Makefile:
 *    SOURCES += cloud_physics.c
 *    HEADERS += cloud_physics.h
 * 
 * 2. Compile with:
 *    gcc -o epacris_enhanced epacris_test_my_v2.c ms_conv_funcs.c ms_radtrans_test.c cloud_physics.c -lm
 * 
 * 3. Make sure all header files are properly included:
 *    - Add #include "cloud_physics.h" to epacris_test_my_v2.c
 *    - Add #include "cloud_physics.h" to ms_conv_funcs.c  
 *    - Add #include "cloud_physics.h" to ms_radtrans_test.c
 */

//=========================================================
//=== Expected Results ===================================
//=========================================================

/**
 * WHAT TO EXPECT:
 * 
 * 1. Cloud Distribution:
 *    - Clouds will be more concentrated near condensation levels
 *    - Less spreading over orders of magnitude in pressure
 *    - Realistic particle sizes (1-100 microns typically)
 * 
 * 2. Sedimentation:
 *    - Larger particles fall faster
 *    - Smaller particles stay suspended longer
 *    - Physically consistent vertical transport
 * 
 * 3. Optical Properties:
 *    - Size-dependent scattering and absorption
 *    - Proper single scattering albedo
 *    - Realistic asymmetry parameters
 * 
 * 4. Mass Conservation:
 *    - Total condensate mass preserved
 *    - No artificial loss or gain
 *    - Warnings for any conservation violations
 * 
 * 5. Diagnostics:
 *    - cloud_diagnostics_*.dat files with particle sizes, fall velocities
 *    - Debug output for key atmospheric layers
 *    - Convergence monitoring
 */ 