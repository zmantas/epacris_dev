/**
 * test_cloud_optics.c
 * 
 * Simple test program to verify cloud optical property tables are read correctly
 * Compile: gcc -o test_cloud_optics test_cloud_optics.c -lm
 */

#include <stdio.h>
#include <stdlib.h>

// Define constants needed by cloud_optics.h
#define NLAMBDA 196
#define zbin 100

#include "cloud_optics.c"

int main() {
    printf("\n");
    printf("========================================\n");
    printf("Testing Cloud Optical Property Reading\n");
    printf("========================================\n\n");
    
    // Read the optical property tables
    read_cloud_optical_tables();
    
    // Print some diagnostic information
    printf("\n========================================\n");
    printf("Verification of Loaded Data\n");
    printf("========================================\n\n");
    
    printf("--- H2O (Water) Optical Properties ---\n");
    printf("Particle size range: %.4f - %.4f micrometers\n", 
           H2O_ice_optics.radius[0], 
           H2O_ice_optics.radius[NPARTICLE_SIZES-1]);
       for (int i = 0; i < NPARTICLE_SIZES; i++) {
           printf("radius[%d] = %f\n", i, H2O_ice_optics.radius[i]);
       }
    printf("Sample values at particle size %.4f um:\n", H2O_ice_optics.radius[10]);
    printf("  Wavelength bin 0:   cross=%.4e, albedo=%.4f, asym=%.4f\n",
           H2O_ice_optics.cross[10][0],
           H2O_ice_optics.albedo[10][0],
           H2O_ice_optics.asym[10][0]);
    printf("  Wavelength bin 50:  cross=%.4e, albedo=%.4f, asym=%.4f\n",
           H2O_ice_optics.cross[10][50],
           H2O_ice_optics.albedo[10][50],
           H2O_ice_optics.asym[10][50]);
    printf("  Wavelength bin 100: cross=%.4e, albedo=%.4f, asym=%.4f\n",
           H2O_ice_optics.cross[10][100],
           H2O_ice_optics.albedo[10][100],
           H2O_ice_optics.asym[10][100]);
    
    printf("\n--- NH3 (Ammonia) Optical Properties ---\n");
    printf("Particle size range: %.4f - %.4f micrometers\n", 
           NH3_ice_optics.radius[0], 
           NH3_ice_optics.radius[NPARTICLE_SIZES-1]);
    printf("Sample values at particle size %.4f um:\n", NH3_ice_optics.radius[10]);
    printf("  Wavelength bin 0:   cross=%.4e, albedo=%.4f, asym=%.4f\n",
           NH3_ice_optics.cross[10][0],
           NH3_ice_optics.albedo[10][0],
           NH3_ice_optics.asym[10][0]);
    printf("  Wavelength bin 50:  cross=%.4e, albedo=%.4f, asym=%.4f\n",
           NH3_ice_optics.cross[10][50],
           NH3_ice_optics.albedo[10][50],
           NH3_ice_optics.asym[10][50]);
    printf("  Wavelength bin 100: cross=%.4e, albedo=%.4f, asym=%.4f\n",
           NH3_ice_optics.cross[10][100],
           NH3_ice_optics.albedo[10][100],
           NH3_ice_optics.asym[10][100]);
    
    printf("\n========================================\n");
    printf("Test completed successfully!\n");
    printf("========================================\n\n");
    
    // Cleanup
    cleanup_cloud_optical_tables();
    
    return 0;
}

