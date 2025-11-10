/**
 * cloud_optics.c
 * 
 * Functions to read and store cloud optical property lookup tables
 * Based on Mie scattering calculations for water and ammonia clouds
 * 
 * Data files from Clouds/CrossP/:
 * - Cross_water_wavelength2.dat      : H2O extinction cross sections
 * - Albedo_water_wavelength2.dat     : H2O single scattering albedo
 * - Geo_water_wavelength2.dat        : H2O asymmetry parameter g
 * - Cross_ammonia_wavelength2.dat    : NH3 extinction cross sections
 * - Albedo_ammonia_wavelength2.dat   : NH3 single scattering albedo
 * - Geo_ammonia_wavelength2.dat      : NH3 asymmetry parameter g
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cloud_optics.h"

// Global optical property tables
CloudOpticalTable H2O_ice_optics;
CloudOpticalTable H2O_liquid_optics;
CloudOpticalTable NH3_ice_optics;

/**
 * Helper function to read a single optical property file
 * 
 * File format (each row is one particle size):
 * - Column 1: Particle radius in micrometers
 * - Columns 2-1388: Optical properties at 1387 wavelengths (total 1388 columns)
 * 
 * We read all 1387 wavelengths and store them at native resolution.
 * Interpolation to RT wavelength grid (NLAMBDA) happens during radiative transfer.
 */
static int read_optical_property_file(const char *filename, 
                                      double radius[NPARTICLE_SIZES],
                                      double data[NPARTICLE_SIZES][NCLOUD_WAVELENGTHS]) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: Cannot open cloud optical property file: %s\n", filename);
        return -1;
    }
    
    printf("Reading cloud optical properties from: %s\n", filename);
    
    // Read data row by row (each row = one particle size)
    for (int i = 0; i < NPARTICLE_SIZES; i++) {
        // Column 1: Read particle radius (micrometers)
        if (fscanf(fp, "%lf", &radius[i]) != 1) {
            printf("ERROR: Failed to read particle radius at row %d in %s\n", i, filename);
            fclose(fp);
            return -1;
        }
        
        // Columns 2-1388: Read all 1387 wavelength values
        for (int j = 0; j < NCLOUD_WAVELENGTHS; j++) {
            if (fscanf(fp, "%lf", &data[i][j]) != 1) {
                printf("ERROR: Failed to read optical data at row %d, wavelength %d in %s\n", 
                       i, j, filename);
                fclose(fp);
                return -1;
            }
        }
    }
    
    fclose(fp);
    printf("Successfully read %d particle sizes x %d wavelengths from %s\n", 
           NPARTICLE_SIZES, NCLOUD_WAVELENGTHS, filename);
    return 0;
}

/**
 * Read all cloud optical property lookup tables
 * Should be called once during initialization in main program
 */
void read_cloud_optical_tables(void) {
    printf("\n");
    printf("=======================================================\n");
    printf("READING CLOUD OPTICAL PROPERTY LOOKUP TABLES\n");
    printf("=======================================================\n");
    
    int status = 0;
    
    // ===== WATER CLOUD OPTICS (H2O) =====
    // Note: Using same tables for ice and liquid - distinguish by temperature later
    printf("\n--- Water Cloud Optical Properties ---\n");
    
    // Read cross sections (extinction) - also populates radius array
    status = read_optical_property_file(
        "Clouds/CrossP/Cross_water_wavelength2.dat",
        H2O_ice_optics.radius,
        H2O_ice_optics.cross
    );
    if (status != 0) {
        printf("FATAL ERROR: Failed to read water cross sections\n");
        exit(1);
    }
    
    // Read albedo (single scattering albedo) - uses same radius array
    status = read_optical_property_file(
        "Clouds/CrossP/Albedo_water_wavelength2.dat",
        H2O_ice_optics.radius,
        H2O_ice_optics.albedo
    );
    if (status != 0) {
        printf("FATAL ERROR: Failed to read water albedo\n");
        exit(1);
    }
    
    // Read asymmetry parameter g - uses same radius array
    status = read_optical_property_file(
        "Clouds/CrossP/Geo_water_wavelength2.dat",
        H2O_ice_optics.radius,
        H2O_ice_optics.asym
    );
    if (status != 0) {
        printf("FATAL ERROR: Failed to read water asymmetry parameter\n");
        exit(1);
    }
    
    // Copy to liquid optics (can differentiate later if needed)
    memcpy(&H2O_liquid_optics, &H2O_ice_optics, sizeof(CloudOpticalTable));
    
    // ===== AMMONIA CLOUD OPTICS (NH3) =====
    printf("\n--- Ammonia Cloud Optical Properties ---\n");
    
    // Read cross sections (extinction) - also populates radius array
    status = read_optical_property_file(
        "Clouds/CrossP/Cross_ammonia_wavelength2.dat",
        NH3_ice_optics.radius,
        NH3_ice_optics.cross
    );
    if (status != 0) {
        printf("FATAL ERROR: Failed to read ammonia cross sections\n");
        exit(1);
    }
    
    // Read albedo - uses same radius array
    status = read_optical_property_file(
        "Clouds/CrossP/Albedo_ammonia_wavelength2.dat",
        NH3_ice_optics.radius,
        NH3_ice_optics.albedo
    );
    if (status != 0) {
        printf("FATAL ERROR: Failed to read ammonia albedo\n");
        exit(1);
    }
    
    // Read asymmetry parameter g - uses same radius array
    status = read_optical_property_file(
        "Clouds/CrossP/Geo_ammonia_wavelength2.dat",
        NH3_ice_optics.radius,
        NH3_ice_optics.asym
    );
    if (status != 0) {
        printf("FATAL ERROR: Failed to read ammonia asymmetry parameter\n");
        exit(1);
    }
    
    printf("\n");
    printf("=======================================================\n");
    printf("CLOUD OPTICAL TABLES SUCCESSFULLY LOADED\n");
    printf("=======================================================\n");
    printf("Particle size range: %.4f - %.4f micrometers\n", 
           H2O_ice_optics.radius[0], H2O_ice_optics.radius[NPARTICLE_SIZES-1]);
    printf("Number of wavelengths in lookup tables: %d\n", NCLOUD_WAVELENGTHS);
    printf("Species loaded: H2O (ice/liquid), NH3 (ice)\n");
    printf("Note: RT grid has NLAMBDA=%d; interpolation needed during RT\n", NLAMBDA);
    printf("=======================================================\n\n");
}

/**
 * Cleanup function (if needed for future memory management)
 * Currently optical tables are static, so no cleanup needed
 */
void cleanup_cloud_optical_tables(void) {
    // No dynamic memory to free currently
    // This function is provided for future extension
}

