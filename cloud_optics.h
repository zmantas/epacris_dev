/**
 * cloud_optics.h
 * 
 * Header file for cloud optical property lookup tables
 * Reads Mie scattering tables from Clouds/CrossP/ directory
 * 
 * Tables structure:
 * - Row dimension: particle radius (41 rows)
 * - Column dimension: wavelength bins (NLAMBDA)
 * - Files: Cross (extinction), Albedo (single scattering), Geo (asymmetry parameter g)
 */

#ifndef CLOUD_OPTICS_H
#define CLOUD_OPTICS_H

#include "constant.h"

// Number of particle size bins in lookup tables
#define NPARTICLE_SIZES 41

// Number of wavelengths in cloud optical property lookup tables
// Note: This is DIFFERENT from NLAMBDA (radiative transfer grid)
// Cloud tables have 1387 wavelengths, RT grid may have 16000
#define NCLOUD_WAVELENGTHS 1387

// Cloud optical property data structures
// Format: [particle_size_index][wavelength_index]
typedef struct {
    double radius[NPARTICLE_SIZES];                              // Particle radii in micrometers
    double cross[NPARTICLE_SIZES][NCLOUD_WAVELENGTHS];          // Cross section (extinction)
    double albedo[NPARTICLE_SIZES][NCLOUD_WAVELENGTHS];         // Single scattering albedo
    double asym[NPARTICLE_SIZES][NCLOUD_WAVELENGTHS];           // Asymmetry parameter g
} CloudOpticalTable;

// Global optical property tables (defined in cloud_optics.c)
extern CloudOpticalTable H2O_ice_optics;
extern CloudOpticalTable H2O_liquid_optics;
extern CloudOpticalTable NH3_ice_optics;

// Function declarations
void read_cloud_optical_tables(void);
void cleanup_cloud_optical_tables(void);

// Future interpolation function (to be implemented in next step)
// void get_cloud_opacity(int layer, int wavelength_idx, int species_id,
//                        double particle_size_um,
//                        double *extinction, double *albedo, double *asymmetry);

#endif // CLOUD_OPTICS_H

