#ifndef CLOUD_PHYSICS_H
#define CLOUD_PHYSICS_H

#include "constant.h"

//=========================================================
//=== Cloud Physics Constants ============================
//=========================================================

// Particle size distribution parameters
#define CLOUD_SIGMA_DEFAULT 2.0        // Log-normal distribution width (dimensionless)
#define PARTICLE_DENSITY_H2O_ICE 917.0 // kg/m³ (ice)
#define PARTICLE_DENSITY_H2O_LIQ 1000.0 // kg/m³ (liquid water)
#define PARTICLE_DENSITY_NH3_ICE 817.0 // kg/m³ (ammonia ice)
#define PARTICLE_DENSITY_CO2_ICE 1563.0 // kg/m³ (dry ice)
#define PARTICLE_DENSITY_CH4_ICE 424.0 // kg/m³ (methane ice)
#define PARTICLE_DENSITY_CO_ICE 789.0  // kg/m³ (carbon monoxide ice)
#define PARTICLE_DENSITY_H2_ICE 86.0   // kg/m³ (hydrogen ice)
#define PARTICLE_DENSITY_N2_ICE 1026.0 // kg/m³ (nitrogen ice)
#define PARTICLE_DENSITY_O2_ICE 1426.0 // kg/m³ (oxygen ice)
#define PARTICLE_DENSITY_H2S_ICE 993.0 // kg/m³ (hydrogen sulfide ice)

// Microphysics parameters
#define MIN_PARTICLE_RADIUS 0.01e-6    // m (minimum 0.01 microns)
#define MAX_PARTICLE_RADIUS 1000.0e-6  // m (maximum 1000 microns)
#define ACCOMMODATION_COEFF 1.0        // Condensation accommodation coefficient
#define MASS_DIFFUSION_COEFF 0.12e-4   // m²/s (mass diffusion coefficient)

// Sedimentation parameters
#define MIN_SEDIMENTATION_VELOCITY 1.0e-8  // m/s (minimum fall velocity)
#define MAX_SEDIMENTATION_VELOCITY 10.0    // m/s (maximum fall velocity)
#define STOKES_CUNNINGHAM_TRANSITION 2.0   // Knudsen number transition

// Cloud distribution parameters
#define CLOUD_DISTRIBUTION_MODE 1      // 0=uniform, 1=Ackerman-Marley, 2=microphysics-based
#define EDDY_DIFFUSION_COEFF 1.0e4     // m²/s (default eddy diffusion coefficient)
#define MIXING_LENGTH_SCALE 1000.0     // m (typical mixing length scale)

//=========================================================
//=== Cloud Physics Structures ===========================
//=========================================================

// Structure to hold particle size distribution information
typedef struct {
    double r_mode;          // Mode radius (m)
    double r_effective;     // Effective radius (m) 
    double sigma;           // Log-normal width parameter
    double number_density;  // Particle number density (m⁻³)
    double mass_density;    // Particle mass density (kg/m³)
    double fall_velocity;   // Sedimentation velocity (m/s)
} ParticleDistribution;

// Structure to hold cloud microphysics parameters
typedef struct {
    double condensation_rate;    // kg/(m³·s)
    double evaporation_rate;     // kg/(m³·s)
    double nucleation_rate;      // particles/(m³·s)
    double coagulation_rate;     // s⁻¹
    double sedimentation_flux;   // kg/(m²·s)
} CloudMicrophysics;

// Structure to hold atmospheric layer properties for cloud calculations
typedef struct {
    double temperature;     // K
    double pressure;        // Pa
    double density;         // kg/m³
    double mean_free_path;  // m
    double viscosity;       // Pa·s
    double scale_height;    // m
    double eddy_diffusion;  // m²/s
} AtmosphericLayer;

//=========================================================
//=== Function Declarations ==============================
//=========================================================

// Core cloud distribution functions
void calculate_cloud_distribution(int layer);
void update_particle_sizes(int layer);
void apply_sedimentation(int layer, double dt);
void calculate_optical_properties(int layer);

// Particle size calculation functions
double calculate_particle_radius_microphysics(int species_id, int layer, double condensate_mass_density);
double calculate_particle_radius_simple(int species_id, double condensate_mass_density, double particle_density);
void calculate_lognormal_distribution(ParticleDistribution* dist, double total_mass, double number_density);

// Sedimentation velocity functions
double calculate_fall_velocity(double particle_radius, double particle_density, AtmosphericLayer* atm_layer);
double calculate_reynolds_correction(double particle_radius, double fall_velocity, AtmosphericLayer* atm_layer);
double calculate_cunningham_correction(double particle_radius, double mean_free_path);

// Cloud microphysics functions
double calculate_condensation_rate(int species_id, int layer, double supersaturation, double surface_area);
double calculate_nucleation_rate(int species_id, int layer, double supersaturation);
double calculate_coagulation_rate(ParticleDistribution* dist, AtmosphericLayer* atm_layer);

// Atmospheric property functions
void calculate_atmospheric_properties(int layer, AtmosphericLayer* atm_layer);
double calculate_mean_free_path(double temperature, double pressure, double molecular_mass);
double calculate_dynamic_viscosity(double temperature);
double calculate_eddy_diffusion_coefficient(int layer);

// Mass conservation functions
void enforce_mass_conservation(int layer);
void redistribute_condensate_mass(int layer, double* mass_change);
void update_vapor_phase_after_sedimentation(int layer, double sedimentation_loss);

// Utility functions
double get_particle_density(int species_id, double temperature);
double interpolate_optical_properties(double particle_radius, int wavelength_index, int species_id);
void write_cloud_diagnostics(int iteration, const char* filename);

// Global arrays for cloud distribution (declared as extern, defined in main file)
extern ParticleDistribution particle_dist[zbin+1][NCONDENSIBLES];
extern CloudMicrophysics cloud_micro[zbin+1][NCONDENSIBLES];
extern AtmosphericLayer atm_properties[zbin+1];

// Global arrays for enhanced cloud optical properties
extern double particle_radius[zbin+1][NCONDENSIBLES];      // Effective radius (m)
extern double particle_number_density[zbin+1][NCONDENSIBLES]; // Number density (m⁻³)
extern double sedimentation_velocity[zbin+1][NCONDENSIBLES];   // Fall velocity (m/s)
extern double cloud_optical_depth[zbin+1][NCONDENSIBLES];     // Optical depth
extern double cloud_single_scattering_albedo[zbin+1][NCONDENSIBLES]; // Single scattering albedo
extern double cloud_asymmetry_parameter[zbin+1][NCONDENSIBLES];      // Asymmetry parameter

#endif // CLOUD_PHYSICS_H 