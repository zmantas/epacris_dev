#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constant.h"
#include "cloud_physics.h"

//=========================================================
//=== Global Arrays for Cloud Physics ====================
//=========================================================

// Arrays for enhanced cloud distribution (defined here, declared as extern in header)
ParticleDistribution particle_dist[zbin+1][NCONDENSIBLES];
CloudMicrophysics cloud_micro[zbin+1][NCONDENSIBLES];
AtmosphericLayer atm_properties[zbin+1];

// Enhanced cloud optical properties
double particle_radius[zbin+1][NCONDENSIBLES];
double particle_number_density[zbin+1][NCONDENSIBLES];
double sedimentation_velocity[zbin+1][NCONDENSIBLES];
double cloud_optical_depth[zbin+1][NCONDENSIBLES];
double cloud_single_scattering_albedo[zbin+1][NCONDENSIBLES];
double cloud_asymmetry_parameter[zbin+1][NCONDENSIBLES];

// External function declarations for saturation pressure functions
extern double ms_psat_h2o(double temp);
extern double ms_psat_nh3(double temp);
extern double ms_psat_co(double temp);
extern double ms_psat_ch4(double temp);
extern double ms_psat_co2(double temp);
extern double ms_psat_h2(double temp);
extern double ms_psat_o2(double temp);
extern double ms_psat_n2(double temp);
extern double ms_psat_h2s(double temp);

// Forward declarations for functions used before definition
void apply_ackerman_marley_distribution(int layer, int species_index);
double get_molecular_mass(int species_id);

//=========================================================
//=== Core Cloud Distribution Functions ==================
//=========================================================

void calculate_cloud_distribution(int layer) {
    /**
     * Main function to calculate physically realistic cloud distribution
     * Implements Ackerman & Marley (2001) + improved microphysics
     * 
     * This replaces the simple ALPHA_RAINOUT approach with proper physics
     */
    
    // Skip if no condensibles or if layer is out of bounds
    if (NCONDENSIBLES == 0 || layer < 1 || layer > zbin) {
        return;
    }
    
    // Calculate atmospheric properties for this layer
    calculate_atmospheric_properties(layer, &atm_properties[layer]);
    
    // Process each condensible species
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        // Get current condensate mass density (kg/m³)
        double condensate_mass_density = clouds[layer][species_id] / MM[layer] * 
                                       atm_properties[layer].density;
        
        if (condensate_mass_density < 1.0e-20) {
            // No condensate - reset all properties
            particle_radius[layer][i] = 0.0;
            particle_number_density[layer][i] = 0.0;
            sedimentation_velocity[layer][i] = 0.0;
            continue;
        }
        
        // STEP 1: Calculate particle size using microphysics
        update_particle_sizes(layer);
        
        // STEP 2: Calculate sedimentation velocity
        double particle_density = get_particle_density(species_id, tl[layer]);
        sedimentation_velocity[layer][i] = calculate_fall_velocity(
            particle_radius[layer][i], particle_density, &atm_properties[layer]);
        
        // STEP 3: Apply Ackerman & Marley cloud distribution
        apply_ackerman_marley_distribution(layer, i);
        
        // STEP 4: Calculate optical properties for radiative transfer
        calculate_optical_properties(layer);
    }
    
    // STEP 5: Enforce mass conservation
    enforce_mass_conservation(layer);
}

void update_particle_sizes(int layer) {
    /**
     * Calculate particle sizes using microphysics approach
     * Based on the particlesizef.c algorithm from your Clouds/ folder
     */
    
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        // Get condensate mass density
        double condensate_mass_density = clouds[layer][species_id] / MM[layer] * 
                                       atm_properties[layer].density;
        
        if (condensate_mass_density < 1.0e-20) {
            particle_radius[layer][i] = 0.0;
            continue;
        }
        
        // Calculate particle radius using microphysics
        particle_radius[layer][i] = calculate_particle_radius_microphysics(
            species_id, layer, condensate_mass_density);
        
        // Calculate number density from mass and particle size
        double particle_density = get_particle_density(species_id, tl[layer]);
        double single_particle_mass = (4.0/3.0) * PI * 
            pow(particle_radius[layer][i], 3.0) * particle_density;
        
        if (single_particle_mass > 0.0) {
            particle_number_density[layer][i] = condensate_mass_density / single_particle_mass;
        } else {
            particle_number_density[layer][i] = 0.0;
        }
        
        // Store in particle distribution structure
        particle_dist[layer][i].r_effective = particle_radius[layer][i];
        particle_dist[layer][i].number_density = particle_number_density[layer][i];
        particle_dist[layer][i].mass_density = condensate_mass_density;
        particle_dist[layer][i].sigma = CLOUD_SIGMA_DEFAULT;
    }
}

void apply_ackerman_marley_distribution(int layer, int species_index) {
    /**
     * Apply Ackerman & Marley (2001) cloud distribution physics
     * This determines how much condensate stays vs falls out
     */
    
    int species_id = CONDENSIBLES[species_index];
    
    // Get atmospheric properties
    double eddy_diffusion = atm_properties[layer].eddy_diffusion;
    double scale_height = atm_properties[layer].scale_height;
    double fall_velocity = sedimentation_velocity[layer][species_index];
    
    // Calculate Ackerman & Marley f_sed parameter
    // f_sed = v_sed * H / K_zz
    double f_sed = (fall_velocity * scale_height) / eddy_diffusion;
    
    // Limit f_sed to reasonable range
    f_sed = fmax(0.001, fmin(f_sed, 100.0));
    
    // Calculate cloud retention factor (what stays in the layer)
    // This replaces the simple ALPHA_RAINOUT with physics-based calculation
    double retention_factor;
    
    if (CLOUD_DISTRIBUTION_MODE == 1) {
        // Ackerman & Marley approach: retention = 1 / (1 + f_sed)
        retention_factor = 1.0 / (1.0 + f_sed);
    } else {
        // Simple uniform distribution (fallback)
        retention_factor = ALPHA_RAINOUT;
    }
    
    // Apply bounds to prevent numerical issues
    retention_factor = fmax(0.01, fmin(retention_factor, 0.99));
    
    // Update cloud amount based on retention
    double original_cloud = clouds[layer][species_id];
    clouds[layer][species_id] = original_cloud * retention_factor;
    
    // Calculate sedimentation flux for mass conservation
    double sedimentation_flux = original_cloud * (1.0 - retention_factor);
    cloud_micro[layer][species_index].sedimentation_flux = sedimentation_flux;
    
    // Debug output for key layers
    if (layer == 60 && species_id == 7) { // Mid-atmosphere H2O
        printf("Ackerman-Marley Layer %d, Species %d:\n", layer, species_id);
        printf("  f_sed = %.6f\n", f_sed);
        printf("  retention_factor = %.6f\n", retention_factor);
        printf("  fall_velocity = %.6e m/s\n", fall_velocity);
        printf("  eddy_diffusion = %.6e m²/s\n", eddy_diffusion);
        printf("  original_cloud = %.6e\n", original_cloud);
        printf("  new_cloud = %.6e\n", clouds[layer][species_id]);
    }
}

void apply_sedimentation(int layer, double dt) {
    /**
     * Apply sedimentation physics to transport condensate downward
     * This handles the vertical redistribution of cloud material
     */
    
    if (layer <= 1) return; // Skip bottom layer
    
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        // Get sedimentation flux from layer above
        double flux_from_above = 0.0;
        if (layer < zbin) {
            flux_from_above = cloud_micro[layer+1][i].sedimentation_flux;
        }
        
        // Calculate layer thickness (convert pressure difference to height)
        double layer_thickness = atm_properties[layer].scale_height * 
            log(P[layer-1] / P[layer]);
        
        // Calculate mass change due to sedimentation
        double mass_gain_from_above = flux_from_above * dt / layer_thickness;
        double mass_loss_downward = cloud_micro[layer][i].sedimentation_flux * dt / layer_thickness;
        
        double net_mass_change = mass_gain_from_above - mass_loss_downward;
        
        // Update cloud abundance
        clouds[layer][species_id] += net_mass_change * MM[layer] / atm_properties[layer].density;
        
        // Ensure non-negative values
        if (clouds[layer][species_id] < 0.0) {
            clouds[layer][species_id] = 0.0;
        }
    }
}

//=========================================================
//=== Particle Size Calculation Functions ===============
//=========================================================

double calculate_particle_radius_microphysics(int species_id, int layer, double condensate_mass_density) {
    /**
     * Calculate particle radius using microphysics approach
     * Based on the particlesizef algorithm from your Clouds/ folder
     * This solves the balance between condensation growth and sedimentation
     */
    
    // Get atmospheric properties
    AtmosphericLayer* atm = &atm_properties[layer];
    
    // Calculate supersaturation
    double partial_pressure = (xx[layer][species_id] / MM[layer]) * pl[layer];
    double saturation_pressure;
    
    // Get saturation pressure for the species
    switch(species_id) {
        case 7:  saturation_pressure = ms_psat_h2o(tl[layer]); break;
        case 9:  saturation_pressure = ms_psat_nh3(tl[layer]); break;
        case 20: saturation_pressure = ms_psat_co(tl[layer]); break;
        case 21: saturation_pressure = ms_psat_ch4(tl[layer]); break;
        case 52: saturation_pressure = ms_psat_co2(tl[layer]); break;
        case 53: saturation_pressure = ms_psat_h2(tl[layer]); break;
        case 54: saturation_pressure = ms_psat_o2(tl[layer]); break;
        case 55: saturation_pressure = ms_psat_n2(tl[layer]); break;
        case 45: saturation_pressure = ms_psat_h2s(tl[layer]); break;
        default: saturation_pressure = 1.0e20; break;
    }
    
    double delta_pressure = partial_pressure - saturation_pressure;
    
    if (delta_pressure <= 0.0 || saturation_pressure >= 1.0e15) {
        // Undersaturated or above critical point
        return 0.0;
    }
    
    // Get molecular mass and particle density
    double molecular_mass = get_molecular_mass(species_id); // g/mol
    double particle_density = get_particle_density(species_id, tl[layer]); // kg/m³
    
    // Calculate derived parameters (following particlesizef.c)
    double scale_height = KB * tl[layer] / (molecular_mass * AMU * GRAVITY);
    double eddy_velocity = atm->eddy_diffusion / scale_height;
    double delta_n = delta_pressure / (KB * tl[layer]);
    
    // Iterative solution for particle size (following particlesizef.c algorithm)
    double Cc = 1.0;  // Cunningham correction factor
    double fa = 1.0;  // Accommodation factor correction
    double sigma = CLOUD_SIGMA_DEFAULT;
    
    double particle_volume = 0.0;
    
    for (int iter = 0; iter < 1000; iter++) {
        // Calculate coefficients
        double cc = -pow(48.0 * PI * PI, 1.0/3.0) * MASS_DIFFUSION_COEFF * 
                   molecular_mass * AMU * fa * delta_n / particle_density * 
                   exp(-pow(log(sigma), 2.0));
        
        double aa = particle_density * GRAVITY / atm->viscosity / 
                   pow(162.0 * PI * PI, 1.0/3.0) / scale_height * Cc * 
                   exp(-pow(log(sigma), 2.0));
        
        double bb = -eddy_velocity / scale_height;
        
        // Solve quadratic equation
        double discriminant = bb * bb - 4.0 * aa * cc;
        if (discriminant < 0.0) break;
        
        double V = pow((-bb + sqrt(discriminant)) / (2.0 * aa), 3.0/2.0);
        double diameter = pow(6.0 * V / PI, 1.0/3.0) * exp(-pow(log(sigma), 2.0));
        
        // Calculate Knudsen number and corrections
        double kn = atm->mean_free_path / diameter;
        double Cc_new = 1.0 + kn * (1.257 + 0.4 * exp(-1.1 / kn));
        double fa_new = (1.0 + kn) / (1.0 + 2.0 * kn * (1.0 + kn) / ACCOMMODATION_COEFF);
        
        // Check convergence
        if (fabs(Cc_new - Cc) + fabs(fa_new - fa) < 0.001) {
            particle_volume = V;
            break;
        }
        
        Cc = Cc_new;
        fa = fa_new;
    }
    
    // Calculate effective radius
    double radius = pow(3.0 * particle_volume / (4.0 * PI), 1.0/3.0) * 
                   exp(-0.5 * pow(log(sigma), 2.0));
    
    // Apply bounds
    radius = fmax(MIN_PARTICLE_RADIUS, fmin(radius, MAX_PARTICLE_RADIUS));
    
    return radius;
}

double calculate_particle_radius_simple(int species_id, double condensate_mass_density, double particle_density) {
    /**
     * Simple particle size calculation for fallback
     * Assumes typical particle number density
     */
    
    if (condensate_mass_density <= 0.0 || particle_density <= 0.0) {
        return 0.0;
    }
    
    // Assume typical particle number density for exoplanet clouds
    double typical_number_density = 1.0e8; // particles/m³
    
    // Calculate particle mass
    double single_particle_mass = condensate_mass_density / typical_number_density;
    
    // Calculate radius from mass and density
    double radius = pow(3.0 * single_particle_mass / (4.0 * PI * particle_density), 1.0/3.0);
    
    // Apply bounds
    radius = fmax(MIN_PARTICLE_RADIUS, fmin(radius, MAX_PARTICLE_RADIUS));
    
    return radius;
}

//=========================================================
//=== Sedimentation Velocity Functions ===================
//=========================================================

double calculate_fall_velocity(double particle_radius, double particle_density, AtmosphericLayer* atm_layer) {
    /**
     * Calculate particle fall velocity using Stokes law with corrections
     * Handles both Stokes and Epstein drag regimes
     */
    
    if (particle_radius <= 0.0 || particle_density <= 0.0) {
        return 0.0;
    }
    
    // Calculate Knudsen number
    double knudsen = atm_layer->mean_free_path / (2.0 * particle_radius);
    
    // Calculate basic Stokes velocity
    double stokes_velocity = (2.0 * GRAVITY * particle_radius * particle_radius * 
                             particle_density) / (9.0 * atm_layer->viscosity);
    
    // Apply Cunningham correction for small particles
    double cunningham_correction = calculate_cunningham_correction(particle_radius, atm_layer->mean_free_path);
    
    // Calculate fall velocity
    double fall_velocity = stokes_velocity * cunningham_correction;
    
    // For very small particles (Epstein regime), use molecular kinetic theory
    if (knudsen > 1.0) {
        double thermal_velocity = sqrt(8.0 * KB * atm_layer->temperature / (PI * 2.33 * AMU)); // Assume H2-dominated
        double epstein_velocity = (4.0 * particle_radius * particle_density * GRAVITY) / 
                                 (3.0 * atm_layer->density * thermal_velocity);
        
        // Blend between regimes
        double blend_factor = 1.0 / (1.0 + knudsen);
        fall_velocity = blend_factor * fall_velocity + (1.0 - blend_factor) * epstein_velocity;
    }
    
    // Apply bounds
    fall_velocity = fmax(MIN_SEDIMENTATION_VELOCITY, fmin(fall_velocity, MAX_SEDIMENTATION_VELOCITY));
    
    return fall_velocity;
}

double calculate_cunningham_correction(double particle_radius, double mean_free_path) {
    /**
     * Calculate Cunningham slip correction factor
     */
    double knudsen = mean_free_path / (2.0 * particle_radius);
    return 1.0 + knudsen * (1.257 + 0.4 * exp(-1.1 / knudsen));
}

//=========================================================
//=== Atmospheric Property Functions =====================
//=========================================================

void calculate_atmospheric_properties(int layer, AtmosphericLayer* atm_layer) {
    /**
     * Calculate atmospheric properties needed for cloud physics
     */
    
    atm_layer->temperature = tl[layer];
    atm_layer->pressure = pl[layer];
    
    // Calculate gas density
    atm_layer->density = (atm_layer->pressure * meanmolecular[layer] * AMU) / 
                        (KB * atm_layer->temperature);
    
    // Calculate mean free path
    atm_layer->mean_free_path = calculate_mean_free_path(atm_layer->temperature, 
                                                        atm_layer->pressure, 
                                                        meanmolecular[layer]);
    
    // Calculate dynamic viscosity
    atm_layer->viscosity = calculate_dynamic_viscosity(atm_layer->temperature);
    
    // Calculate scale height
    atm_layer->scale_height = (KB * atm_layer->temperature) / 
                             (meanmolecular[layer] * AMU * GRAVITY);
    
    // Calculate eddy diffusion coefficient
    atm_layer->eddy_diffusion = calculate_eddy_diffusion_coefficient(layer);
}

double calculate_mean_free_path(double temperature, double pressure, double molecular_mass) {
    /**
     * Calculate mean free path of gas molecules
     */
    double collision_cross_section = 2.0e-19; // m² (typical for H2)
    double number_density = pressure / (KB * temperature);
    return 1.0 / (sqrt(2.0) * collision_cross_section * number_density);
}

double calculate_dynamic_viscosity(double temperature) {
    /**
     * Calculate dynamic viscosity using Sutherland's law
     * Assumes H2-dominated atmosphere
     */
    double T_ref = 273.15; // K
    double mu_ref = 8.76e-6; // Pa·s
    double S = 72.0; // K (Sutherland constant for H2)
    
    return mu_ref * pow(temperature / T_ref, 1.5) * (T_ref + S) / (temperature + S);
}

double calculate_eddy_diffusion_coefficient(int layer) {
    /**
     * Calculate eddy diffusion coefficient
     * Can be enhanced with more sophisticated models later
     */
    
    // Simple parameterization based on pressure
    double pressure_bar = pl[layer] / 1.0e5;
    double K_zz = EDDY_DIFFUSION_COEFF;
    
    // Scale with pressure (typical for exoplanet atmospheres)
    if (pressure_bar < 0.1) {
        K_zz *= pow(pressure_bar / 0.1, -0.5); // Increase at low pressure
    }
    
    return K_zz;
}

//=========================================================
//=== Mass Conservation Functions ========================
//=========================================================

void enforce_mass_conservation(int layer) {
    /**
     * Ensure mass conservation after cloud distribution calculations
     * This is critical for maintaining physical consistency
     */
    
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        // Get current vapor and cloud amounts
        double vapor_amount = xx[layer][species_id];
        double cloud_amount = clouds[layer][species_id];
        double total_amount = vapor_amount + cloud_amount;
        
        // Check for negative values
        if (vapor_amount < 0.0) {
            printf("WARNING: Negative vapor amount for species %d in layer %d: %.6e\n", 
                   species_id, layer, vapor_amount);
            xx[layer][species_id] = 0.0;
        }
        
        if (cloud_amount < 0.0) {
            printf("WARNING: Negative cloud amount for species %d in layer %d: %.6e\n", 
                   species_id, layer, cloud_amount);
            clouds[layer][species_id] = 0.0;
        }
        
        // Ensure total abundance doesn't exceed physical limits
        double max_abundance = MM[layer] * 0.1; // Maximum 10% by number
        if (total_amount > max_abundance) {
            double scale_factor = max_abundance / total_amount;
            xx[layer][species_id] *= scale_factor;
            clouds[layer][species_id] *= scale_factor;
            
            printf("WARNING: Scaled down abundances for species %d in layer %d by factor %.3f\n", 
                   species_id, layer, scale_factor);
        }
    }
}

//=========================================================
//=== Utility Functions ==================================
//=========================================================

double get_particle_density(int species_id, double temperature) {
    /**
     * Get particle density based on species and temperature
     * Returns density in kg/m³
     */
    
    switch(species_id) {
        case 7:  // H2O
            return (temperature < 273.15) ? PARTICLE_DENSITY_H2O_ICE : PARTICLE_DENSITY_H2O_LIQ;
        case 9:  // NH3
            return PARTICLE_DENSITY_NH3_ICE;
        case 20: // CO
            return PARTICLE_DENSITY_CO_ICE;
        case 21: // CH4
            return PARTICLE_DENSITY_CH4_ICE;
        case 52: // CO2
            return PARTICLE_DENSITY_CO2_ICE;
        case 53: // H2
            return PARTICLE_DENSITY_H2_ICE;
        case 54: // O2
            return PARTICLE_DENSITY_O2_ICE;
        case 55: // N2
            return PARTICLE_DENSITY_N2_ICE;
        case 45: // H2S
            return PARTICLE_DENSITY_H2S_ICE;
        default:
            return 1000.0; // Default density
    }
}

double get_molecular_mass(int species_id) {
    /**
     * Get molecular mass in g/mol
     */
    
    switch(species_id) {
        case 7:  return 18.015; // H2O
        case 9:  return 17.031; // NH3
        case 20: return 28.010; // CO
        case 21: return 16.043; // CH4
        case 52: return 44.010; // CO2
        case 53: return 2.016;  // H2
        case 54: return 31.998; // O2
        case 55: return 28.014; // N2
        case 45: return 34.081; // H2S
        default: return 18.015; // Default to H2O
    }
}

void calculate_optical_properties(int layer) {
    /**
     * Calculate optical properties for radiative transfer
     * This will interface with your existing Mie scattering code
     */
    
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        if (particle_radius[layer][i] <= 0.0) {
            // No particles
            cloud_optical_depth[layer][i] = 0.0;
            cloud_single_scattering_albedo[layer][i] = 1.0;
            cloud_asymmetry_parameter[layer][i] = 0.0;
            continue;
        }
        
        // Calculate optical properties based on particle size
        double radius_microns = particle_radius[layer][i] * 1.0e6;
        
        // Simple parameterization for optical properties
        if (species_id == 7) { // H2O
            if (tl[layer] < 273.15) {
                // Ice
                cloud_single_scattering_albedo[layer][i] = 0.99;
                cloud_asymmetry_parameter[layer][i] = 0.85;
            } else {
                // Liquid
                cloud_single_scattering_albedo[layer][i] = 0.95;
                cloud_asymmetry_parameter[layer][i] = 0.80;
            }
        } else {
            // Other species - use defaults
            cloud_single_scattering_albedo[layer][i] = 0.90;
            cloud_asymmetry_parameter[layer][i] = 0.75;
        }
        
        // Calculate optical depth (simplified)
        double number_density = particle_number_density[layer][i];
        double cross_section = PI * particle_radius[layer][i] * particle_radius[layer][i];
        double layer_thickness = atm_properties[layer].scale_height * 
                               log(P[layer-1] / P[layer]);
        
        cloud_optical_depth[layer][i] = number_density * cross_section * layer_thickness;
    }
}

void write_cloud_diagnostics(int iteration, const char* filename) {
    /**
     * Write cloud diagnostics to file for analysis
     */
    
    FILE* fp = fopen(filename, "w");
    if (!fp) return;
    
    fprintf(fp, "# Cloud Physics Diagnostics - Iteration %d\n", iteration);
    fprintf(fp, "# Layer  Species  Radius(um)  Number_Density(m-3)  Fall_Velocity(m/s)  Cloud_VMR  Vapor_VMR\n");
    
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < NCONDENSIBLES; i++) {
            int species_id = CONDENSIBLES[i];
            
            fprintf(fp, "%d  %d  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                   layer, species_id,
                   particle_radius[layer][i] * 1.0e6, // Convert to microns
                   particle_number_density[layer][i],
                   sedimentation_velocity[layer][i],
                   clouds[layer][species_id] / MM[layer],
                   xx[layer][species_id] / MM[layer]);
        }
    }
    
    fclose(fp);
} 