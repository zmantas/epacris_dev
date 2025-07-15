// Enhanced cloud physics implementation - SIMPLE VERSION
// This file is included directly in epacris_test_my_v2.c following EPACRIS single-include architecture

//=========================================================
//=== Simple Enhanced Cloud Physics ======================
//=========================================================

void enhanced_rainout_physics(int layer)
{
    // Only run if enhanced cloud physics is enabled
    // Hardcode constants to avoid include order issues
    int cloud_physics_mode = 0;  // 0=Graham+2021 only, 1=Ackerman-Marley, 2=Hybrid
    
    if (cloud_physics_mode == 0) {
        return; // Graham+2021 only mode - do nothing
    }
    
    // For now, just apply a simple size-dependent retention factor
    // This is much simpler than my previous overcomplicated approach!
    
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        // Get current cloud amount
        double cloud_vmr = clouds[layer][species_id] / MM[layer];
        
        if (cloud_vmr < 1.0e-20) {
            // No clouds - nothing to do
            particle_radius_um[layer][i] = 0.0;
            fall_velocity_ms[layer][i] = 0.0;
            cloud_retention[layer][i] = 1.0;
            continue;
        }
        
        // Simple particle size estimate based on temperature and pressure
        double particle_size_um = calculate_simple_particle_size(species_id, tl[layer], pl[layer]);
        particle_radius_um[layer][i] = particle_size_um;
        
        // Simple fall velocity estimate
        double fall_velocity = calculate_simple_fall_velocity(species_id, particle_size_um, tl[layer], pl[layer]);
        fall_velocity_ms[layer][i] = fall_velocity;
        
        // Simple retention factor based on fall velocity
        double retention = calculate_simple_retention(fall_velocity, tl[layer], pl[layer]);
        cloud_retention[layer][i] = retention;
        
        // Apply the retention factor to clouds
        clouds[layer][species_id] *= retention;
    }
    
    // Optional: Write diagnostics for one layer only
    if (layer == 60 && cloud_physics_mode > 0) {
        write_simple_diagnostics(layer);
    }
}

double calculate_simple_particle_size(int species_id, double temperature, double pressure)
{
    // Very simple particle size estimation
    // Larger particles at lower temperatures, smaller at higher pressures
    
    double base_size = 10.0; // microns - reasonable default
    
    if (species_id == 7) {       // H2O
        if (temperature < 273.0) {
            base_size = 50.0;    // Ice crystals - larger
        } else {
            base_size = 10.0;    // Liquid droplets - smaller
        }
    } else if (species_id == 9) { // NH3
        base_size = 20.0;        // NH3 ice
    } else if (species_id == 21) { // CH4
        base_size = 30.0;        // CH4 ice
    } else if (species_id == 52) { // CO2
        base_size = 100.0;       // CO2 ice (dry ice) - larger
    }
    
    // Simple temperature and pressure dependence
    double temp_factor = 300.0 / temperature; // Larger at lower T
    double pressure_factor = 1.0e5 / pressure; // Larger at lower P
    
    double final_size = base_size * sqrt(temp_factor) * sqrt(pressure_factor);
    
    // Apply reasonable bounds
    if (final_size < 1.0) final_size = 1.0;     // Minimum 1 micron
    if (final_size > 1000.0) final_size = 1000.0; // Maximum 1 mm
    
    return final_size;
}

double calculate_simple_fall_velocity(int species_id, double radius_um, double temperature, double pressure)
{
    // Simple Stokes law fall velocity
    double radius_m = radius_um * 1.0e-6; // Convert to meters
    
    // Simple particle density estimates (kg/m³)
    double particle_density = 1000.0; // Default
    if (species_id == 7) {       // H2O
        particle_density = (temperature < 273.0) ? 917.0 : 1000.0;
    } else if (species_id == 9) { // NH3
        particle_density = 817.0;
    } else if (species_id == 21) { // CH4
        particle_density = 424.0;
    } else if (species_id == 52) { // CO2
        particle_density = 1563.0;
    }
    
    // Simple atmospheric properties
    double gas_density = pressure / (287.0 * temperature); // Approximate air density
    double viscosity = 1.8e-5 * pow(temperature / 273.0, 0.7); // Simple viscosity
    
    // Stokes velocity
    double gravity = GRAVITY * MASS_PLANET / (RADIUS_PLANET * RADIUS_PLANET);
    double stokes_velocity = (2.0 * radius_m * radius_m * particle_density * gravity) / (9.0 * viscosity);
    
    // Apply reasonable bounds
    if (stokes_velocity < 1.0e-6) stokes_velocity = 1.0e-6; // 1 μm/s minimum
    if (stokes_velocity > 10.0) stokes_velocity = 10.0;     // 10 m/s maximum
    
    return stokes_velocity;
}

double calculate_simple_retention(double fall_velocity, double temperature, double pressure)
{
    // Simple retention calculation inspired by Ackerman & Marley
    
    // Estimate atmospheric scale height
    double scale_height = (KBOLTZMANN * temperature) / (29.0 * AMU * GRAVITY); // Assuming ~29 amu mean molecular weight
    
    // Simple eddy diffusion coefficient
    double eddy_diffusion = 1.0e6;  // m²/s (typical for exoplanets)
    
    // Calculate settling parameter
    double f_sed = (fall_velocity * scale_height) / eddy_diffusion;
    
    // Retention factor: what fraction stays suspended
    double retention = 1.0 / (1.0 + f_sed);
    
    // Apply bounds
    if (retention < 0.01) retention = 0.01; // Minimum 1% retention
    if (retention > 0.99) retention = 0.99; // Maximum 99% retention
    
    return retention;
}

void write_simple_diagnostics(int layer)
{
    FILE *fp = fopen("simple_cloud_diagnostics.dat", "w");
    if (fp == NULL) return;
    
    fprintf(fp, "# Simple Enhanced Cloud Physics Diagnostics - Layer %d\n", layer);
    fprintf(fp, "# Mode: %d, Temperature: %.1f K, Pressure: %.2e Pa\n", 
            0, tl[layer], pl[layer]);
    fprintf(fp, "# Species  Radius_um  Fall_vel_ms  Retention  Cloud_VMR\n");
    
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        double cloud_vmr = clouds[layer][species_id] / MM[layer];
        
        if (cloud_vmr > 1.0e-20) {
            fprintf(fp, "%d  %.2f  %.4e  %.3f  %.4e\n",
                   species_id,
                   particle_radius_um[layer][i],
                   fall_velocity_ms[layer][i],
                   cloud_retention[layer][i],
                   cloud_vmr);
        }
    }
    
    fclose(fp);
    printf("Enhanced cloud diagnostics written for layer %d\n", layer);
} 