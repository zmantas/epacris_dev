#ifndef CLOUD_PHYSICS_SIMPLE_H
#define CLOUD_PHYSICS_SIMPLE_H

//=========================================================
//=== Function Prototypes ================================
//=========================================================

// Main enhanced cloud physics function
void enhanced_rainout_physics(int layer);

// Simple helper functions
double calculate_simple_particle_size(int species_id, double temperature, double pressure);
double calculate_simple_fall_velocity(int species_id, double radius_um, double temperature, double pressure);
double calculate_simple_retention(double fall_velocity, double temperature, double pressure);
void write_simple_diagnostics(int layer);

#endif // CLOUD_PHYSICS_SIMPLE_H 