# Cloud Optical Property Implementation

## Status: ✅ Data Reading Complete

This document summarizes the cloud optical property implementation for EPACRIS radiative transfer.

---

## What Has Been Implemented

### 1. **New Files Created**

- **`cloud_optics.h`** - Header file with data structures and function declarations
- **`cloud_optics.c`** - Implementation for reading Mie scattering lookup tables
- **`test_cloud_optics.c`** - Standalone test program to verify data reading

### 2. **Modified Files**

- **`global_temp.h`** - Added cloud optics function declarations
- **`epacris_test_my_v2.c`** - Added calls to read and cleanup cloud optical tables

### 3. **Data Structures**

```c
typedef struct {
    double radius[NPARTICLE_SIZES];           // 41 particle radii (micrometers)
    double cross[NPARTICLE_SIZES][NLAMBDA];   // Cross section (extinction)
    double albedo[NPARTICLE_SIZES][NLAMBDA];  // Single scattering albedo
    double asym[NPARTICLE_SIZES][NLAMBDA];    // Asymmetry parameter g
} CloudOpticalTable;
```

Three global tables are available:
- `H2O_ice_optics` - Water ice cloud properties
- `H2O_liquid_optics` - Water liquid cloud properties (currently copy of ice)
- `NH3_ice_optics` - Ammonia ice cloud properties

---

## Data Files Read

From `Clouds/CrossP/` directory:

### Water Clouds
- `Cross_water_wavelength2.dat` - Extinction cross sections
- `Albedo_water_wavelength2.dat` - Single scattering albedo
- `Geo_water_wavelength2.dat` - Asymmetry parameter (g)

### Ammonia Clouds
- `Cross_ammonia_wavelength2.dat` - Extinction cross sections
- `Albedo_ammonia_wavelength2.dat` - Single scattering albedo
- `Geo_ammonia_wavelength2.dat` - Asymmetry parameter (g)

**Format**: 41 rows (particle sizes) × 196 columns (wavelength bins)

---

## Current Integration

### Reading Phase (Initialization)

In `epacris_test_my_v2.c` around line 910:
```c
// Read cloud optical property lookup tables
printf("Reading cloud optical properties\n");
read_cloud_optical_tables();
```

### Cleanup Phase (End of Program)

Around line 1072:
```c
cleanup_cloud_optical_tables();  // Clean up cloud optical tables
```

---

## What Still Needs to Be Done (Next Steps)

### **Step 2: Interpolation Functions**

Create functions to interpolate optical properties based on:
- Particle size from `ms_conv_funcs.c` (currently stored in `particle_radius_um[layer][species]`)
- Wavelength index from radiative transfer loop
- Species ID (H2O vs NH3, ice vs liquid)

**Example function signature:**
```c
void get_cloud_opacity(int layer, int wavelength_idx, int species_id,
                      double particle_size_um,
                      double *extinction, double *albedo, double *asymmetry);
```

### **Step 3: Integration with Radiative Transfer**

In `ms_radtrans_test.c`:

1. **Allocate arrays** (before wavelength loop):
```c
double **cH2O, **aH2O, **gH2O;   // Water: cross-section, albedo, g
double **cNH3, **aNH3, **gNH3;   // NH3: cross-section, albedo, g
// Allocate as [zbin+1][NLAMBDA]
```

2. **Call interpolation** (inside wavelength loop, for each layer):
```c
for (j=1; j<=zbin; j++) {
    // Get particle size from ms_conv_funcs.c
    double particle_size = particle_radius_um[j][species_idx];
    
    // Interpolate optical properties
    get_cloud_opacity(j, i, SPECIES_H2O, particle_size,
                     &cH2O[j][i], &aH2O[j][i], &gH2O[j][i]);
}
```

3. **Apply to optical depth calculation** (currently commented out around line 251-255):
```c
// Uncomment and fix:
wa[j] += cH2O[j1][i]*(1.0-aH2O[j1][i])/MM[j1];  // Cloud absorption
wa[j] += cNH3[j1][i]*(1.0-aNH3[j1][i])/MM[j1];
ws[j] += cH2O[j1][i]*aH2O[j1][i]/MM[j1];        // Cloud scattering
ws[j] += cNH3[j1][i]*aNH3[j1][i]/MM[j1];
```

4. **Apply to asymmetry factor** (around line 292-293):
```c
g[j] += cH2O[j1][i]*aH2O[j1][i]*gH2O[j1][i]/MM[j1];
g[j] += cNH3[j1][i]*aNH3[j1][i]*gNH3[j1][i]/MM[j1];
```

---

## Testing

### Verify Data Reading
```bash
cd /mnt/d/JPL_models/epacris_dev
gcc -o test_cloud_optics test_cloud_optics.c -lm
./test_cloud_optics
```

### Expected Output
```
CLOUD OPTICAL TABLES SUCCESSFULLY LOADED
Particle size range: 0.0100 - [max_size] micrometers
Number of wavelength bins: 196
Species loaded: H2O (ice/liquid), NH3 (ice)
```

---

## Data Flow Diagram

```
Initialization (epacris_test_my_v2.c):
  read_cloud_optical_tables()
    ↓
  Loads: H2O_ice_optics, H2O_liquid_optics, NH3_ice_optics
    ↓
  Stored in memory for entire runtime

Radiative Transfer (ms_radtrans_test.c):
  For each wavelength i, layer j:
    ↓
  particle_size = particle_radius_um[j][species_idx]  (from ms_conv_funcs.c)
    ↓
  get_cloud_opacity(j, i, species, particle_size, &c, &a, &g)
    ↓
  Interpolate from lookup tables
    ↓
  Apply to wa[], ws[], g[] arrays
    ↓
  Calculate optical depth τ[j]
    ↓
  2-stream radiative transfer solver
```

---

## Notes

- **Particle sizes** calculated in `ms_conv_funcs.c` during convection
- **Cloud masses** stored in `clouds[layer][species_id]` array
- Current implementation reads tables once at initialization - no runtime file I/O
- Lookup tables use log-spaced particle size grid for efficient interpolation
- Each optical property (cross, albedo, g) depends on both particle size AND wavelength

---

## Contact

For questions about this implementation, refer to:
- Original Mie code: `Mie_code/` directory
- Cloud interpolation reference: `Clouds/cloud_interp.m` (MATLAB version)
- Particle size calculation: `ms_conv_funcs.c` lines 2600-2620

