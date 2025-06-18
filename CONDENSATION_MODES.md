# Dynamic Condensation Detection in EPACRIS

This document explains the new dynamic condensation detection feature that allows automatic identification of condensible species based on atmospheric conditions.

## Overview

The EPACRIS model now supports three modes for determining which species should be treated as condensible (Xv) versus dry (Xd):

1. **Manual Mode (CONDENSATION_MODE = 0)**: Use predefined lists (original behavior)
2. **Automatic Mode (CONDENSATION_MODE = 1)**: Dynamically detect condensibles based on saturation
3. **Hybrid Mode (CONDENSATION_MODE = 2)**: Start with manual list, validate and expand automatically

## Configuration

### In `Input/conv_test/AlphaAb.h`:

```c
#define CONDENSATION_MODE 1  // Choose your mode (0, 1, or 2)

// Timing control for when detection occurs
#define CONDENSATION_TIMING 1  // 0 = Once before RC loop
                              // 1 = Every RC iteration (dynamic)
                              // 2 = Before loop + every N iterations
#define DETECTION_FREQUENCY 3  // For TIMING=2: detect every N iterations

// Manual mode settings
#define NCONDENSIBLES_MANUAL 2
#define CONDENSIBLES_MANUAL (int[]){7,9}  // H2O, NH3
#define ALPHA_RAINOUT_MANUAL (double[]){0.1,0.1}

// Automatic detection parameters
#define SATURATION_THRESHOLD 0.01  // 1% saturation ratio threshold
#define TEMP_RANGE_CHECK 1         // Check temperature ranges
#define MAX_CONDENSIBLES 20        // Maximum species that can condense
```

## Detection Timing Options

### CONDENSATION_TIMING = 0 (Static Detection)
- **When**: Once before the radiative-convective loop starts
- **Pros**: Fast, stable, no convergence issues
- **Cons**: Misses species that become condensible as temperature profile evolves
- **Best for**: Initial testing, known atmospheric compositions

### CONDENSATION_TIMING = 1 (Dynamic Detection)
- **When**: Every radiative-convective iteration
- **Pros**: Most physically accurate, adapts to changing temperature profile
- **Cons**: Computational overhead, potential convergence issues if list changes frequently
- **Best for**: Exploratory studies, atmospheres with large temperature changes

### CONDENSATION_TIMING = 2 (Periodic Detection)
- **When**: Before loop + every N iterations (set by DETECTION_FREQUENCY)
- **Pros**: Balance between accuracy and stability
- **Cons**: May miss rapid changes between detection points
- **Best for**: Most practical applications, good compromise

## How It Works

### Manual Mode (0)
- Uses the predefined `CONDENSIBLES_MANUAL` list
- Behavior identical to original EPACRIS
- Best for: Known atmospheric compositions, validation runs

### Automatic Mode (1)
- Scans all atmospheric layers for species approaching saturation
- Checks species: H2O(7), NH3(9), CO(20), CH4(21), CO2(52), H2(53), O2(54), N2(55)
- Includes species with saturation ratio > `SATURATION_THRESHOLD`
- Sets default alpha values (10% retention)
- Best for: Exploratory studies, unknown compositions

### Hybrid Mode (2)
- Starts with manual list as baseline
- Adds any additional species detected automatically
- Preserves manual alpha values, uses defaults for new species
- Best for: Extending known compositions, comprehensive studies

## Detection Criteria

A species is considered condensible if:

1. **Saturation pressure is finite**: Not the 1e+20 "above critical" case
2. **Saturation ratio > threshold**: `partial_pressure/saturation_pressure > SATURATION_THRESHOLD`
3. **Temperature in range**: Within condensation temperature range (if `TEMP_RANGE_CHECK = 1`)

## Species Temperature Ranges

| Species | ID | Solid Phase | Liquid Phase | Critical Point |
|---------|----|-----------|--------------|--------------| 
| H2O     | 7  | < 273K    | 273-647K     | > 647K       |
| NH3     | 9  | < 195K    | 195-405K     | > 405K       |
| CO      | 20 | < 68K     | 68-134K      | > 134K       |
| CH4     | 21 | < 91K     | 91-190K      | > 190K       |
| CO2     | 52 | < 217K    | 217-304K     | > 304K       |
| H2      | 53 | < 14K     | 14-33K       | > 33K        |
| O2      | 54 | < 54K     | 54-155K      | > 155K       |
| N2      | 55 | < 63K     | 63-126K      | > 126K       |

## Output

The model will print diagnostic information:

```
CONDENSATION MODE: Automatic - Will detect condensibles dynamically
DETECTING CONDENSIBLE SPECIES ACROSS ATMOSPHERE...
  Species 7 detected as condensible
  Species 52 detected as condensible
FINAL CONDENSIBLES LIST: 2 species
Species IDs: 7 52 
Alpha values: 0.10 0.10 
```

## Usage Examples

### For Hot Jupiter (T > 1000K):
- Automatic mode might detect: None (all above critical points)
- Manual mode recommended with high-T species

### For Temperate Planet (T ~ 200-400K):
- Automatic mode might detect: H2O, CO2, NH3
- Hybrid mode recommended to ensure H2O is included

### For Cold Planet (T < 100K):
- Automatic mode might detect: Most species
- Manual mode recommended to avoid over-condensation

## Functions Added

### In `ms_conv_funcs.c`:

- `initialize_condensibles_mode()`: Sets up condensibles based on mode
- `detect_condensibles_atmosphere()`: Scans atmosphere for condensibles  
- `check_species_condensible()`: Tests individual species
- `update_condensibles_list()`: Updates list for specific layer

## Integration

The detection runs once at the start of the climate calculation, after the helium calculation but before the main radiative-convective loop. This ensures the condensibles list is established before any adiabatic calculations.

## Recommendations

1. **Start with Manual Mode** for validation against known results
2. **Use Automatic Mode** for exploratory studies of new atmospheres
3. **Use Hybrid Mode** when you want to ensure specific species are included but also discover others
4. **Adjust SATURATION_THRESHOLD** based on your precision requirements (0.01 = 1% of saturation)
5. **Monitor the diagnostic output** to understand what species are being detected

## Troubleshooting

- **No species detected**: Lower `SATURATION_THRESHOLD` or check temperature ranges
- **Too many species detected**: Raise `SATURATION_THRESHOLD` or use manual mode
- **Unexpected species**: Check partial pressures and saturation curves for those species
- **Performance issues**: Automatic detection adds minimal overhead, but manual mode is fastest 