/*Function to read in cross section data */
/* input cross section in m^2 */
/* input wavelength in microns */
/* output mean cross section in cm^2 */
/* Interpolation to specified T and P */

#include <math.h>
#include <time.h>  // Added for timing
#include <string.h>
#include "constant.h"

#define MAX_CACHED_SPECIES 50  // Maximum number of species to cache

// External references to variables defined in the main file
extern const char* species[];  // External reference to species array
extern double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO;
extern double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4;
extern double **opacNH3;
extern double **opacC2H2, **opacC2H4, **opacC2H6, **opacHCN, **opacCH2O2, **opacHNO3;
extern double **opacN2O, **opacN2, **opacNO, **opacNO2, **opacOCS;
extern double tl[], pl[];  // External reference to temperature and pressure arrays
// zbin is defined in a header as a macro, not a variable
// Use CROSSHEADING_STR which is the string version of the CROSSHEADING macro
extern char CROSSHEADING_STR[];  // External reference to opacity file path

// Define NUM_SPECIES macro
#define NUM_SPECIES (sizeof(species) / sizeof(species[0]))

void readcross(char Fname[], double **xsc);
void reinterpolate_opacities_by_file(double **xsc, char Fname[]);
void read_all_opacities();
void cleanup_opacity_cache(void);
void reinterpolate_all_opacities();


// Structure to hold opacity data for one species
struct OpacityData {
	char filename[1024];
	char species_name[10];  // Add species name field
	double ***opac;
	double wave[NLAMBDA];
	double temp[NTEMP];
	double pres[NPRESSURE];
	int is_loaded;
};

// Static storage for opacity data - one cache per species
static struct OpacityData opacity_caches[MAX_CACHED_SPECIES] = {0};
static int num_cached_species = 0;

// Helper function to find or create a cache entry for a species
static struct OpacityData* get_cache_entry(const char* filename) {
	// Extract species name from filename
	char species_name[10];
	const char* last_slash = strrchr(filename, '/');
	const char* opac_prefix = strstr(last_slash ? last_slash + 1 : filename, "opac");
	if (opac_prefix) {
		strncpy(species_name, opac_prefix + 4, 9);  // Skip "opac" prefix
		species_name[9] = '\0';  // Ensure null termination
	} else {
		strcpy(species_name, "unknown");
	}
	
	// First try to find existing cache by species name
	for (int i = 0; i < num_cached_species; i++) {
		if (strcmp(opacity_caches[i].species_name, species_name) == 0) {
			return &opacity_caches[i];
		}
	}
	
	// If not found and we have room, create new cache
	if (num_cached_species < MAX_CACHED_SPECIES) {
		struct OpacityData* new_cache = &opacity_caches[num_cached_species++];
		strcpy(new_cache->filename, filename);
		strcpy(new_cache->species_name, species_name);
		new_cache->opac = NULL;
		new_cache->is_loaded = 0;
		return new_cache;
	}
	
	// If we're out of cache slots, reuse the first one
	printf("Warning: Reusing first opacity cache slot\n");
	struct OpacityData* cache = &opacity_caches[0];
	if (cache->opac != NULL) {
		free_f3tensor(cache->opac, 0, NLAMBDA-1, 0, NTEMP-1, 0, NPRESSURE-1);
	}
	strcpy(cache->filename, filename);
	strcpy(cache->species_name, species_name);
	cache->is_loaded = 0;
	return cache;
}

// Debug function to print some sample opacity values
void print_debug_opacities(double **xsc, const char *stage, const char* filename) {
	struct OpacityData* cache = get_cache_entry(filename);
	printf("\nDEBUG %s for %s:\n", stage, filename);
	printf("Sample opacity values at z=1, lambda indices 0, NLAMBDA/2, NLAMBDA-1:\n");
	printf("lambda[0] = %.3e, xsc[1][0] = %.3e\n", cache->wave[0], xsc[1][0]);
	printf("lambda[%d] = %.3e, xsc[1][%d] = %.3e\n", 
		   NLAMBDA/2, cache->wave[NLAMBDA/2], 
		   NLAMBDA/2, xsc[1][NLAMBDA/2]);
	printf("lambda[%d] = %.3e, xsc[1][%d] = %.3e\n", 
		   NLAMBDA-1, cache->wave[NLAMBDA-1], 
		   NLAMBDA-1, xsc[1][NLAMBDA-1]);
	
	// Print some temperature and pressure values
	printf("\nTemperature range: %.1f to %.1f K\n", 
		   cache->temp[0], cache->temp[NTEMP-1]);
	printf("Pressure range: %.1e to %.1e Pa\n", 
		   cache->pres[0], cache->pres[NPRESSURE-1]);
	printf("Current T,P at z=1: %.1f K, %.1e Pa\n", tl[1], pl[1]);
}


// Read all opacities from the species list
void read_all_opacities() {
    char crossfile[1024];

    printf("Reading opacities for %d species\n", (int)NUM_SPECIES); 
    for (int i = 0; i < NUM_SPECIES; i++) {
        //printf("Reading opac%s.dat\n", species[i]);
        strcpy(crossfile, CROSSHEADING_STR);
        strcat(crossfile, "opac");
        strcat(crossfile, species[i]);
        strcat(crossfile, ".dat");
        
        // Use switch to map species name to correct variables
        switch(species[i][0]) {
            case 'C':
                if (strcmp(species[i], "C2H2") == 0) {
                    readcross(crossfile, opacC2H2);
                    // planckmean(MeanC2H2, SMeanC2H2, opacC2H2);
                } else if (strcmp(species[i], "C2H4") == 0) {
                    readcross(crossfile, opacC2H4);
                    // planckmean(MeanC2H4, SMeanC2H4, opacC2H4);
                } else if (strcmp(species[i], "C2H6") == 0) {
                    readcross(crossfile, opacC2H6);
                    // planckmean(MeanC2H6, SMeanC2H6, opacC2H6);
                } else if (strcmp(species[i], "CH2O2") == 0) {
                    readcross(crossfile, opacCH2O2);
                    // planckmean(MeanCH2O2, SMeanCH2O2, opacCH2O2);
                } else if (strcmp(species[i], "CH4") == 0) {
                    readcross(crossfile, opacCH4);
                    // planckmean(MeanCH4, SMeanCH4, opacCH4);
                } else if (strcmp(species[i], "CO2") == 0) {
                    readcross(crossfile, opacCO2);
                    // planckmean(MeanCO2, SMeanCO2, opacCO2);
                } else if (strcmp(species[i], "CO") == 0) {
                    readcross(crossfile, opacCO);
                    // planckmean(MeanCO, SMeanCO, opacCO);
                }
                break;
            case 'H':
                if (strcmp(species[i], "H2CO") == 0) {
                    readcross(crossfile, opacH2CO);
                    // planckmean(MeanH2CO, SMeanH2CO, opacH2CO);
                } else if (strcmp(species[i], "H2O2") == 0) {
                    readcross(crossfile, opacH2O2);
                    // planckmean(MeanH2O2, SMeanH2O2, opacH2O2);
                } else if (strcmp(species[i], "H2O") == 0) {
                    readcross(crossfile, opacH2O);
                    // planckmean(MeanH2O, SMeanH2O, opacH2O);
                } else if (strcmp(species[i], "H2S") == 0) {
                    readcross(crossfile, opacH2S);
                    // planckmean(MeanH2S, SMeanH2S, opacH2S);
                } else if (strcmp(species[i], "HCN") == 0) {
                    readcross(crossfile, opacHCN);
                    // planckmean(MeanHCN, SMeanHCN, opacHCN);
                } else if (strcmp(species[i], "HNO3") == 0) {
                    readcross(crossfile, opacHNO3);
                    // planckmean(MeanHNO3, SMeanHNO3, opacHNO3);
                } else if (strcmp(species[i], "HO2") == 0) {
                    readcross(crossfile, opacHO2);
                    // planckmean(MeanHO2, SMeanHO2, opacHO2);
                }
                break;
            case 'N':
                if (strcmp(species[i], "N2O") == 0) {
                    readcross(crossfile, opacN2O);
                    // planckmean(MeanN2O, SMeanN2O, opacN2O);
                } else if (strcmp(species[i], "N2") == 0) {
                    readcross(crossfile, opacN2);
                    // planckmean(MeanN2, SMeanN2, opacN2);
                } else if (strcmp(species[i], "NH3") == 0) {
                    readcross(crossfile, opacNH3);
                    // planckmean(MeanNH3, SMeanNH3, opacNH3);
                } else if (strcmp(species[i], "NO2") == 0) {
                    readcross(crossfile, opacNO2);
                    // planckmean(MeanNO2, SMeanNO2, opacNO2);
                } else if (strcmp(species[i], "NO") == 0) {
                    readcross(crossfile, opacNO);
                    // planckmean(MeanNO, SMeanNO, opacNO);
                }
                break;
            case 'O':
                if (strcmp(species[i], "O2") == 0) {
                    readcross(crossfile, opacO2);
                    // planckmean(MeanO2, SMeanO2, opacO2);
                } else if (strcmp(species[i], "O3") == 0) {
                    readcross(crossfile, opacO3);
                    // planckmean(MeanO3, SMeanO3, opacO3);
                } else if (strcmp(species[i], "OCS") == 0) {
                    readcross(crossfile, opacOCS);
                    // planckmean(MeanOCS, SMeanOCS, opacOCS);
                } else if (strcmp(species[i], "OH") == 0) {
                    readcross(crossfile, opacOH);
                    // planckmean(MeanOH, SMeanOH, opacOH);
                }
                break;
            case 'S':
                if (strcmp(species[i], "SO2") == 0) {
                    readcross(crossfile, opacSO2);
                    // planckmean(MeanSO2, SMeanSO2, opacSO2);
                }
                break;
        }
    }
}

// Function to clean up opacity cache
void cleanup_opacity_cache(void) {
    for (int i = 0; i < num_cached_species; i++) {
        if (opacity_caches[i].is_loaded && opacity_caches[i].opac != NULL) {
            free_f3tensor(opacity_caches[i].opac, 0, NLAMBDA-1, 0, NTEMP-1, 0, NPRESSURE-1);
            opacity_caches[i].is_loaded = 0;
            opacity_caches[i].opac = NULL;
        }
    }
    num_cached_species = 0;
}

// Function to reinterpolate all opacities
void reinterpolate_all_opacities() {
    char crossfile[1024];
    
    //printf("Reinterpolating opacities for %d species\n", (int)NUM_SPECIES); 
    for (int i = 0; i < NUM_SPECIES; i++) {
        //printf("Reinterpolating %s opacities\n", species[i]);
        
        // Build the filename
        strcpy(crossfile, CROSSHEADING_STR);
        strcat(crossfile, "opac");
        strcat(crossfile, species[i]);
        strcat(crossfile, ".dat");
        
        // Use switch to map species name to correct variables
        switch(species[i][0]) {
            case 'C':
                if (strcmp(species[i], "C2H2") == 0) {
                    reinterpolate_opacities_by_file(opacC2H2, crossfile);
                } else if (strcmp(species[i], "C2H4") == 0) {
                    reinterpolate_opacities_by_file(opacC2H4, crossfile);
                } else if (strcmp(species[i], "C2H6") == 0) {
                    reinterpolate_opacities_by_file(opacC2H6, crossfile);
                } else if (strcmp(species[i], "CH2O2") == 0) {
                    reinterpolate_opacities_by_file(opacCH2O2, crossfile);
                } else if (strcmp(species[i], "CH4") == 0) {
                    reinterpolate_opacities_by_file(opacCH4, crossfile);
                } else if (strcmp(species[i], "CO2") == 0) {
                    reinterpolate_opacities_by_file(opacCO2, crossfile);
                } else if (strcmp(species[i], "CO") == 0) {
                    reinterpolate_opacities_by_file(opacCO, crossfile);
                }
                break;
            case 'H':
                if (strcmp(species[i], "H2CO") == 0) {
                    reinterpolate_opacities_by_file(opacH2CO, crossfile);
                } else if (strcmp(species[i], "H2O2") == 0) {
                    reinterpolate_opacities_by_file(opacH2O2, crossfile);
                } else if (strcmp(species[i], "H2O") == 0) {
                    reinterpolate_opacities_by_file(opacH2O, crossfile);
                } else if (strcmp(species[i], "H2S") == 0) {
                    reinterpolate_opacities_by_file(opacH2S, crossfile);
                } else if (strcmp(species[i], "HCN") == 0) {
                    reinterpolate_opacities_by_file(opacHCN, crossfile);
                } else if (strcmp(species[i], "HNO3") == 0) {
                    reinterpolate_opacities_by_file(opacHNO3, crossfile);
                } else if (strcmp(species[i], "HO2") == 0) {
                    reinterpolate_opacities_by_file(opacHO2, crossfile);
                }
                break;
            case 'N':
                if (strcmp(species[i], "N2O") == 0) {
                    reinterpolate_opacities_by_file(opacN2O, crossfile);
                } else if (strcmp(species[i], "N2") == 0) {
                    reinterpolate_opacities_by_file(opacN2, crossfile);
                } else if (strcmp(species[i], "NH3") == 0) {
                    reinterpolate_opacities_by_file(opacNH3, crossfile);
                } else if (strcmp(species[i], "NO2") == 0) {
                    reinterpolate_opacities_by_file(opacNO2, crossfile);
                } else if (strcmp(species[i], "NO") == 0) {
                    reinterpolate_opacities_by_file(opacNO, crossfile);
                }
                break;
            case 'O':
                if (strcmp(species[i], "O2") == 0) {
                    reinterpolate_opacities_by_file(opacO2, crossfile);
                } else if (strcmp(species[i], "O3") == 0) {
                    reinterpolate_opacities_by_file(opacO3, crossfile);
                } else if (strcmp(species[i], "OCS") == 0) {
                    reinterpolate_opacities_by_file(opacOCS, crossfile);
                } else if (strcmp(species[i], "OH") == 0) {
                    reinterpolate_opacities_by_file(opacOH, crossfile);
                }
                break;
            case 'S':
                if (strcmp(species[i], "SO2") == 0) {
                    reinterpolate_opacities_by_file(opacSO2, crossfile);
                }
                break;
        }
    }
}


void readcross(char Fname[], double **xsc)
{
	int i, j, k;
	FILE *fim;
	double wave1[NLAMBDA]; /* in nm */
	double presdummy, cross[NLAMBDA];
	double pfitting, tfitting;
	
	// Get the cache entry for this species
	struct OpacityData* cache = get_cache_entry(Fname);

	// Check if we need to load the data
	if (!cache->is_loaded) {
		printf("\nReading opacity file: %s\n", Fname);
		
		// Allocate opacity tensor if first time
		if (cache->opac == NULL) {
			cache->opac = f3tensor(0, NLAMBDA-1, 0, NTEMP-1, 0, NPRESSURE-1);
		}
		
		fim = fopen(Fname, "r");
		if (!fim) {
			printf("Error: Could not open opacity file %s\n", Fname);
			return;
		}
		
		/* Header Lines */
		for (i=0; i<NTEMP-1; i++) {
			fscanf(fim, "%lf", cache->temp+i);
		}
		i=NTEMP-1;
		fscanf(fim, "%lf\n", cache->temp+i);
		for (i=0; i<NPRESSURE-1; i++) {
			fscanf(fim, "%le", cache->pres+i);
		}
		i=NPRESSURE-1;
		fscanf(fim, "%le\n", cache->pres+i);
		/* Read in data */
		for (i=0; i<NLAMBDA; i++) {
			fscanf(fim, "%le\n", cache->wave+i);
			for (j=0; j<NPRESSURE; j++) {
				fscanf(fim, "%le", &presdummy);
				for (k=0; k<NTEMP-1; k++) {
					fscanf(fim, "%le", &cache->opac[i][k][j]);
				}
				k=NTEMP-1;
				fscanf(fim, "%le\n", &cache->opac[i][k][j]);
			}
		}
		fclose(fim);
		cache->is_loaded = 1;
		
		// Print debug info after reading
		//printf("\nSuccessfully read opacity file. Cache status:\n");
		//printf("is_loaded = %d\n", cache->is_loaded);
		//printf("opac pointer = %p\n", (void*)cache->opac);
	} else {
		printf("\nUsing cached opacity data from: %s\n", cache->filename);
	}
	
	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = cache->wave[i]*1.0E+9; /* convert to nm */
	}
	
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i] > cache->pres[NPRESSURE-1]) {
				pfitting = cache->pres[NPRESSURE-1];
			} else if (pl[i] < cache->pres[0]) {
				pfitting = cache->pres[0];
			} else {
				pfitting = pl[i];
			}
			if (tl[i] > cache->temp[NTEMP-1]) {
				tfitting = cache->temp[NTEMP-1];
			} else if (tl[i] < cache->temp[0]) {
				tfitting = cache->temp[0];
			} else {
				tfitting = tl[i];
			}
			cross[j] = fmax(Interpolation2D(tfitting, pfitting, 
										  cache->temp, NTEMP, 
										  cache->pres, NPRESSURE, 
										  *(cache->opac+j))*1.0E+4, 0.0);
			xsc[i][j] = cross[j];
		}
	}

	// Print debug info after interpolation
	//print_debug_opacities(xsc, "After initial interpolation", Fname);
}

// New function that takes a filename parameter
void reinterpolate_opacities_by_file(double **xsc, char Fname[]) {
	// Find the cache entry for this species
	struct OpacityData* cache = get_cache_entry(Fname);
	
	if (!cache || !cache->is_loaded) {
		printf("Error: No opacity data loaded for %s\n", Fname);
		return;
	}

    //debug
	// printf("\nReinterpolating opacities from cached data: %s (Species: %s)\n", 
	// 	   Fname, cache->species_name);

	int i, j;
	double wave1[NLAMBDA]; /* in nm */
	double pfitting, tfitting;
	double cross[NLAMBDA];

	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = cache->wave[i]*1.0E+9; /* convert to nm */
	}
	
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i] > cache->pres[NPRESSURE-1]) {
				pfitting = cache->pres[NPRESSURE-1];
			} else if (pl[i] < cache->pres[0]) {
				pfitting = cache->pres[0];
			} else {
				pfitting = pl[i];
			}
			if (tl[i] > cache->temp[NTEMP-1]) {
				tfitting = cache->temp[NTEMP-1];
			} else if (tl[i] < cache->temp[0]) {
				tfitting = cache->temp[0];
			} else {
				tfitting = tl[i];
			}
			cross[j] = fmax(Interpolation2D(tfitting, pfitting, 
										  cache->temp, NTEMP, 
										  cache->pres, NPRESSURE, 
										  *(cache->opac+j))*1.0E+4, 0.0);
			xsc[i][j] = cross[j];
		}
	}

	// Print debug info after reinterpolation
	//print_debug_opacities(xsc, "After reinterpolation", Fname);
}
