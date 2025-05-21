/*Function to read in cross section data */
/* input cross section in m^2 */
/* input wavelength in microns */
/* output mean cross section in cm^2 */
/* Interpolation to specified T and P */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constant.h"

#define MAX_SPECIES_CACHE 50  // Maximum number of species to cache

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

// Function prototypes
void readcross(char Fname[], double **xsc);
void reinterpolate_opacities_by_file(double **xsc, char Fname[]);
void read_all_opacities();
void cleanup_opacity_cache(void);
void reinterpolate_all_opacities();

// Simple structure to hold opacity data for one species
typedef struct {
	char filename[256];          // File path
	char species_name[32];       // Species name (e.g., "H2O")
	double ***opac;              // 3D opacity data array [wavelength][temperature][pressure]
	double wave[NLAMBDA];        // Wavelength array
	double temp[NTEMP];          // Temperature array
	double pres[NPRESSURE];      // Pressure array
	int is_loaded;               // Whether data is loaded
} OpacityFile;

// Static array of all opacity files we support
static OpacityFile opacity_files[MAX_SPECIES_CACHE] = {0};
static int num_cached_files = 0;

// Helper function to get an opacity file by species name or filename
static OpacityFile* get_opacity_file(const char* filename) {
	// Extract species name from filename
	const char* last_slash = strrchr(filename, '/');
	const char* file_base = last_slash ? last_slash + 1 : filename;
	const char* opac_prefix = strstr(file_base, "opac");
	
	char species_name[32] = {0};
	if (opac_prefix) {
		strncpy(species_name, opac_prefix + 4, sizeof(species_name) - 1);
		// Remove .dat extension if present
		char* dot = strchr(species_name, '.');
		if (dot) *dot = '\0';
	} else {
		strcpy(species_name, "unknown");
	}
	
	// Check if we already have this file cached
	for (int i = 0; i < num_cached_files; i++) {
		if (strcmp(opacity_files[i].species_name, species_name) == 0) {
			return &opacity_files[i];
		}
	}
	
	// If not found and we have room, create a new entry
	if (num_cached_files < MAX_SPECIES_CACHE) {
		OpacityFile* new_file = &opacity_files[num_cached_files++];
		strcpy(new_file->filename, filename);
		strcpy(new_file->species_name, species_name);
		new_file->opac = NULL;
		new_file->is_loaded = 0;
		return new_file;
	}
	
	// If we're out of cache slots, return the first one (least recently used)
	printf("Warning: Opacity cache full, reusing first slot\n");
	return &opacity_files[0];
}

// Read all opacities from the species list
void read_all_opacities() {
	char crossfile[1024];

	printf("Reading opacities for %d species\n", (int)NUM_SPECIES); 
	for (int i = 0; i < NUM_SPECIES; i++) {
		strcpy(crossfile, CROSSHEADING_STR);
		strcat(crossfile, "opac");
		strcat(crossfile, species[i]);
		strcat(crossfile, ".dat");
		
		// Use switch to map species name to correct variables
		switch(species[i][0]) {
			case 'C':
				if (strcmp(species[i], "C2H2") == 0) {
					readcross(crossfile, opacC2H2);
				} else if (strcmp(species[i], "C2H4") == 0) {
					readcross(crossfile, opacC2H4);
				} else if (strcmp(species[i], "C2H6") == 0) {
					readcross(crossfile, opacC2H6);
				} else if (strcmp(species[i], "CH2O2") == 0) {
					readcross(crossfile, opacCH2O2);
				} else if (strcmp(species[i], "CH4") == 0) {
					readcross(crossfile, opacCH4);
				} else if (strcmp(species[i], "CO2") == 0) {
					readcross(crossfile, opacCO2);
				} else if (strcmp(species[i], "CO") == 0) {
					readcross(crossfile, opacCO);
				}
				break;
			case 'H':
				if (strcmp(species[i], "H2CO") == 0) {
					readcross(crossfile, opacH2CO);
				} else if (strcmp(species[i], "H2O2") == 0) {
					readcross(crossfile, opacH2O2);
				} else if (strcmp(species[i], "H2O") == 0) {
					readcross(crossfile, opacH2O);
				} else if (strcmp(species[i], "H2S") == 0) {
					readcross(crossfile, opacH2S);
				} else if (strcmp(species[i], "HCN") == 0) {
					readcross(crossfile, opacHCN);
				} else if (strcmp(species[i], "HNO3") == 0) {
					readcross(crossfile, opacHNO3);
				} else if (strcmp(species[i], "HO2") == 0) {
					readcross(crossfile, opacHO2);
				}
				break;
			case 'N':
				if (strcmp(species[i], "N2O") == 0) {
					readcross(crossfile, opacN2O);
				} else if (strcmp(species[i], "N2") == 0) {
					readcross(crossfile, opacN2);
				} else if (strcmp(species[i], "NH3") == 0) {
					readcross(crossfile, opacNH3);
				} else if (strcmp(species[i], "NO2") == 0) {
					readcross(crossfile, opacNO2);
				} else if (strcmp(species[i], "NO") == 0) {
					readcross(crossfile, opacNO);
				}
				break;
			case 'O':
				if (strcmp(species[i], "O2") == 0) {
					readcross(crossfile, opacO2);
				} else if (strcmp(species[i], "O3") == 0) {
					readcross(crossfile, opacO3);
				} else if (strcmp(species[i], "OCS") == 0) {
					readcross(crossfile, opacOCS);
				} else if (strcmp(species[i], "OH") == 0) {
					readcross(crossfile, opacOH);
				}
				break;
			case 'S':
				if (strcmp(species[i], "SO2") == 0) {
					readcross(crossfile, opacSO2);
				}
				break;
		}
	}
}

// Function to clean up opacity cache
void cleanup_opacity_cache(void) {
	printf("Cleaning up opacity cache...\n");
	
	for (int i = 0; i < num_cached_files; i++) {
		if (opacity_files[i].is_loaded && opacity_files[i].opac != NULL) {
			free_f3tensor(opacity_files[i].opac, 0, NLAMBDA-1, 0, NTEMP-1, 0, NPRESSURE-1);
			opacity_files[i].opac = NULL;
			opacity_files[i].is_loaded = 0;
			printf("Freed opacity cache for %s\n", opacity_files[i].species_name);
		}
	}
	
	num_cached_files = 0;
	printf("Opacity cache cleanup complete\n");
}

// Function to reinterpolate all opacities
void reinterpolate_all_opacities() {
	char crossfile[1024];
	
	for (int i = 0; i < NUM_SPECIES; i++) {
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

// Main function to read opacity data for a specific file
void readcross(char Fname[], double **xsc)
{
	int i, j, k;
	FILE *fim;
	double wave1[NLAMBDA]; /* in nm */
	double presdummy, cross[NLAMBDA];
	double pfitting, tfitting;
	
	// Get the cache entry for this species
	OpacityFile* opac_file = get_opacity_file(Fname);

	// Check if we need to load the data
	if (!opac_file->is_loaded) {
		printf("Reading opacity file: %s\n", Fname);
		
		// Allocate opacity tensor if first time
		if (opac_file->opac == NULL) {
			opac_file->opac = f3tensor(0, NLAMBDA-1, 0, NTEMP-1, 0, NPRESSURE-1);
		}
		
		fim = fopen(Fname, "r");
		if (!fim) {
			printf("Error: Could not open opacity file %s\n", Fname);
			// Initialize with zeros if file not found
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					xsc[i][j] = 0.0;
				}
			}
			return;
		}
		
		/* Header Lines */
		for (i=0; i<NTEMP-1; i++) {
			fscanf(fim, "%lf", opac_file->temp+i);
		}
		i=NTEMP-1;
		fscanf(fim, "%lf\n", opac_file->temp+i);
		
		for (i=0; i<NPRESSURE-1; i++) {
			fscanf(fim, "%le", opac_file->pres+i);
		}
		i=NPRESSURE-1;
		fscanf(fim, "%le\n", opac_file->pres+i);
		
		/* Read in data */
		for (i=0; i<NLAMBDA; i++) {
			fscanf(fim, "%le\n", opac_file->wave+i);
			for (j=0; j<NPRESSURE; j++) {
				fscanf(fim, "%le", &presdummy);
				for (k=0; k<NTEMP-1; k++) {
					fscanf(fim, "%le", &opac_file->opac[i][k][j]);
				}
				k=NTEMP-1;
				fscanf(fim, "%le\n", &opac_file->opac[i][k][j]);
			}
		}
		
		fclose(fim);
		opac_file->is_loaded = 1;
	} else {
		printf("Using cached opacity data for: %s\n", opac_file->species_name);
	}
	
	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = opac_file->wave[i]*1.0E+9; /* convert to nm */
	}
	
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i] > opac_file->pres[NPRESSURE-1]) {
				pfitting = opac_file->pres[NPRESSURE-1];
			} else if (pl[i] < opac_file->pres[0]) {
				pfitting = opac_file->pres[0];
			} else {
				pfitting = pl[i];
			}
			
			if (tl[i] > opac_file->temp[NTEMP-1]) {
				tfitting = opac_file->temp[NTEMP-1];
			} else if (tl[i] < opac_file->temp[0]) {
				tfitting = opac_file->temp[0];
			} else {
				tfitting = tl[i];
			}
			
			cross[j] = fmax(Interpolation2D(tfitting, pfitting, 
										  opac_file->temp, NTEMP, 
										  opac_file->pres, NPRESSURE, 
										  *(opac_file->opac+j))*1.0E+4, 0.0);
			xsc[i][j] = cross[j];
		}
	}
}

// Function to reinterpolate opacity data without re-reading the file
void reinterpolate_opacities_by_file(double **xsc, char Fname[]) {
	// Find the cache entry for this species
	OpacityFile* opac_file = get_opacity_file(Fname);
	
	if (!opac_file->is_loaded || !opac_file->opac) {
		// If data not loaded, load it
		readcross(Fname, xsc);
		return;
	}

	int i, j;
	double wave1[NLAMBDA]; /* in nm */
	double pfitting, tfitting;
	double cross[NLAMBDA];

	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = opac_file->wave[i]*1.0E+9; /* convert to nm */
	}
	
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i] > opac_file->pres[NPRESSURE-1]) {
				pfitting = opac_file->pres[NPRESSURE-1];
			} else if (pl[i] < opac_file->pres[0]) {
				pfitting = opac_file->pres[0];
			} else {
				pfitting = pl[i];
			}
			
			if (tl[i] > opac_file->temp[NTEMP-1]) {
				tfitting = opac_file->temp[NTEMP-1];
			} else if (tl[i] < opac_file->temp[0]) {
				tfitting = opac_file->temp[0];
			} else {
				tfitting = tl[i];
			}
			
			cross[j] = fmax(Interpolation2D(tfitting, pfitting, 
										  opac_file->temp, NTEMP, 
										  opac_file->pres, NPRESSURE, 
										  *(opac_file->opac+j))*1.0E+4, 0.0);
			xsc[i][j] = cross[j];
		}
	}
}
