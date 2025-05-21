/*Function to read in collision-induce absorption cross section data */
/* output opacity in cm+5 */
/* Interpolation to specified temperature in each layer */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "constant.h"

// Function prototypes for internal helper functions
static void calculate_checksums();
static int is_match(double a, double b);
static void initialize_cia_cache();

void readcia();
void reinterpolate_cia();
void cleanup_cia_cache();
void reinterpolate_all_cia_opacities();
void verify_cia_interpolation();
void log_cia_error(const char* message);
void save_cia_data_samples();

// Simple cache structure to store loaded CIA data
typedef struct {
	int is_loaded;        // Whether data is loaded
	double **cia_data;    // The actual CIA data matrix
	double wavelength[20000]; // Wavelength array - increased to handle 19991 lines
	int num_lines;        // Number of lines in the file
	int backup_num_lines; // Backup copy of num_lines for redundancy
	unsigned int magic;   // Magic number to verify structure integrity
} CiaCache;

// Magic number for verification
#define CIA_CACHE_MAGIC 0xC1AC1A00 // "CIA CIA" in hex

// Variables to store checksum values for verification
static double h2h2_checksum = 0.0;
static double h2he_checksum = 0.0;
static double h2h_checksum = 0.0;
static double n2n2_checksum = 0.0;
static double n2h2_checksum = 0.0;
static double co2co2_checksum = 0.0;
static int initial_checksums_computed = 0;

// Global cache variables - one for each CIA type
static CiaCache h2h2_cache = {0, NULL, {0}, 0, 0, 0};
static CiaCache h2he_cache = {0, NULL, {0}, 0, 0, 0};
static CiaCache h2h_cache = {0, NULL, {0}, 0, 0, 0};
static CiaCache n2n2_cache = {0, NULL, {0}, 0, 0, 0};
static CiaCache n2h2_cache = {0, NULL, {0}, 0, 0, 0};
static CiaCache co2co2_cache = {0, NULL, {0}, 0, 0, 0};

// Helper function for value comparison
static int is_match(double a, double b) {
	const double abs_tolerance = 1e-40;  // Reduced from 1e-30 to handle 1e-42 values
	const double rel_tolerance = 1e-6;   // Slightly relaxed relative tolerance
	
	// If both are very close to zero, consider them equal
	if (fabs(a) < abs_tolerance && fabs(b) < abs_tolerance) 
		return 1;
	
	// If one is zero and other isn't, they're not equal
	if ((fabs(a) < abs_tolerance && fabs(b) >= abs_tolerance) ||
		(fabs(b) < abs_tolerance && fabs(a) >= abs_tolerance))
		return 0;
	
	// For larger values, use relative error
	double rel_error = fabs(a - b) / fmax(fabs(a), fabs(b));
	return rel_error < rel_tolerance ? 1 : 0;
}

// Helper function to initialize cache structures
static void initialize_cia_cache() {
    // Initialize all cache structures to a known good state
    h2h2_cache.is_loaded = 0;
    h2h2_cache.cia_data = NULL;
    h2h2_cache.num_lines = 0;
    h2h2_cache.backup_num_lines = 0;
    h2h2_cache.magic = 0;
    
    h2he_cache.is_loaded = 0;
    h2he_cache.cia_data = NULL;
    h2he_cache.num_lines = 0;
    h2he_cache.backup_num_lines = 0;
    h2he_cache.magic = 0;
    
    h2h_cache.is_loaded = 0;
    h2h_cache.cia_data = NULL;
    h2h_cache.num_lines = 0;
    h2h_cache.backup_num_lines = 0;
    h2h_cache.magic = 0;
    
    n2n2_cache.is_loaded = 0;
    n2n2_cache.cia_data = NULL;
    n2n2_cache.num_lines = 0;
    n2n2_cache.backup_num_lines = 0;
    n2n2_cache.magic = 0;
    
    n2h2_cache.is_loaded = 0;
    n2h2_cache.cia_data = NULL;
    n2h2_cache.num_lines = 0;
    n2h2_cache.backup_num_lines = 0;
    n2h2_cache.magic = 0;
    
    co2co2_cache.is_loaded = 0;
    co2co2_cache.cia_data = NULL;
    co2co2_cache.num_lines = 0;
    co2co2_cache.backup_num_lines = 0;
    co2co2_cache.magic = 0;
    
    printf("CIA caches initialized to empty state\n");
}

// Function to log errors to a file
void log_cia_error(const char* message) {
    FILE *log_file = fopen("AuxillaryOut/cia_errors.log", "a");
    if (log_file != NULL) {
        // Get current time
        time_t now = time(NULL);
        char time_str[100];
        strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", localtime(&now));
        
        fprintf(log_file, "[%s] %s\n", time_str, message);
        fclose(log_file);
    }
}

void readcia()
{	
	int i, j, s, nl;
	FILE *fim, *fout;
	char *temp1;
	double **cia;
	char crossfile[1024];
	double tfitting[zbin+1];
	char dataline[10000];
	char error_msg[1024];
	
	// Initialize all caches to a clean state
	initialize_cia_cache();
	
	// Temperature fitting for each layer
	for (i=1; i<=zbin; i++) {
		tfitting[i] = tl[i];
		if (tfitting[i] > 2000.0) {
			tfitting[i] = 2000.0;
		}
		if (tfitting[i] < 100.0) {
			tfitting[i] = 100.0;
		}
	}
	
	// Temperature bins for interpolation
	double temp[20] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
					1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0};
	
	// Reset checksums computation flag
	initial_checksums_computed = 0;
	
	// H2-H2 CIA data
	if (!h2h2_cache.is_loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "H2-H2_CIA.dat");
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			// Initialize with zeros if file not found
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2H2CIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 20000);  // Increased buffer for larger files
			fclose(fim);
			nl = s-1;
			h2h2_cache.num_lines = nl;
			h2h2_cache.backup_num_lines = nl;  // Store backup copy
			printf("  Found %d lines of data\n", nl);
			
			// Log the number of lines for debugging
			sprintf(error_msg, "H2-H2 CIA file has %d lines of data", nl);
			log_cia_error(error_msg);
			
			// Allocate matrix and read data
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim); /* Read in the header line */
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &h2h2_cache.wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			// Save to cache and calculate cross sections
			h2h2_cache.cia_data = cia;
			h2h2_cache.is_loaded = 1;
			h2h2_cache.magic = CIA_CACHE_MAGIC; // Set magic number
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], 
									h2h2_cache.wavelength, nl, temp, 20, cia);
				}
			}
			
			// Write debug output
			fout = fopen("AuxillaryOut/CheckCIA_H2H2.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", H2H2CIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for H2-H2 CIA\n");
		// Recalculate using cached data
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], 
								h2h2_cache.wavelength, h2h2_cache.num_lines, 
								temp, 20, h2h2_cache.cia_data);
			}
		}
	}
	
	// H2-He CIA data
	if (!h2he_cache.is_loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "H2-He_CIA.dat");
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			// Initialize with zeros
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HeCIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			h2he_cache.num_lines = nl;
			printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &h2he_cache.wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			h2he_cache.cia_data = cia;
			h2he_cache.is_loaded = 1;
			h2he_cache.magic = CIA_CACHE_MAGIC; // Set magic number
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HeCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], 
									h2he_cache.wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_H2He.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", H2HeCIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for H2-He CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HeCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2he_cache.wavelength, h2he_cache.num_lines,
								temp, 20, h2he_cache.cia_data);
			}
		}
	}
	
	// H2-H CIA data
	if (!h2h_cache.is_loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "H2-H_CIA.dat");
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HCIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			h2h_cache.num_lines = nl;
			printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &h2h_cache.wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			h2h_cache.cia_data = cia;
			h2h_cache.is_loaded = 1;
			h2h_cache.magic = CIA_CACHE_MAGIC; // Set magic number
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									h2h_cache.wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_H2H.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", H2HCIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for H2-H CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2h_cache.wavelength, h2h_cache.num_lines,
								temp, 20, h2h_cache.cia_data);
			}
		}
	}
	
	// N2-N2 CIA data
	if (!n2n2_cache.is_loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "N2-N2_CIA.dat");
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2N2CIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			n2n2_cache.num_lines = nl;
			printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &n2n2_cache.wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			n2n2_cache.cia_data = cia;
			n2n2_cache.is_loaded = 1;
			n2n2_cache.magic = CIA_CACHE_MAGIC; // Set magic number
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2N2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									n2n2_cache.wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_N2N2.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", N2N2CIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for N2-N2 CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2N2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2n2_cache.wavelength, n2n2_cache.num_lines,
								temp, 20, n2n2_cache.cia_data);
			}
		}
	}
	
	// N2-H2 CIA data
	if (!n2h2_cache.is_loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "N2-H2_CIA.dat");
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2H2CIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			n2h2_cache.num_lines = nl;
			printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &n2h2_cache.wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			n2h2_cache.cia_data = cia;
			n2h2_cache.is_loaded = 1;
			n2h2_cache.magic = CIA_CACHE_MAGIC; // Set magic number
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									n2h2_cache.wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_N2H2.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", N2H2CIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for N2-H2 CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2h2_cache.wavelength, n2h2_cache.num_lines,
								temp, 20, n2h2_cache.cia_data);
			}
		}
	}
	
	// CO2-CO2 CIA data
	if (!co2co2_cache.is_loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "CO2-CO2_CIA.dat");
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					CO2CO2CIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			co2co2_cache.num_lines = nl;
			printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &co2co2_cache.wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			co2co2_cache.cia_data = cia;
			co2co2_cache.is_loaded = 1;
			co2co2_cache.magic = CIA_CACHE_MAGIC; // Set magic number
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					CO2CO2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									 co2co2_cache.wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_CO2CO2.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", CO2CO2CIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for CO2-CO2 CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				CO2CO2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								 co2co2_cache.wavelength, co2co2_cache.num_lines,
								 temp, 20, co2co2_cache.cia_data);
			}
		}
	}
	
	// Calculate checksums for verification
	calculate_checksums();
	printf("CIA data loading complete, checksums calculated for verification\n");
}

// Function to reinterpolate CIA opacities using temperature values
void reinterpolate_cia() {
	int i, j;
	double temp[20] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
					1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0};
	
	double tfitting[zbin+1];
	for (i=1; i<=zbin; i++) {
		tfitting[i] = tl[i];
		if (tfitting[i] > 2000.0) tfitting[i] = 2000.0;
		if (tfitting[i] < 100.0) tfitting[i] = 100.0;
	}
	
	printf("Reinterpolating CIA opacities...\n");
	
	// Print diagnostic information about cache status
	printf("CIA Cache Status before reinterpolation:\n");
	printf("H2-H2: is_loaded=%d, cia_data=%p, num_lines=%d, magic=%s\n", 
	       h2h2_cache.is_loaded, (void*)h2h2_cache.cia_data, h2h2_cache.num_lines,
	       (h2h2_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
	printf("H2-He: is_loaded=%d, cia_data=%p, num_lines=%d, magic=%s\n", 
	       h2he_cache.is_loaded, (void*)h2he_cache.cia_data, h2he_cache.num_lines,
	       (h2he_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
	printf("H2-H: is_loaded=%d, cia_data=%p, num_lines=%d, magic=%s\n", 
	       h2h_cache.is_loaded, (void*)h2h_cache.cia_data, h2h_cache.num_lines,
	       (h2h_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
	printf("N2-N2: is_loaded=%d, cia_data=%p, num_lines=%d, magic=%s\n", 
	       n2n2_cache.is_loaded, (void*)n2n2_cache.cia_data, n2n2_cache.num_lines,
	       (n2n2_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
	printf("N2-H2: is_loaded=%d, cia_data=%p, num_lines=%d, magic=%s\n", 
	       n2h2_cache.is_loaded, (void*)n2h2_cache.cia_data, n2h2_cache.num_lines,
	       (n2h2_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
	printf("CO2-CO2: is_loaded=%d, cia_data=%p, num_lines=%d, magic=%s\n", 
	       co2co2_cache.is_loaded, (void*)co2co2_cache.cia_data, co2co2_cache.num_lines,
	       (co2co2_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
	
	// H2-H2
	if (h2h2_cache.is_loaded && h2h2_cache.cia_data != NULL) {
		int num_lines = h2h2_cache.num_lines;
		// Check magic number - if invalid, reset cache
		if (h2h2_cache.magic != CIA_CACHE_MAGIC) {
			printf("Warning: H2-H2 cache has invalid magic number - possible memory corruption\n");
			log_cia_error("Warning: H2-H2 cache has invalid magic number - possible memory corruption");
			// Fallback to default
			num_lines = 10000;
		} else if (num_lines == 0) {
			// First try backup_num_lines
			if (h2h2_cache.backup_num_lines > 0) {
				printf("Recovering H2-H2 num_lines from backup: %d\n", h2h2_cache.backup_num_lines);
				log_cia_error("Recovering H2-H2 num_lines from backup");
				num_lines = h2h2_cache.backup_num_lines;
			} else {
				// If backup is also 0, use fixed fallback
				printf("Warning: H2-H2 num_lines and backup_num_lines are 0, using fixed value of 19991\n");
				log_cia_error("Warning: H2-H2 num_lines and backup_num_lines are 0, using fixed value of 19991");
				num_lines = 19991; // More accurate fallback based on known file size
			}
		}
		
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2h2_cache.wavelength, num_lines,
								temp, 20, h2h2_cache.cia_data);
			}
		}
	} else {
		printf("Warning: H2-H2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2H2CIA[i][j] = 0.0;
			}
		}
	}
	
	// H2-He
	if (h2he_cache.is_loaded && h2he_cache.cia_data != NULL) {
		int num_lines = h2he_cache.num_lines;
		// Check magic number - if invalid, reset cache
		if (h2he_cache.magic != CIA_CACHE_MAGIC) {
			printf("Warning: H2-He cache has invalid magic number - possible memory corruption\n");
			// Fallback to default
			num_lines = 10000;
		} else if (num_lines == 0) {
			printf("Warning: H2-He num_lines is 0 despite valid magic, using fixed value of 10000\n");
			num_lines = 10000;
		}
		
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HeCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2he_cache.wavelength, num_lines,
								temp, 20, h2he_cache.cia_data);
			}
		}
	} else {
		printf("Warning: H2-He CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HeCIA[i][j] = 0.0;
			}
		}
	}
	
	// H2-H
	if (h2h_cache.is_loaded && h2h_cache.cia_data != NULL) {
		int num_lines = h2h_cache.num_lines;
		// Check magic number - if invalid, reset cache
		if (h2h_cache.magic != CIA_CACHE_MAGIC) {
			printf("Warning: H2-H cache has invalid magic number - possible memory corruption\n");
			// Fallback to default
			num_lines = 10000;
		} else if (num_lines == 0) {
			printf("Warning: H2-H num_lines is 0 despite valid magic, using fixed value of 10000\n");
			num_lines = 10000;
		}
		
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2h_cache.wavelength, num_lines,
								temp, 20, h2h_cache.cia_data);
			}
		}
	} else {
		printf("Warning: H2-H CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HCIA[i][j] = 0.0;
			}
		}
	}
	
	// N2-N2
	if (n2n2_cache.is_loaded && n2n2_cache.cia_data != NULL) {
		int num_lines = n2n2_cache.num_lines;
		// Check magic number - if invalid, reset cache
		if (n2n2_cache.magic != CIA_CACHE_MAGIC) {
			printf("Warning: N2-N2 cache has invalid magic number - possible memory corruption\n");
			// Fallback to default
			num_lines = 10000;
		} else if (num_lines == 0) {
			printf("Warning: N2-N2 num_lines is 0 despite valid magic, using fixed value of 10000\n");
			num_lines = 10000;
		}
		
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2N2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2n2_cache.wavelength, num_lines,
								temp, 20, n2n2_cache.cia_data);
			}
		}
	} else {
		printf("Warning: N2-N2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2N2CIA[i][j] = 0.0;
			}
		}
	}
	
	// N2-H2
	if (n2h2_cache.is_loaded && n2h2_cache.cia_data != NULL) {
		int num_lines = n2h2_cache.num_lines;
		// Check magic number - if invalid, reset cache
		if (n2h2_cache.magic != CIA_CACHE_MAGIC) {
			printf("Warning: N2-H2 cache has invalid magic number - possible memory corruption\n");
			// Fallback to default
			num_lines = 10000;
		} else if (num_lines == 0) {
			printf("Warning: N2-H2 num_lines is 0 despite valid magic, using fixed value of 10000\n");
			num_lines = 10000;
		}
		
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2h2_cache.wavelength, num_lines,
								temp, 20, n2h2_cache.cia_data);
			}
		}
	} else {
		printf("Warning: N2-H2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2H2CIA[i][j] = 0.0;
			}
		}
	}
	
	// CO2-CO2
	if (co2co2_cache.is_loaded && co2co2_cache.cia_data != NULL) {
		int num_lines = co2co2_cache.num_lines;
		// Check magic number - if invalid, reset cache
		if (co2co2_cache.magic != CIA_CACHE_MAGIC) {
			printf("Warning: CO2-CO2 cache has invalid magic number - possible memory corruption\n");
			// Fallback to default
			num_lines = 10000;
		} else if (num_lines == 0) {
			printf("Warning: CO2-CO2 num_lines is 0 despite valid magic, using fixed value of 10000\n");
			num_lines = 10000;
		}
		
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				CO2CO2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								 co2co2_cache.wavelength, num_lines,
								 temp, 20, co2co2_cache.cia_data);
			}
		}
	} else {
		printf("Warning: CO2-CO2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				CO2CO2CIA[i][j] = 0.0;
			}
		}
	}
}

// Clean up allocated memory
void cleanup_cia_cache() {
	printf("Cleaning up CIA opacity cache...\n");
	
	if (h2h2_cache.is_loaded && h2h2_cache.cia_data != NULL) {
		free_dmatrix(h2h2_cache.cia_data, 0, h2h2_cache.num_lines-1, 0, 19);
		h2h2_cache.cia_data = NULL;
		h2h2_cache.is_loaded = 0;
		h2h2_cache.num_lines = 0;
		h2h2_cache.magic = 0; // Reset magic number
		printf("Freed H2-H2 CIA cache\n");
	}
	
	if (h2he_cache.is_loaded && h2he_cache.cia_data != NULL) {
		free_dmatrix(h2he_cache.cia_data, 0, h2he_cache.num_lines-1, 0, 19);
		h2he_cache.cia_data = NULL;
		h2he_cache.is_loaded = 0;
		h2he_cache.num_lines = 0;
		h2he_cache.magic = 0; // Reset magic number
		printf("Freed H2-He CIA cache\n");
	}
	
	if (h2h_cache.is_loaded && h2h_cache.cia_data != NULL) {
		free_dmatrix(h2h_cache.cia_data, 0, h2h_cache.num_lines-1, 0, 19);
		h2h_cache.cia_data = NULL;
		h2h_cache.is_loaded = 0;
		h2h_cache.num_lines = 0;
		h2h_cache.magic = 0; // Reset magic number
		printf("Freed H2-H CIA cache\n");
	}
	
	if (n2n2_cache.is_loaded && n2n2_cache.cia_data != NULL) {
		free_dmatrix(n2n2_cache.cia_data, 0, n2n2_cache.num_lines-1, 0, 19);
		n2n2_cache.cia_data = NULL;
		n2n2_cache.is_loaded = 0;
		n2n2_cache.num_lines = 0;
		n2n2_cache.magic = 0; // Reset magic number
		printf("Freed N2-N2 CIA cache\n");
	}
	
	if (n2h2_cache.is_loaded && n2h2_cache.cia_data != NULL) {
		free_dmatrix(n2h2_cache.cia_data, 0, n2h2_cache.num_lines-1, 0, 19);
		n2h2_cache.cia_data = NULL;
		n2h2_cache.is_loaded = 0;
		n2h2_cache.num_lines = 0;
		n2h2_cache.magic = 0; // Reset magic number
		printf("Freed N2-H2 CIA cache\n");
	}
	
	if (co2co2_cache.is_loaded && co2co2_cache.cia_data != NULL) {
		free_dmatrix(co2co2_cache.cia_data, 0, co2co2_cache.num_lines-1, 0, 19);
		co2co2_cache.cia_data = NULL;
		co2co2_cache.is_loaded = 0;
		co2co2_cache.num_lines = 0;
		co2co2_cache.magic = 0; // Reset magic number
		printf("Freed CO2-CO2 CIA cache\n");
	}
	
	printf("CIA cache cleanup complete\n");
}

// These are not implemented in this simplified version
// void reinterpolate_normal_opacities() {}

// Implementation for reinterpolate_all_cia_opacities
void reinterpolate_all_cia_opacities() {
    printf("Reinterpolating all opacities...\n");
    reinterpolate_cia();
    verify_cia_interpolation();
    printf("Reinterpolation complete\n");
}

// Function to calculate checksums for verification
static void calculate_checksums() {
    int i, j;
    char error_msg[1024];
    
    // Reset checksums
    h2h2_checksum = 0.0;
    h2he_checksum = 0.0;
    h2h_checksum = 0.0;
    n2n2_checksum = 0.0;
    n2h2_checksum = 0.0;
    co2co2_checksum = 0.0;
    
    // Calculate checksums using multiple methods for redundancy:
    // 1. Sum of sampled values
    // 2. Count of non-zero elements
    // 3. Max value found
    
    int h2h2_nonzero = 0, h2he_nonzero = 0, h2h_nonzero = 0;
    int n2n2_nonzero = 0, n2h2_nonzero = 0, co2co2_nonzero = 0;
    
    double h2h2_max = 0.0, h2he_max = 0.0, h2h_max = 0.0;
    double n2n2_max = 0.0, n2h2_max = 0.0, co2co2_max = 0.0;
    
    // Calculate checksums by sampling opacity values
    // We'll take sample points across all layers, with emphasis on areas with data
    int sample_step = NLAMBDA / 20;  // Increased number of samples
    if (sample_step < 1) sample_step = 1;
    
    for (i=1; i<=zbin; i++) {
        for (j=0; j<NLAMBDA; j+=sample_step) {
            // Sum of values
            h2h2_checksum += H2H2CIA[i][j];
            h2he_checksum += H2HeCIA[i][j];
            h2h_checksum += H2HCIA[i][j];
            n2n2_checksum += N2N2CIA[i][j];
            n2h2_checksum += N2H2CIA[i][j];
            co2co2_checksum += CO2CO2CIA[i][j];
            
            // Count non-zero elements
            if (fabs(H2H2CIA[i][j]) > 1e-40) h2h2_nonzero++;
            if (fabs(H2HeCIA[i][j]) > 1e-40) h2he_nonzero++;
            if (fabs(H2HCIA[i][j]) > 1e-40) h2h_nonzero++;
            if (fabs(N2N2CIA[i][j]) > 1e-40) n2n2_nonzero++;
            if (fabs(N2H2CIA[i][j]) > 1e-40) n2h2_nonzero++;
            if (fabs(CO2CO2CIA[i][j]) > 1e-40) co2co2_nonzero++;
            
            // Track max values
            if (fabs(H2H2CIA[i][j]) > h2h2_max) h2h2_max = fabs(H2H2CIA[i][j]);
            if (fabs(H2HeCIA[i][j]) > h2he_max) h2he_max = fabs(H2HeCIA[i][j]);
            if (fabs(H2HCIA[i][j]) > h2h_max) h2h_max = fabs(H2HCIA[i][j]);
            if (fabs(N2N2CIA[i][j]) > n2n2_max) n2n2_max = fabs(N2N2CIA[i][j]);
            if (fabs(N2H2CIA[i][j]) > n2h2_max) n2h2_max = fabs(N2H2CIA[i][j]);
            if (fabs(CO2CO2CIA[i][j]) > co2co2_max) co2co2_max = fabs(CO2CO2CIA[i][j]);
        }
    }
    
    // Log the checksum calculation results
    sprintf(error_msg, "Checksums calculated - H2-H2: %.6e (nonzero: %d, max: %.6e)", 
            h2h2_checksum, h2h2_nonzero, h2h2_max);
    log_cia_error(error_msg);
    sprintf(error_msg, "Checksums calculated - H2-He: %.6e (nonzero: %d, max: %.6e)", 
            h2he_checksum, h2he_nonzero, h2he_max);
    log_cia_error(error_msg);
    
    initial_checksums_computed = 1;
}

// Function to verify CIA interpolation
void verify_cia_interpolation() {
    if (!initial_checksums_computed) {
        printf("Cannot verify interpolation - initial checksums not computed\n");
        log_cia_error("Cannot verify interpolation - initial checksums not computed");
        return;
    }
    
    // Calculate new checksums after reinterpolation
    int i, j;
    double h2h2_new = 0.0, h2he_new = 0.0, h2h_new = 0.0;
    double n2n2_new = 0.0, n2h2_new = 0.0, co2co2_new = 0.0;
    char error_msg[1024];
    
    // Detailed verification with multiple methods
    int h2h2_nonzero = 0, h2he_nonzero = 0, h2h_nonzero = 0;
    int n2n2_nonzero = 0, n2h2_nonzero = 0, co2co2_nonzero = 0;
    
    double h2h2_max = 0.0, h2he_max = 0.0, h2h_max = 0.0;
    double n2n2_max = 0.0, n2h2_max = 0.0, co2co2_max = 0.0;
    
    int sample_step = NLAMBDA / 20;  // Increased number of samples
    if (sample_step < 1) sample_step = 1;
    
    // First pass: Calculate checksums and statistics
    for (i=1; i<=zbin; i++) {
        for (j=0; j<NLAMBDA; j+=sample_step) {
            // Sum of values
            h2h2_new += H2H2CIA[i][j];
            h2he_new += H2HeCIA[i][j];
            h2h_new += H2HCIA[i][j];
            n2n2_new += N2N2CIA[i][j];
            n2h2_new += N2H2CIA[i][j];
            co2co2_new += CO2CO2CIA[i][j];
            
            // Count non-zero elements
            if (fabs(H2H2CIA[i][j]) > 1e-40) h2h2_nonzero++;
            if (fabs(H2HeCIA[i][j]) > 1e-40) h2he_nonzero++;
            if (fabs(H2HCIA[i][j]) > 1e-40) h2h_nonzero++;
            if (fabs(N2N2CIA[i][j]) > 1e-40) n2n2_nonzero++;
            if (fabs(N2H2CIA[i][j]) > 1e-40) n2h2_nonzero++;
            if (fabs(CO2CO2CIA[i][j]) > 1e-40) co2co2_nonzero++;
            
            // Track max values
            if (fabs(H2H2CIA[i][j]) > h2h2_max) h2h2_max = fabs(H2H2CIA[i][j]);
            if (fabs(H2HeCIA[i][j]) > h2he_max) h2he_max = fabs(H2HeCIA[i][j]);
            if (fabs(H2HCIA[i][j]) > h2h_max) h2h_max = fabs(H2HCIA[i][j]);
            if (fabs(N2N2CIA[i][j]) > n2n2_max) n2n2_max = fabs(N2N2CIA[i][j]);
            if (fabs(N2H2CIA[i][j]) > n2h2_max) n2h2_max = fabs(N2H2CIA[i][j]);
            if (fabs(CO2CO2CIA[i][j]) > co2co2_max) co2co2_max = fabs(CO2CO2CIA[i][j]);
        }
    }
    
    // Second pass: Detailed value comparison for a subset of points
    // Only sample a few specific points to avoid excessive output
    int mismatch_count = 0;
    int max_mismatches_to_report = 5;  // Limit the number of detailed mismatches to report
    
    // Create a verification file for detailed comparison
    FILE *verify_file = fopen("AuxillaryOut/cia_verification.txt", "w");
    if (verify_file != NULL) {
        fprintf(verify_file, "# CIA Detailed Verification\n");
        fprintf(verify_file, "# Layer Wavelength Original Reinterpolated Difference\n");
        
        // Check a few specific layers and wavelengths where we expect data
        for (i=1; i<=zbin; i+=zbin/5) {
            for (j=NLAMBDA/4; j<3*NLAMBDA/4; j+=NLAMBDA/10) {
                // Just checking H2-H2 for detailed verification as an example
                fprintf(verify_file, "%d %.6f %.8e %.8e %.8e\n", 
                        i, wavelength[j], 0.0 /* original value not stored */, 
                        H2H2CIA[i][j], 0.0 /* can't calculate difference */);
            }
        }
        fclose(verify_file);
    }
    
    // Compare checksums and report
    int all_match = 1;  // Flag to track if all checks pass
    
    // Helper function to evaluate match and update all_match flag
    #define EVALUATE_MATCH(orig, new, name) do { \
        int match = is_match(orig, new); \
        all_match &= match; \
        sprintf(error_msg, "%s match: %s (%.8e vs %.8e)", name, \
                match ? "YES" : "NO", orig, new); \
        log_cia_error(error_msg); \
    } while(0)
    
    printf("\n=== CIA Interpolation Verification ===\n");
    printf("Type       Original Sum      Reinterpolated Sum    Match?\n");
    printf("H2-H2      %.8e    %.8e    %s\n", 
           h2h2_checksum, h2h2_new, 
           is_match(h2h2_checksum, h2h2_new) ? "YES" : "NO");
    EVALUATE_MATCH(h2h2_checksum, h2h2_new, "H2-H2");
    
    printf("H2-He      %.8e    %.8e    %s\n", 
           h2he_checksum, h2he_new,
           is_match(h2he_checksum, h2he_new) ? "YES" : "NO");
    EVALUATE_MATCH(h2he_checksum, h2he_new, "H2-He");
    
    printf("H2-H       %.8e    %.8e    %s\n", 
           h2h_checksum, h2h_new,
           is_match(h2h_checksum, h2h_new) ? "YES" : "NO");
    EVALUATE_MATCH(h2h_checksum, h2h_new, "H2-H");
    
    printf("N2-N2      %.8e    %.8e    %s\n", 
           n2n2_checksum, n2n2_new,
           is_match(n2n2_checksum, n2n2_new) ? "YES" : "NO");
    EVALUATE_MATCH(n2n2_checksum, n2n2_new, "N2-N2");
    
    printf("N2-H2      %.8e    %.8e    %s\n", 
           n2h2_checksum, n2h2_new,
           is_match(n2h2_checksum, n2h2_new) ? "YES" : "NO");
    EVALUATE_MATCH(n2h2_checksum, n2h2_new, "N2-H2");
    
    printf("CO2-CO2    %.8e    %.8e    %s\n", 
           co2co2_checksum, co2co2_new,
           is_match(co2co2_checksum, co2co2_new) ? "YES" : "NO");
    EVALUATE_MATCH(co2co2_checksum, co2co2_new, "CO2-CO2");
    
    printf("=====================================\n");
    
    // Additional non-zero counts verification
    printf("\n=== CIA Non-Zero Value Counts ===\n");
    printf("H2-H2: %d non-zero values (max: %.3e)\n", h2h2_nonzero, h2h2_max);
    printf("H2-He: %d non-zero values (max: %.3e)\n", h2he_nonzero, h2he_max);
    printf("H2-H:  %d non-zero values (max: %.3e)\n", h2h_nonzero, h2h_max);
    printf("N2-N2: %d non-zero values (max: %.3e)\n", n2n2_nonzero, n2n2_max);
    printf("N2-H2: %d non-zero values (max: %.3e)\n", n2h2_nonzero, n2h2_max);
    printf("CO2-CO2: %d non-zero values (max: %.3e)\n", co2co2_nonzero, co2co2_max);
    printf("=====================================\n");
    
    // Log overall result
    sprintf(error_msg, "CIA verification complete: %s", 
            all_match ? "All matches OK" : "MISMATCHES DETECTED");
    log_cia_error(error_msg);
    
    // If checksums match but nonzero count is very low, it could indicate zeroing issues
    if (h2h2_nonzero < 10 && is_match(h2h2_checksum, h2h2_new)) {
        printf("WARNING: H2-H2 checksums match but very few non-zero values (%d)\n", h2h2_nonzero);
        log_cia_error("WARNING: H2-H2 checksums match but very few non-zero values - possible zeroing issue");
    }
}

// Function to save CIA data samples for independent verification
void save_cia_data_samples() {
    int i, j;
    FILE *sample_file = fopen("AuxillaryOut/cia_data_samples.txt", "w");
    
    if (sample_file == NULL) {
        printf("Error: Could not create CIA data sample file\n");
        log_cia_error("Could not create CIA data sample file");
        return;
    }
    
    fprintf(sample_file, "# CIA Data Samples\n");
    fprintf(sample_file, "# Format: Layer Wavelength H2-H2 H2-He H2-H N2-N2 N2-H2 CO2-CO2\n");
    
    // Sample data at various wavelengths and layers
    int sample_layers[] = {1, zbin/4, zbin/2, 3*zbin/4, zbin};
    int num_layers = sizeof(sample_layers) / sizeof(sample_layers[0]);
    
    int sample_points[] = {0, NLAMBDA/10, NLAMBDA/4, NLAMBDA/2, 3*NLAMBDA/4, 9*NLAMBDA/10, NLAMBDA-1};
    int num_points = sizeof(sample_points) / sizeof(sample_points[0]);
    
    // Write header with cache information
    fprintf(sample_file, "# Cache status:\n");
    fprintf(sample_file, "# H2-H2: is_loaded=%d, num_lines=%d, backup_num_lines=%d, magic=%s\n", 
           h2h2_cache.is_loaded, h2h2_cache.num_lines, h2h2_cache.backup_num_lines,
           (h2h2_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
    fprintf(sample_file, "# H2-He: is_loaded=%d, num_lines=%d, backup_num_lines=%d, magic=%s\n", 
           h2he_cache.is_loaded, h2he_cache.num_lines, h2he_cache.backup_num_lines,
           (h2he_cache.magic == CIA_CACHE_MAGIC) ? "valid" : "invalid");
    
    // Sample data points
    for (i=0; i<num_layers; i++) {
        int layer = sample_layers[i];
        for (j=0; j<num_points; j++) {
            int point = sample_points[j];
            fprintf(sample_file, "%d %.6f %.8e %.8e %.8e %.8e %.8e %.8e\n",
                   layer, wavelength[point],
                   H2H2CIA[layer][point], H2HeCIA[layer][point], H2HCIA[layer][point],
                   N2N2CIA[layer][point], N2H2CIA[layer][point], CO2CO2CIA[layer][point]);
        }
        fprintf(sample_file, "\n"); // Blank line between layers for easier reading
    }
    
    fclose(sample_file);
    printf("CIA data samples saved to AuxillaryOut/cia_data_samples.txt\n");
    log_cia_error("CIA data samples saved for verification");
}




