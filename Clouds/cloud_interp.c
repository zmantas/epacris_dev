#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double interp1(double* x, double* y, double xi) {
    int size = sizeof(x) / sizeof(x[0])
    if (xi <= x[0]) {
        return y[0];
    } else if (xi >= x[size - 1]) {
        return y[size - 1];
    }

    int i;
    for (i = 0; i < size - 1; i++) {
        if (xi >= x[i] && xi <= x[i + 1]) {
            // Perform linear interpolation
            double t = (xi - x[i]) / (x[i + 1] - x[i]);
            return y[i] + t * (y[i + 1] - y[i]);
        }
    }
}

int main() {
    FILE *file;
    int i, j;
    int rows = 16000; // Adjust as necessary
    int cols = 21; // Adjust as necessary

    double H2OL_r[rows];
    double H2OL_c[rows][cols];
    double H2OL_a[rows][cols];
    double H2OL_g[rows][cols];

    double H2OI_r[rows];
    double H2OI_c[rows][cols];
    double H2OI_a[rows][cols];
    double H2OI_g[rows][cols];

    // Read CrossP/cross_water_wavelength2.dat
    file = fopen("CrossP/cross_water_wavelength2.dat", "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    for (i = 0; i < rows; i++) {
        fscanf(file, "%lf", &H2OL_r[i]);
        for (j = 0; j < cols; j++) {
            fscanf(file, "%lf", &H2OL_c[i][j]);
        }
    }
    fclose(file);

    // Read CrossP/albedo_water_wavelength2.dat
    file = fopen("CrossP/albedo_water_wavelength2.dat", "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            fscanf(file, "%lf", &H2OL_a[i][j]);
        }
    }
    fclose(file);

    // Read CrossP/geo_water_wavelength2.dat
    file = fopen("CrossP/geo_water_wavelength2.dat", "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            fscanf(file, "%lf", &H2OL_g[i][j]);
        }
    }
    fclose(file);

    // Read CrossP/cross_water_wavelength2.dat for H2OI
    file = fopen("CrossP/cross_water_wavelength2.dat", "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    for (i = 0; i < rows; i++) {
        fscanf(file, "%lf", &H2OI_r[i]);
        for (j = 0; j < cols; j++) {
            fscanf(file, "%lf", &H2OI_c[i][j]);
        }
    }
    fclose(file);

    // Read CrossP/albedo_water_wavelength2.dat for H2OI
    file = fopen("CrossP/albedo_water_wavelength2.dat", "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            fscanf(file, "%lf", &H2OI_a[i][j]);
        }
    }
    fclose(file);

    // Read CrossP/geo_water_wavelength2.dat for H2OI
    file = fopen("CrossP/geo_water_wavelength2.dat", "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            fscanf(file, "%lf", &H2OI_g[i][j]);
        }
    }
    fclose(file);
    int zl;
    int SIZE;
    int j, i;
    double PI = 3.1416
    double sig = 2.0;
    double r2, r0, VP;
    double cloudmden[zl];
    double cloudden[zl];
    double particlemsize[zl];
    double particlesize[lzl];
    double tl[zl];
    double NH3I_r[SIZE], NH3I_c[SIZE], NH3I_a[SIZE], NH3I_g[SIZE];
    double H2OI_r[SIZE], H2OI_c[SIZE], H2OI_a[SIZE], H2OI_g[SIZE];
    double H2OL_r[SIZE], H2OL_c[SIZE], H2OL_a[SIZE], H2OL_g[SIZE];
    double croa[zl][SIZE], alba[zl][SIZE], geoa[zl][SIZE];
    double crow[zl][SIZE], albw[zl][SIZE], geow[zl][SIZE];

    for (j = 0; j < zl; j++) {
        r2 = particlemsize[j];
        if (cloudmden[j] < 1e-12) {
            for (i = 0; i < SIZE; i++) {
                croa[j][i] = 0;
                alba[j][i] = 1;
                geoa[j][i] = 0;
            }
        } else {
            r0 = r2 * exp(-log(sig) * log(sig));
            VP = 4 * PI / 3 * pow(r2 * 1.0E-6 * exp(0.5 * log(sig) * log(sig)), 3) * 1.0E+6 * 0.87;   %g
            for (i = 0; i < SIZE; i++) {
                double log_r0 = log10(fmax(0.01, fmin(r0, 100)));
                croa[j][i] = cloudmden[j] / VP * 1.0E-3 * pow(10, interp1(log10(NH3I_r), log10(NH3I_c), log10(fmax(0.01, fmin(r0, 100)))));
                alba[j][i] = interp1(log10(NH3I_r), NH3I_a, log10(fmax(0.01, fmin(r0, 100))));
                geoa[j][i] = interp1(log10(NH3I_r), NH3I_g, log10(fmax(0.01, fmin(r0, 100))));
            }
        }

        r2 = particlesize[j];
        if (cloudden[j] < 1e-12) {
            for (i = 0; i < SIZE; i++) {
                crow[j][i] = 0;
                albw[j][i] = 1;
                geow[j][i] = 0;
            }
        } else {
            r0 = r2 * exp(-log(sig) * log(sig));
            if (tl[j] < 273.16) { // ice
                VP = 4 * PI / 3 * pow(r2 * 1.0E-6 * exp(0.5 * log(sig) * log(sig)), 3) * 1.0E+6 * 0.92;
                for (i = 0; i < SIZE; i++) {
                    double log_r0 = log10(fmax(0.01, fmin(r0, 100)));
                    crow[j][i] = cloudden[j] / VP * 1.0E-3 * pow(10, interp1(log10(H2OI_r), log10(H2OI_c), log10(fmax(0.01, fmin(r0, 100)))));
                    albw[j][i] = interp1(log10(H2OI_r), H2OI_a, log10(fmax(0.01, fmin(r0, 100))));
                    geow[j][i] = interp1(log10(H2OI_r), H2OI_g, log10(fmax(0.01, fmin(r0, 100))));
                }
            } else { // liquid
                VP = 4 * PI / 3 * pow(r2 * 1.0E-6 * exp(0.5 * log(sig) * log(sig)), 3) * 1.0E+6 * 1.0;
                for (i = 0; i < SIZE; i++) {
                    double log_r0 = log10(fmax(0.01, fmin(r0, 100)));
                    crow[j][i] = cloudden[j] / VP * 1.0E-3 * pow(10, interp1(log10(H2OL_r), log10(H2OL_c), log10(fmax(0.01, fmin(r0, 100)))));
                    albw[j][i] = interp1(log10(H2OL_r), H2OL_a, log10(fmax(0.01, fmin(r0, 100))));
                    geow[j][i] = interp1(log10(H2OL_r), H2OL_g, log10(fmax(0.01, fmin(r0, 100))));
                }
            }
        }
    }

}


