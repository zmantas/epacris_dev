#include <stdio.h>
#include <math.h>

// Define constants
#define HPLANCK 6.626068E-34
#define CLIGHT 299792458
#define KB 1.3806503E-23
#define AMU 1.66053886E-27
#define AVOGADRO 6.022141E+23
#define RGAS 8.314472
#define GRAVITY 6.674E-11
#define PI 3.1416         

// Main function to calculate particle sizes
void particlesizef(double g, double T, double P, double M, double MM, double KE, double deltaP, double *r0, double *r1, double *r2, double *VP) {

    double H = KB * T / M / AMU / g;
    double u = KE / H;
    double mu = 8.76E-6 * (293.85 + 72) / (293.85 + 72) * pow(T / 293.85, 1.5);
    double lambda = 2 * mu / P / pow(8 * M * 1.0E-3 / PI / RGAS / T, 0.5);
    double deltan = deltaP / KB / T;

    // Droplet parameters
    double rho = 1.0e3; // kg m^-3
    double acc = 1.0;

    // Mass diffusion coefficient
    double D = 0.12e-4;

    // Particle Size and Number
    double Cc0= 1.0;
    double fa = 1.0;
    double sig = 2.0;

    // Iteration parameters
    double Vs;
    for (int dump = 1; dump <= 1e3; ++dump) {
        double cc = -(pow(48.0 * PI * PI, 1.0 / 3.0) * D * MM * AMU * fa * deltan / rho * exp(-pow(log(sig), 2.0)));
        double aa = rho * g / mu / pow(162.0 * PI * PI, 1.0 / 3.0) / H * Cc0* exp(-pow(log(sig), 2.0));
        double bb = -u / H;

        double V = pow((-bb + sqrt(bb * bb - 4.0 * aa * cc)) / 2.0 / aa, 3.0 / 2.0);
        double d1 = pow(6.0 * V / PI, 1.0 / 3.0) * exp(-pow(log(sig), 2.0));

        double kn = lambda / d1;
        double Cc1 = 1.0 + kn * (1.257 + 0.4 * exp(-1.1 / kn));
        double fa1 = (1.0 + kn) / (1.0 + 2.0 * kn * (1.0 + kn) / acc);
        
        if (fabs(Cc1 - Cc0) + fabs(fa1 - fa) < 0.001) {
            Vs = V;
            break;
        } else {
            Cc0= Cc1;
            fa = fa1;
        }
    }
    *r0 = pow((3.0 * Vs) / (4.0 * PI), 1.0 / 3.0) * exp(-1.5 * pow(log(sig), 2.0)) * 1.0E+6;
    *r1 = pow((3.0 * Vs) / (4.0 * PI), 1.0 / 3.0) * exp(-pow(log(sig), 2.0)) * 1.0E+6;
    *r2 = pow((3.0 * Vs) / (4.0 * PI), 1.0 / 3.0) * exp(-0.5 * pow(log(sig), 2.0)) * 1.0E+6;
    *VP = Vs * 1.0E+6;
    
}

int main() {

#define g 9.8
#define T 200.0
#define P 1000.0
#define M 2.0
#define MM 18.0
#define KE 1.0
#define deltaP 40.

    double r0, r1, r2, VP;

    // Example call to particlesizef function
    // Values for g, T, P, M, MM, KE, deltaP need to be defined based on the specific use case
    particlesizef(g, T, P, M, MM, KE, deltaP, &r0, &r1, &r2, &VP);

    // Output the results
    printf("r0: %f, r1: %f, r2: %f, VP: %e \n", r0, r1, r2, VP);

    return 0;
}

