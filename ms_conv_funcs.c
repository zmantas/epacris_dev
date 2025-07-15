#include <math.h>
#include <stdbool.h>  // For bool type
#include "constant.h"
#include "condensed_heat.h"  // Include new condensed heat capacity functions

// Note: Global variables for dynamic condensibles management are defined in the main file
// NCONDENSIBLES and CONDENSIBLES[] are declared in global_temp.h
// ALPHA_RAINOUT is a single constant value defined in AlphaAb.h

// NEW: Iterative equilibrium approach for condensation
void ms_iterative_condensation_equilibrium(int layer, double* convergence_metric);

//=========================================================
//=== Function identifiers ================================

void ms_adiabat(int lay, double lapse[], double xxHe, double* cp, double saturation_ratios[]); // calculate lapse rate (NO rainout)
void ms_rainout(int lay, double* mass_loss_ratio); // apply rainout physics
void apply_enhanced_cloud_physics(int layer, double gravity); // enhanced cloud physics (Ackerman & Marley)
void apply_realistic_sedimentation(double gravity); // apply realistic sedimentation with mass conservation
void apply_equilibrium_cloud_distribution(double gravity); // apply Ackerman & Marley (2001) exponential cloud redistribution using real particle sizes
void get_particle_properties(int species_id, double temperature, double *density, double *accommodation_coeff, double *molecular_mass); // get species-specific properties
void ms_conv_check(double tempb[],double P[],double lapse[],int isconv[],int* ncl,int* nrl);   // Checks each pair of layers for convective stability
void ms_temp_adj(double tempb[],double P[],double lapse[],int isconv[], double cp[], int ncreg[], double pot_temp[]); // Adjust temperatures due to convection

// New functions for dynamic condensation detection
void initialize_condensibles_mode(); // Initialize condensibles based on CONDENSATION_MODE
int check_species_condensible(int species_id, double temp, double pressure, double partial_pressure); // Check if species should be condensible
void update_condensibles_list(int layer); // Update condensibles list for a specific layer
void detect_condensibles_atmosphere(); // Detect condensibles across entire atmosphere
void report_condensibles_changes(int iteration); // Report changes in condensibles list

//=== Helper function identifiers =========================
double ms_psat_h2o(double temp); //calc saturation pressure of water
double ms_psat_nh3(double temp); //calc saturation pressure
double ms_psat_co(double temp); //calculate saturation pressure of carbon monoxide
double ms_psat_ch4(double temp); //calculate saturation pressure of methane
double ms_psat_co2(double temp); //calculate saturation pressure of carbon dioxide
double ms_psat_h2(double temp); //calculate saturation pressure of hydrogen
double ms_psat_o2(double temp); //calculate saturation pressure of oxygen
double ms_psat_n2(double temp); //calculate saturation pressure of nitrogen
double ms_psat_h2s(double temp); //calculate saturation pressure of hydrogen sulfide
double ms_latent(int mol, double temp); //calc saturation pressure of water

//=========================================================
//=== Functions ===========================================
void ms_adiabat(int lay, double lapse[], double xxHe, double* cp, double saturation_ratios[])
{
// Calculate adiabat acc. to Graham+2021 Eq.(1) 
// The formulation below should be valid for dry, moist, and multi-species condensate conditions in dilute and non-dilute cases.
//
// Condensibles by species number from chemistry:
// H2O : 7
// NH3 : 9
// CO  : 20
// CH4 : 21
// CO2 : 52
// H2  : 53
// O2  : 54
// N2  : 55
// Cp avail: Air(0.8 N2, 0.2 O2), CO2, He, N2, NH3, CH4, H2, O2, CO, H2O

    //LOCAL parameters
    int i,j;
    //input file    int NCONDENSIBLES = 1; //number of potential condensibles; here H2O

    //variables for adiabatic calculation
    double Xd,cp_d; // VMR of non-condensible gas, heat capacity of non-condensible gas
    double Xv[NCONDENSIBLES],Xvold,cp_v[NCONDENSIBLES]; //condensibles mol fractions in vapor form and heat cap. [Xvold is the old value of Xv]
    double Xc[NCONDENSIBLES],cp_c[NCONDENSIBLES]; //condensibles fractions in condensed form (Xc) and heat cap. (cp_c)
    double alpha[NCONDENSIBLES], beta[NCONDENSIBLES], latent[NCONDENSIBLES]; // Latent heat and B=L/RT for each condensible
    
    // MINIMAL FIX: Initialize only beta array to prevent false latent heat effects
    // (Other arrays will be properly set in the loop below)
    for (int init_i = 0; init_i < NCONDENSIBLES; init_i++) {
        beta[init_i] = 0.0;  // Critical: ensure no false latent heat
    }

    //variables for condensation
    double psat[NCONDENSIBLES]; //saturation pressures
    int saturated[NCONDENSIBLES]; //check which molecules are available for condensation
    double dp[NCONDENSIBLES]; //condensed fraction

    double cpxx; //cp * xx for all dry species for the cp_d calculation
    double MM_cp; // MM for all species we have a cp for

    double reservoirs[NCONDENSIBLES]; //keeping track of what was rained out + was available at surface at beginning

    double lapse_num, lapse_denom; //numerator and denominator of lapse rate
    double cp_num, cp_denom; //numerator and denominator of heat capacity
    double sum_beta_xv, big_sum_denom_num_left_term, big_sum_denom_num_right_term; //sum of beta*Xv and big sum of denominator terms of lapse rate

    //START function
    //print number of condensibles species
    // if(lay==1) printf("%s\n","--- Adiabat calculation ---");
    // if(lay==1) printf("%s\t%d\n","Number of condensibles species:",NCONDENSIBLES);

    Xd = 1.0; //mole fraction of non-condensible gas
    cp_d = 0.0; //heat capacity of non-condensible gas

    //Assign mol fractions:
    for (i=0; i<NCONDENSIBLES; i++)
    {
        // if(lay==1) printf("%s%d%s\t%d\n","CONDENSIBLE[",i,"] = ",CONDENSIBLES[i]);

        //if(lay==1) printf("%s%d%s\t%e\n","Xv[",i,"] = ",Xv[i]); (debugging)

        //calculate saturation pressure for each condensible species
        if(CONDENSIBLES[i]==7) psat[i] = ms_psat_h2o(tl[lay]);  // Condenses ~273K (liquid), ~180K (ice)
        if(CONDENSIBLES[i]==9) psat[i] = ms_psat_nh3(tl[lay]);  // Condenses ~195K (triple point)
        // These are added and not checked for sources
        if(CONDENSIBLES[i]==20) psat[i] = ms_psat_co(tl[lay]);  // Condenses ~68K (triple point)
        if(CONDENSIBLES[i]==21) psat[i] = ms_psat_ch4(tl[lay]); // Condenses ~90K (triple point)
        if(CONDENSIBLES[i]==52) psat[i] = ms_psat_co2(tl[lay]); // Condenses ~195-216K (dry ice)
        if(CONDENSIBLES[i]==53) psat[i] = ms_psat_h2(tl[lay]);  // Condenses ~14K (triple point)
        if(CONDENSIBLES[i]==54) psat[i] = ms_psat_o2(tl[lay]);  // Condenses ~54K (triple point)
        if(CONDENSIBLES[i]==55) psat[i] = ms_psat_n2(tl[lay]);  // Condenses ~63K (triple point)
        if(CONDENSIBLES[i]==45) psat[i] = ms_psat_h2s(tl[lay]);  // Condenses ~188K (triple point)
        /* Additional species that could be included:
        if(CONDENSIBLES[i]==43) psat[i] = ms_psat_so2(tl[lay]);  // Condenses ~198K (triple point)
        if(CONDENSIBLES[i]==99) psat[i] = ms_psat_s8(tl[lay]);   // Condenses ~388K (solid form)
        if(CONDENSIBLES[i]==27) psat[i] = ms_psat_c2h2(tl[lay]); // Condenses ~192K (triple point)
        if(CONDENSIBLES[i]==29) psat[i] = ms_psat_c2h4(tl[lay]); // Condenses ~104K (triple point)
        if(CONDENSIBLES[i]==31) psat[i] = ms_psat_c2h6(tl[lay]); // Condenses ~90K (triple point)
        if(CONDENSIBLES[i]==37) psat[i] = ms_psat_hcn(tl[lay]);  // Condenses ~260K (triple point)
        if(CONDENSIBLES[i]==100) psat[i] = ms_psat_kcl(tl[lay]); // Condenses ~1000-1100K (hot Jupiters)
        if(CONDENSIBLES[i]==101) psat[i] = ms_psat_nacl(tl[lay]); // Condenses ~1000-1100K (hot Jupiters)
        if(CONDENSIBLES[i]==103) psat[i] = ms_psat_na2s(tl[lay]); // Condenses ~700-800K (hot Jupiters)
        if(CONDENSIBLES[i]==104) psat[i] = ms_psat_mgsiO3(tl[lay]); // Condenses ~1500-1600K (very hot)
        if(CONDENSIBLES[i]==105) psat[i] = ms_psat_mg2siO4(tl[lay]); // Condenses ~1600-1700K (very hot)
        if(CONDENSIBLES[i]==106) psat[i] = ms_psat_fe(tl[lay]);  // Condenses ~1600-1700K (very hot)
        if(CONDENSIBLES[i]==107) psat[i] = ms_psat_al2O3(tl[lay]); // Condenses ~1700-1800K (very hot)
        */


        // Start condensation calculation
        Xc[i] = clouds[lay][CONDENSIBLES[i]]/MM[lay]; //preexisting clouds taken into account
        Xv[i] = xx[lay][CONDENSIBLES[i]]/MM[lay]; //gas fraction of condensibles

        // Calculate total available condensible including condensed stuff
        double Xtotal = Xv[i] + Xc[i];  // Total condensible (gas + cloud)

        // Calculate equilibrium partitioning
        double Xv_sat = psat[i] / pl[lay];  // Maximum gas phase at saturation
        double Xc_equilibrium = fmax(0.0, Xtotal - Xv_sat);  // Equilibrium cloud amount

        // CORRECT: Update phases while conserving total mass
        if (Xtotal > Xv_sat) {
            // Supersaturated: some should be cloud
            Xv[i] = Xv_sat;
            Xc[i] = Xc_equilibrium;
        } else {
            // Undersaturated: all should be gas  
            Xv[i] = Xtotal;
            Xc[i] = 0.0;
        }

        Xd -= Xv[i] + Xc[i];

        // //calculate condensation rate for each condensible species (super-saturation)
        // dp[i] = pl[lay]*Xv[i] - psat[i];
        // //debugging if(lay==1 || dp[i]>0.0) printf("%s%d%s%e\t%s%d%s%e\n","pl[",lay,"] = ",pl[lay],"psat[",i,"] = ",psat[i]);
       
        // if(dp[i]>=0.0) //if saturated, then add to clouds and reduce gas fraction
        // {
        //     Xc[i] += (1.0 - psat[i]/(pl[lay]*Xv[i]))*Xv[i]; //condensed fraction
        //     Xvold = Xv[i]; //copy old vapor fraction
        //     // adjust vapor fraction to be in equilibrium with the saturation pressure
        //     Xv[i] = psat[i]/(pl[lay]); // here for surface layer: assumes "infinite" ocean reservoir for all condensibles (changed from Xv[i] *= psat[i]/(pl[lay]*Xv[i]) to Xv[i] = psat[i]/pl[lay]))
        //     //printf("%s%d\t%s%d%s%e%s%e\t%s%d%s%e\n","Saturation in layer ",lay," Xv[",CONDENSIBLES[i],"] was ",Xvold, " now ",Xv[i]," Xc[",CONDENSIBLES[i],"] = ",Xc[i]);
        // }
        // else //ensure unsaturated treatment, i.e. dry adiabat, without any latent heat release is applied
        // {
        //     //Old implementation with reset
        //     //Xc[i] = 0.0;
        //     //Xv[i] = 0.0; this was originally in the code, but I dont think this reset is needed
        //     // Xv is used for lapse rate calculation
        //     // MASS CONSERVATION FIX: Add cloud mass back to vapor phase before removing the cloud
        //     if (Xc[i] > 0.0) {
        //         // Only log significant evaporation events
        //         if (Xc[i] > 1.0e-10) {
        //             printf("Cloud evaporation: Layer %d, Species %d, Cloud VMR %.2e being returned to gas phase\n", 
        //                    lay, CONDENSIBLES[i], Xc[i]);
        //             printf("  Before evaporation: Xv[%d] = %.6e, xx[%d][%d] = %.6e\n", 
        //                    i, Xv[i], lay, CONDENSIBLES[i], xx[lay][CONDENSIBLES[i]]);
        //         }
        //         // Add the cloud mass back to vapor phase
        //         Xv[i] += Xc[i];
        //         // Reset cloud to zero
        //         Xc[i] = 0.0;
                
        //         if (Xc[i] > 1.0e-10) {
        //             printf("  After adding cloud to vapor: Xv[%d] = %.6e\n", i, Xv[i]);
        //         }
        //     } else {
        //         // No cloud mass to evaporate
        //         Xc[i] = 0.0;
        //     }
        // }
        // remaining dry mol fraction is equal to 1 - sum of all vapor and condensed fractions

    //if(lay==1 || dp[i]>0.0) printf("%s%e\n","Xd = ",Xd);
    }
    
    //*=========Heat capacities, condensate retention, latent heat=============*//
    //=======================================================================*//

    // CORRECTED APPROACH: Calculate heat capacities separately for each phase
    // following Graham+2021 methodology exactly
    
    // STEP 1: Calculate heat capacity of TRULY DRY (non-condensible) species only
    double cpxx_dry = 0.0;  // Heat capacity contribution from dry species
    double MM_dry = 0.0;    // Number density of dry species
    
    // Add helium (always dry)
    cpxx_dry += xxHe * HeHeat(tl[lay]);
    MM_dry += xxHe;
    
    // Add other species that are NOT in the condensibles list
    // This ensures we only count truly non-condensible species as "dry"
    
    // Check each species to see if it's condensible
    bool is_condensible[NSP+1] = {false}; // Initialize all to false
    for (i=0; i<NCONDENSIBLES; i++) {
        is_condensible[CONDENSIBLES[i]] = true;
    }
    
    // Add non-condensible species to dry calculation
    int species_list[] = {7, 9, 20, 21, 45, 52, 53, 54, 55}; // All species with heat capacity functions
    int num_species = 9;
    
    for (int s=0; s<num_species; s++) {
        int species_id = species_list[s];
        if (!is_condensible[species_id]) {
            // This species is truly dry (non-condensible)
            switch(species_id) {
                case 7:  cpxx_dry += xx[lay][7] * H2OHeat(tl[lay]); MM_dry += xx[lay][7]; break;
                case 9:  cpxx_dry += xx[lay][9] * NH3Heat(tl[lay]); MM_dry += xx[lay][9]; break;
                case 20: cpxx_dry += xx[lay][20] * COHeat(tl[lay]); MM_dry += xx[lay][20]; break;
                case 21: cpxx_dry += xx[lay][21] * CH4Heat(tl[lay]); MM_dry += xx[lay][21]; break;
                case 45: cpxx_dry += xx[lay][45] * H2SHeat(tl[lay]); MM_dry += xx[lay][45]; break;
                case 52: cpxx_dry += xx[lay][52] * CO2Heat(tl[lay]); MM_dry += xx[lay][52]; break;
                case 53: cpxx_dry += xx[lay][53] * H2Heat(tl[lay]); MM_dry += xx[lay][53]; break;
                case 54: cpxx_dry += xx[lay][54] * O2Heat(tl[lay]); MM_dry += xx[lay][54]; break;
                case 55: cpxx_dry += xx[lay][55] * N2Heat(tl[lay]); MM_dry += xx[lay][55]; break;
            }
        }
    }
    
    // STEP 2: Calculate average heat capacity of dry species
    if (MM_dry > 0.0) {
        cp_d = cpxx_dry / MM_dry;
    } else {
        // No dry species - use a reasonable default (this shouldn't happen normally)
        cp_d = 29.0; // Approximate cp for diatomic gases
        printf("WARNING: No dry species found in layer %d, using default cp_d = %.1f\n", lay, cp_d);
    }
    
    // STEP 3: Calculate heat capacities for condensible species
    for (i=0; i<NCONDENSIBLES; i++)
    {
        // Set vapor and condensed heat capacities for each condensible species
        if(CONDENSIBLES[i]==7)  {
            cp_v[i] = H2OHeat(tl[lay]); 
            cp_c[i] = H2O_liquid_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==9)  {
            cp_v[i] = NH3Heat(tl[lay]); 
            cp_c[i] = NH3_liquid_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==20) {
            cp_v[i] = COHeat(tl[lay]);  
            cp_c[i] = CO_condensed_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==21) {
            cp_v[i] = CH4Heat(tl[lay]); 
            cp_c[i] = CH4_condensed_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==52) {
            cp_v[i] = CO2Heat(tl[lay]); 
            cp_c[i] = CO2_solid_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==53) {
            cp_v[i] = H2Heat(tl[lay]);  
            cp_c[i] = H2_condensed_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==54) {
            cp_v[i] = O2Heat(tl[lay]);  
            cp_c[i] = O2_condensed_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==55) {
            cp_v[i] = N2Heat(tl[lay]);  
            cp_c[i] = N2_condensed_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else if(CONDENSIBLES[i]==45) {
            cp_v[i] = H2SHeat(tl[lay]);  
            cp_c[i] = H2S_condensed_heat_capacity(tl[lay]);
            alpha[i] = ALPHA_RAINOUT;
        }
        else {
            // Unknown condensible species - use defaults
            cp_v[i] = 30.0;  // Default vapor heat capacity
            cp_c[i] = 30.0;  // Default condensed heat capacity  
            alpha[i] = ALPHA_RAINOUT;
            printf("WARNING: Unknown condensible species %d in layer %d\n", CONDENSIBLES[i], lay);
        }
    
        // Calculate latent heat of the condensible species
        latent[i] = ms_latent(CONDENSIBLES[i],tl[lay]);
        
        // CORRECT APPROACH - following JANUS logic exactly
        // β is ONLY applied when condensation is actively occurring (P ≥ P_sat)
        // When P < P_sat, species is treated as dry (β = 0, no latent heat effect)
        double partial_pressure = Xv[i] * pl[lay];
        
        // Get critical temperature for this species
        double T_crit = 1e6;  // Default very high temperature (no condensation)
        if(CONDENSIBLES[i]==7)  T_crit = 647.1;   // H2O critical temperature
        if(CONDENSIBLES[i]==9)  T_crit = 405.5;   // NH3 critical temperature
        if(CONDENSIBLES[i]==20) T_crit = 132.92;  // CO critical temperature
        if(CONDENSIBLES[i]==21) T_crit = 190.44;  // CH4 critical temperature
        if(CONDENSIBLES[i]==52) T_crit = 304.2;   // CO2 critical temperature
        if(CONDENSIBLES[i]==53) T_crit = 33.2;    // H2 critical temperature
        if(CONDENSIBLES[i]==54) T_crit = 154.54;  // O2 critical temperature
        if(CONDENSIBLES[i]==55) T_crit = 126.2;   // N2 critical temperature
        if(CONDENSIBLES[i]==45) T_crit = 373.3;   // H2S critical temperature
        
        if (partial_pressure >= psat[i] && tl[lay] < T_crit) {
            // CONDENSATION OCCURRING: Apply positive β (latent heat release)
            beta[i] = latent[i] / R_GAS / tl[lay];
        } else {
            // NO CONDENSATION: Treat as dry species (β = 0)
            // This includes: P < P_sat (undersaturated) OR T > T_critical
            beta[i] = 0.0;
        }
       // This is an example of a possible physics based implementation with fall velocity
        // can also do based on pressure.
        // // PHYSICALLY-BASED RAINOUT CALCULATION
        // // Calculate retention fraction alpha[i] based on fall velocity vs mixing time
        // if (Xc[i] > 0.0) {
        //     // Typical droplet/crystal sizes and fall velocities
        //     double particle_radius, fall_velocity, mixing_time, fall_time;
        //     double layer_thickness = 1000.0; // meters (typical layer thickness)
        //     double atmospheric_density = pl[lay] / (R_GAS * tl[lay]) * 0.001; // kg/m³
            
        //     if (CONDENSIBLES[i] == 7) { // H2O
        //         if (tl[lay] < 273.0) {
        //             // Ice crystals
        //             particle_radius = 50e-6;  // 50 microns
        //             fall_velocity = 0.5;      // m/s (typical for ice crystals)
        //         } else {
        //             // Liquid droplets  
        //             particle_radius = 10e-6;  // 10 microns
        //             fall_velocity = 0.01;     // m/s (small droplets fall slowly)
        //         }
        //     } else if (CONDENSIBLES[i] == 52) { // CO2 (dry ice)
        //         particle_radius = 100e-6;     // 100 microns
        //         fall_velocity = 1.0;          // m/s (dense CO2 ice)
        //     } else {
        //         // Default for other species
        //         particle_radius = 20e-6;      // 20 microns
        //         fall_velocity = 0.1;          // m/s
        //     }
            
        //     // Adjust fall velocity for atmospheric density
        //     fall_velocity *= sqrt(atmospheric_density / 1.225); // Scale with air density
            
        //     // Calculate time scales
        //     fall_time = layer_thickness / fall_velocity;      // Time to fall through layer
        //     mixing_time = 3600.0;  // Assumed convective mixing time (1 hour)
            
        //     // Retention fraction: what stays vs what falls out
        //     alpha[i] = mixing_time / (mixing_time + fall_time);
            
        //     // Apply bounds
        //     if (alpha[i] < 0.01) alpha[i] = 0.01;  // Minimum 1% retention
        //     if (alpha[i] > 0.99) alpha[i] = 0.99;  // Maximum 99% retention
            
        //     if (lay==1) printf("Rainout calc - Species %d: radius=%.1e m, vfall=%.3f m/s, tfall=%.0f s, alpha=%.3f\n",
        //            CONDENSIBLES[i], particle_radius, fall_velocity, fall_time, alpha[i]);
        // } else {
        //     alpha[i] = 1.0;  // No condensate, so retention is irrelevant
        // }
        
        //if(lay==1)  printf("  After: cpxx=%.3e, MM_cp=%.3e\n", cpxx, MM_cp);
    }
    
    // // Debug output for heat capacity calculation
    // if(lay==1) {
    //     printf("HEAT CAPACITY CALCULATION - Layer %d:\n", lay);
    //     printf("  Dry species: MM_dry = %.3e, cpxx_dry = %.3e, cp_d = %.3f J/(mol·K)\n", 
    //            MM_dry, cpxx_dry, cp_d);
    //     for (i=0; i<NCONDENSIBLES; i++) {
    //         printf("  Species %d: cp_v = %.3f, cp_c = %.3f J/(mol·K)\n", 
    //                CONDENSIBLES[i], cp_v[i], cp_c[i]);
    //     }
    // }

    // EQ 1 from Graham+2021 d ln T / d ln P = (x_d + Σx_{v,i}) / 
                //(x_d * [c_d*x_d + Σ(x_{v,i}*(c_{v,i} - R*β_i + R*β_i²) + α_i*x_{c,i}*c_{c,i})] / 
                // [R*(x_d + Σβ_i*x_{v,i})] + Σβ_i*x_{v,i})
    
    // Calculate adiabatic lapse rate
    lapse_num = Xd;
    sum_beta_xv = 0.0;
    big_sum_denom_num_left_term = cp_d * Xd;
    big_sum_denom_num_right_term = 0.0;
    cp_num = cp_d*Xd;
    cp_denom = Xd;

    for (i=0; i<NCONDENSIBLES; i++) 
    {
        lapse_num += Xv[i];
        sum_beta_xv += beta[i]*Xv[i];
        big_sum_denom_num_right_term += Xv[i]*(cp_v[i] - R_GAS*beta[i] + R_GAS*beta[i]*beta[i]) + alpha[i]*Xc[i]*cp_c[i];
        cp_num +=  Xv[i]*cp_v[i] + alpha[i]*Xc[i]*cp_c[i];
        cp_denom += Xv[i];  // CORRECT: Graham's paper denominator is (x_d + Σx_{v,i}) only
        

        // if (Xc[i] > 1.0e-10 && lay==120) {  // Debug output for any condensation, regardless of layer
        //     printf("HEAT CAPACITY DEBUG - Layer %d, Species %d:\n", lay, CONDENSIBLES[i]);
        //     printf("  Xv[%d] = %.6e, Xc[%d] = %.6e, alpha[%d] = %.3f\n", i, Xv[i], i, Xc[i], i, alpha[i]);
        //     printf("  cp_v[%d] = %.3f J/(mol·K), cp_c[%d] = %.3f J/(mol·K)\n", i, cp_v[i], i, cp_c[i]);
        //     printf("  Vapor contribution to cp_num: %.6e\n", Xv[i]*cp_v[i]);
        //     printf("  Condensate contribution to cp_num: %.6e\n", alpha[i]*Xc[i]*cp_c[i]);
        //     printf("  Running cp_num = %.6e, cp_denom = %.6e\n", cp_num, cp_denom);
        //     printf("  Current effective cp = %.3f J/(mol·K)\n", cp_num/cp_denom);
        // }
    }

    lapse_denom = Xd * (big_sum_denom_num_left_term + big_sum_denom_num_right_term) / (R_GAS*(Xd + sum_beta_xv)) + sum_beta_xv;

    lapse[lay] = lapse_num / lapse_denom;

    *cp = cp_num / cp_denom; //return to climate module

    //if(lay==1) printf("%s%e\n","lapse_num = ",lapse_num);
    //if(lay==1) printf("%s%e\n","lapse_denom = ",lapse_denom);
    //if(lay==1) printf("%s%e\n","sum_beta_xv = ",sum_beta_xv);
    //if(lay==1) printf("%s%e\n","big_sum_denom_num = ",big_sum_denom_num);
//    if(lay==1) printf("%s%e\n","lapse = ",lapse[lay]);

    //propagate compositional changes back to main program 
    // this is where the rainout needs to happen
    
    // NO RAINOUT HERE - rainout is now handled separately in ms_rainout()
    // Just set the gas and cloud phases based on current equilibrium
    for (i=0; i<NCONDENSIBLES; i++) {
        xx[lay][CONDENSIBLES[i]] = Xv[i] * MM[lay];
        clouds[lay][CONDENSIBLES[i]] = Xc[i] * MM[lay];
        
        // Calculate and store saturation ratios for plotting
        // saturation_ratio = partial_pressure / saturation_pressure
        double partial_pressure = Xv[i] * pl[lay];
        if (psat[i] > 0.0) {
            saturation_ratios[i] = partial_pressure / psat[i];
        } else {
            saturation_ratios[i] = 0.0;
        }
    }
    
    // if (lay == zbin && Xv[0] > 0.1) {  // Check only top layer with significant H2O
    //     printf("LAPSE RATE DEBUG - Layer %d:\n", lay);
    //     printf("  Xd = %.6f, Xv[H2O] = %.6f, Xc[H2O] = %.6f\n", Xd, Xv[0], Xc[0]);
    //     printf("  beta[H2O] = %.6f, latent[H2O] = %.3e J/mol\n", beta[0], latent[0]);
    //     printf("  sum_beta_xv = %.6f\n", sum_beta_xv);
    //     printf("  lapse_num = %.6f\n", lapse_num);
    //     printf("  lapse_denom = %.6f\n", lapse_denom);
    //     printf("  Final lapse rate = %.6f\n", lapse[lay]);
    //     printf("  Temperature = %.1f K, Pressure = %.3e Pa\n", tl[lay], pl[lay]);
    //     printf("\n");
    // }
}// END: void ms_adiabat()
//*********************************************************
//*********************************************************
void ms_rainout(int lay, double* mass_loss_ratio)
{
    // Apply rainout physics once per RC iteration
    // This function handles the Graham+2021 rainout approach with proper mass conservation
    
    int i, j;
    
    // Get cloud and vapor amounts for condensible species
    double Xv[NCONDENSIBLES], Xc[NCONDENSIBLES];
    double alpha[NCONDENSIBLES]; // retention factors (must match ms_adiabat values)
    
    for (i=0; i<NCONDENSIBLES; i++) {
        Xv[i] = xx[lay][CONDENSIBLES[i]]/MM[lay]; // vapor phase mole fraction
        Xc[i] = clouds[lay][CONDENSIBLES[i]]/MM[lay]; // cloud phase mole fraction
        
        // Use config-defined alpha values - consistent with ms_adiabat
        alpha[i] = ALPHA_RAINOUT;
    }
    
    // Store original cloud amounts before modification
    double Xc_original[NCONDENSIBLES];
    for (i=0; i<NCONDENSIBLES; i++) {
        Xc_original[i] = Xc[i];
        //printf("RAINOUT: Layer %d, Species %d, Xc_original[%d] = %.6e\n", 
        //       lay, CONDENSIBLES[i], i, Xc_original[i]);
    }

    // STEP 1: Update condensed species amounts after rainout
    for (i=0; i<NCONDENSIBLES; i++) {
        // adjust condensable fraction with alpha
        Xc[i] = alpha[i] * Xc[i];
        xx[lay][CONDENSIBLES[i]] = Xv[i] * MM[lay];
        clouds[lay][CONDENSIBLES[i]] = Xc[i] * MM[lay];

        //printf("RAINOUT: Layer %d, Species %d, Xc[%d] = %.6e, Xv[%d] = %.6e, clouds[%d] = %.6e\n", 
        //       lay, CONDENSIBLES[i], i, Xc[i], i, Xv[i], i, clouds[lay][CONDENSIBLES[i]]);
    }

        
    // STEP 2: Update non-condensible species based on conservation approach
    if (PRESSURE_CONSERVATION == 1) {
        // REALISTIC APPROACH (Closed System):
        // Non-condensibles maintain absolute abundance, mole fractions increase naturally
        // No scaling needed - they're already physically correct
        // Atmospheric mass loss will be tracked for pressure adjustment
        
        // Calculate total mass loss for pressure adjustment
        double original_total_mass = 0.0;
        double new_total_mass = 0.0;
        double remaining_fraction = 0.0;
        // Sum all current species abundances
        for (j=1; j<=NSP; j++) {
            new_total_mass += xx[lay][j];
        }
        
        // Estimate original mass before rainout (approximate)
        original_total_mass = new_total_mass / remaining_fraction;
        
        // Store mass loss ratio for pressure adjustment (global variable or pass back)
        *mass_loss_ratio = new_total_mass / original_total_mass;
        
        if (lay==62) { // Debug output for top layer
            printf("REALISTIC RAINOUT: Layer %d - Mass loss ratio: %.6f\n", lay, *mass_loss_ratio);
            printf("  Original mass estimate: %.6e, New mass: %.6e\n", original_total_mass, new_total_mass);
        }
    } else {
        // OPEN SYSTEM APPROACH: Maintain Graham+2021 constraint x_d + x_v + x_c = 1.0
        // Layer-specific scaling: only species with clouds in THIS layer are exempt
        // All other species (including condensibles without clouds here) get scaled
        
        // Calculate total mole fraction lost from cloud rainout in THIS layer
        double total_mole_fraction_lost = 0.0;
        for (i=0; i<NCONDENSIBLES; i++) {
            total_mole_fraction_lost += (1.0 - alpha[i]) * Xc_original[i];
        }

        if (total_mole_fraction_lost > 0.0) {
            // Identify which species actually have clouds in THIS layer
            bool has_clouds_in_layer[NSP+1] = {false};  // Initialize all to false
            for (i=0; i<NCONDENSIBLES; i++) {
                if (clouds[lay][CONDENSIBLES[i]] > 1.0e-20) {  // Small threshold for numerical precision
                    has_clouds_in_layer[CONDENSIBLES[i]] = true;
                }
            }
            
            // Calculate current mole fraction of species that will NOT be scaled
            // (only those with actual clouds in this layer)
            double unscaled_mf = 0.0;
            for (j=1; j<=NSP; j++) {
                if (has_clouds_in_layer[j]) {
                    unscaled_mf += xx[lay][j]/MM[lay];
                }
            }
            
            // Calculate current mole fraction of species that WILL be scaled
            // (everything else, including condensibles without clouds here)
            double scalable_mf = 0.0;
            for (j=1; j<=NSP; j++) {
                if (!has_clouds_in_layer[j]) {
                    scalable_mf += xx[lay][j]/MM[lay];
                }
            }
            
            // Calculate new cloud mole fraction (after rainout)
            double new_cloud_mf = 0.0;
            for (i=0; i<NCONDENSIBLES; i++) {
                new_cloud_mf += clouds[lay][CONDENSIBLES[i]]/MM[lay];
            }
            
            // Scale factor: remaining space divided by what needs to fit
            double scale_factor = (1.0 - unscaled_mf - new_cloud_mf) / scalable_mf;
            
            if (lay==60) { // Debug output for mid-atmosphere layer
                printf("OPEN SYSTEM RAINOUT: Layer %d\n", lay);
                printf("  Cloud mole fraction lost: %.6e\n", total_mole_fraction_lost);
                printf("  Unscaled species mole fraction: %.6e\n", unscaled_mf);
                printf("  Scalable species mole fraction: %.6e\n", scalable_mf);
                printf("  New cloud mole fraction: %.6e\n", new_cloud_mf);
                printf("  Scale factor: %.6f\n", scale_factor);
                
                // Show which species have clouds in this layer
                printf("  Species with clouds in this layer:");
                for (j=1; j<=NSP; j++) {
                    if (has_clouds_in_layer[j]) {
                        printf(" %d", j);
                    }
                }
                printf("\n");
            }

            // Scale only species that do NOT have clouds in this layer
            for (j=1; j<=NSP; j++) {
                if (!has_clouds_in_layer[j]) {
                    xx[lay][j] *= scale_factor;
                }
                // Species with clouds in this layer remain unchanged
            }
            
            if (lay==60) { // Verification
                // Check that mole fractions sum to 1.0
                double total_mole_fraction = 0.0;
                
                // Gas phase mole fractions
                for (j=1; j<=NSP; j++) {
                    total_mole_fraction += xx[lay][j]/MM[lay];
                }
                
                // Condensed phase mole fractions (after rainout)
                for (i=0; i<NCONDENSIBLES; i++) {
                    total_mole_fraction += clouds[lay][CONDENSIBLES[i]]/MM[lay];
                }
                
                printf("  Verification: Total mole fraction = %.6f (should be 1.0)\n", total_mole_fraction);
                
                // Show breakdown
                double scaled_mf = 0.0, unscaled_mf_check = 0.0, cloud_mf = 0.0;
                for (j=1; j<=NSP; j++) {
                    if (has_clouds_in_layer[j]) {
                        unscaled_mf_check += xx[lay][j]/MM[lay];
                    } else {
                        scaled_mf += xx[lay][j]/MM[lay];
                    }
                }
                for (i=0; i<NCONDENSIBLES; i++) {
                    cloud_mf += clouds[lay][CONDENSIBLES[i]]/MM[lay];
                }
                printf("  Breakdown: Scaled=%.6f, Unscaled=%.6f, Clouds=%.6f\n", 
                       scaled_mf, unscaled_mf_check, cloud_mf);
            }
        }
        
        // For open system, no pressure adjustment needed (pressure stays fixed)
        *mass_loss_ratio = 1.0;
    }
    
    // NOTE: Enhanced cloud physics is now handled separately
    // This function only handles Graham's basic alpha rainout
    // Enhanced cloud physics and sedimentation are handled in separate functions
    
}// END: void ms_rainout()

//=========================================================
//=== Enhanced Cloud Physics Function ====================
//=========================================================

// SPECIES-SPECIFIC PARTICLE PROPERTIES FUNCTION
// Returns particle density [kg/m³], accommodation coefficient [dimensionless], and molecular mass [AMU]
// based on condensible species ID and temperature
void get_particle_properties(int species_id, double temperature, double *density, double *accommodation_coeff, double *molecular_mass) {
    switch (species_id) {
        case 7: // H₂O (Water)
            *density = (temperature < 273.0) ? 917.0 : 1000.0; // Ice vs liquid water [kg/m³]
            *accommodation_coeff = 0.8; // Accommodation coefficient for water vapor [dimensionless]
            *molecular_mass = 18.015; // Molecular mass [AMU]
            break;
            
        case 9: // NH₃ (Ammonia)
            *density = (temperature < 195.0) ? 817.0 : 682.0; // Solid vs liquid ammonia [kg/m³]
            *accommodation_coeff = 0.9; // High sticking probability for NH₃ [dimensionless]
            *molecular_mass = 17.031; // Molecular mass [AMU]
            break;
            
        case 21: // CH₄ (Methane)
            *density = (temperature < 90.7) ? 519.0 : 422.0; // Solid vs liquid methane [kg/m³]
            *accommodation_coeff = 0.7; // Moderate sticking for CH₄ [dimensionless]
            *molecular_mass = 16.043; // Molecular mass [AMU]
            break;
            
        case 52: // CO₂ (Carbon dioxide)
            *density = (temperature < 216.5) ? 1563.0 : 1178.0; // Dry ice vs liquid CO₂ [kg/m³]
            *accommodation_coeff = 0.85; // Good sticking for CO₂ [dimensionless]
            *molecular_mass = 44.010; // Molecular mass [AMU]
            break;
            
        case 45: // H₂S (Hydrogen sulfide)
            *density = (temperature < 187.7) ? 993.0 : 949.0; // Solid vs liquid H₂S [kg/m³]
            *accommodation_coeff = 0.9; // High sticking for H₂S [dimensionless]
            *molecular_mass = 34.081; // Molecular mass [AMU]
            break;
            
        case 20: // CO (Carbon monoxide)
            *density = (temperature < 68.0) ? 929.0 : 789.0; // Solid vs liquid CO [kg/m³]
            *accommodation_coeff = 0.6; // Lower sticking for light molecules [dimensionless]
            *molecular_mass = 28.010; // Molecular mass [AMU]
            break;
            
        case 55: // N₂ (Nitrogen)
            *density = (temperature < 63.1) ? 1026.0 : 808.0; // Solid vs liquid N₂ [kg/m³]
            *accommodation_coeff = 0.6; // Lower sticking for light molecules [dimensionless]
            *molecular_mass = 28.014; // Molecular mass [AMU]
            break;
            
        case 54: // O₂ (Oxygen)
            *density = (temperature < 54.4) ? 1426.0 : 1141.0; // Solid vs liquid O₂ [kg/m³]
            *accommodation_coeff = 0.6; // Lower sticking for light molecules [dimensionless]
            *molecular_mass = 31.998; // Molecular mass [AMU]
            break;
            
        case 53: // H₂ (Hydrogen)
            *density = (temperature < 14.0) ? 86.8 : 70.8; // Solid vs liquid H₂ [kg/m³]
            *accommodation_coeff = 0.5; // Very low sticking for H₂ [dimensionless]
            *molecular_mass = 2.016; // Molecular mass [AMU]
            break;
            
        default: // Unknown species - use water-like properties
            *density = 1000.0; // Default to liquid water density [kg/m³]
            *accommodation_coeff = 0.8; // Default accommodation coefficient [dimensionless]
            *molecular_mass = 18.0; // Default molecular mass [AMU]
            break;
    }
}

//particlesizef function from Kyo
// PARTICLE SIZE CALCULATION USING MICROPHYSICAL BALANCE
// This function calculates equilibrium particle sizes by balancing:
// 1) Condensational growth (mass diffusion from supersaturated vapor)
// 2) Gravitational settling (particles fall due to gravity)
// 3) Eddy diffusion (turbulent mixing opposes settling)
//
// Core Physics: Particles grow until fall velocity balances with diffusion
// Key Equation: u = v_fall * H / K_zz (dimensionless fall parameter)
// Where: H = scale height [m], K_zz = eddy diffusion coefficient [m²/s]
void particlesizef_local(double g, double T, double P, double mean_molecular_mass, int condensible_species_id, double Kzz, double deltaP, double *r0, double *r1, double *r2, double *VP, double *effective_settling_velocity, double *scale_height, double *retention_factor) {
    // Local constants to avoid conflicts with main code
    double KB_local = 1.3806503E-23; // Boltzmann constant [J/K]
    double AMU_local = 1.66053886E-27; // Atomic mass unit [kg]
    double RGAS_local = 8.314472; // Universal gas constant [J/(mol·K)]
    double PI_local = 3.1416; // Pi

    // GET SPECIES-SPECIFIC PROPERTIES
    double molecular_mass_condensible; // [AMU]
    double rho, acc; // Particle density [kg/m³] and accommodation coefficient [dimensionless]
    
    // Get all species-specific properties from unified function
    get_particle_properties(condensible_species_id, T, &rho, &acc, &molecular_mass_condensible);

    // ATMOSPHERIC SCALE HEIGHT CALCULATION
    // H = k_B * T / (m_avg * g) where m_avg = mean_molecular_mass * AMU
    // Units: [J/K] * [K] / ([AMU] * [kg/AMU] * [m/s²]) = [m]
    double H = KB_local * T / mean_molecular_mass / AMU_local / g; // Atmospheric scale height [m]
    
    // DIMENSIONLESS FALL PARAMETER
    // u = K_zz / H represents ratio of diffusion to settling
    // Higher u = more diffusion relative to settling = smaller particles
    double u = Kzz / H; // Dimensionless eddy diffusion parameter
    
    // ATMOSPHERIC VISCOSITY (temperature-dependent)
    // Sutherland's formula: μ = μ₀ * (T/T₀)^1.5 * (T₀ + S)/(T + S)
    
    // ATMOSPHERIC COMPOSITION SELECTION
    // Choose atmospheric type: 0=H2, 1=Air, 2=CO2, 3=N2, 4=Ar, 5=He
    int atmosphere_type = 0;  // DEFAULT: H2 atmosphere (common for gas giants/exoplanets)
    
    // Sutherland constants for different atmospheric compositions
    // Values from Sutherland (1893), COMSOL, and atmospheric physics literature
    double mu0, T0, S_sutherland;
    
    if (atmosphere_type == 0) {
        // HYDROGEN (H2) - Default for gas giant/exoplanet atmospheres
        mu0 = 8.411E-6;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 97.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 1) {
        // AIR (N2/O2 mix) - Earth-like atmospheres
        mu0 = 1.716E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K] 
        S_sutherland = 110.4;  // Sutherland constant [K]
    } else if (atmosphere_type == 2) {
        // CARBON DIOXIDE (CO2) - Venus-like or early Mars
        mu0 = 1.370E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 222.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 3) {
        // NITROGEN (N2) - Titan-like atmospheres
        mu0 = 1.663E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 107.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 4) {
        // ARGON (Ar) - Some planetary atmospheres
        mu0 = 2.125E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 114.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 5) {
        // HELIUM (He) - Helium-rich atmospheres
        mu0 = 1.865E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 79.4;   // Sutherland constant [K]
    } else {
        // DEFAULT: Fall back to H2 if unknown type
        mu0 = 8.411E-6;    // H2 values
        T0 = 273.15;
        S_sutherland = 97.0;
    }
    
    // Calculate atmospheric viscosity using Sutherland's formula
    double mu = mu0 * pow(T / T0, 1.5) * (T0 + S_sutherland) / (T + S_sutherland); // Viscosity [Pa·s]
    
    // MEAN FREE PATH OF GAS MOLECULES
    // λ = 2μ / (P * √(8*mean_molecular_mass/(πRT))) - kinetic theory result
    // Units: [Pa·s] / ([Pa] * √([kg/mol]/([J/(mol·K)] * [K]))) = [m]
    double lambda = 2 * mu / P / pow(8 * mean_molecular_mass * 1.0E-3 / PI_local / RGAS_local / T, 0.5); // Mean free path [m]
    
    // NUMBER DENSITY OF CONDENSING MOLECULES
    // From ideal gas law: n = ΔP / (k_B * T)
    // where ΔP = P_partial - P_saturation (supersaturation pressure)
    // Units: [Pa] / ([J/K] * [K]) = [molecules/m³]
    double deltan = deltaP / KB_local / T; // Excess number density [molecules/m³]

    // MASS DIFFUSION COEFFICIENT
    // Approximation for molecular diffusion in gas
    double D = 0.12e-4; // Mass diffusion coefficient [m²/s]

    // PARTICLE SIZE DISTRIBUTION PARAMETERS
    double Cc0= 1.0; // Initial Cunningham slip correction factor [dimensionless]
    double fa = 1.0; // Initial ventilation factor [dimensionless]  
    double sig = 2.0; // Log-normal size distribution width parameter [dimensionless]

    // ITERATIVE SOLUTION FOR EQUILIBRIUM PARTICLE SIZE
    // CORRECTED IMPLEMENTATION FOLLOWING HU+2019 METHODOLOGY
    // Solve for particle volume V by balancing:
    // 1) Condensational growth: ∝ D * deltan / ρ
    // 2) Effective gravitational settling: ∝ max[settling_term - u, 0]  
    // 3) Turbulent diffusion: ∝ u / H
    double Vs;
    for (int dump = 1; dump <= 1e3; ++dump) {
        // CONDENSATION TERM: Rate of volume growth from vapor deposition
        // cc = -48^(1/3) * π^(2/3) * D * molecular_mass_condensible * fa * Δn / ρ * exp(-ln²(σ))
        // Negative because we're solving the quadratic equation
        // This term should be negligible if deltan is low
        double cc = -(pow(48.0 * PI_local * PI_local, 1.0 / 3.0) * D * molecular_mass_condensible * AMU_local * fa * deltan / rho * exp(-pow(log(sig), 2.0)));
        
        // SETTLING TERM: Gravitational fall velocity coefficient  
        // aa = ρ * g / (μ * 162^(1/3) * π^(2/3) * H) * Cc * exp(-ln²(σ))
        // From Stokes law: v_fall = 2*r²*ρ*g*Cc/(9*μ) where r ∝ V^(1/3)
        // This represents the settling velocity term before the max[... - u, 0] operation
        double aa = rho * g / mu / pow(162.0 * PI_local * PI_local, 1.0 / 3.0) / H * Cc0* exp(-pow(log(sig), 2.0));
        
        // DIFFUSION TERM: Turbulent mixing coefficient
        // bb = -u/H = -K_zz/H² (opposes settling)
        double bb = -u / H;

        // SOLVE QUADRATIC EQUATION: aa*V^(2/3) + bb*V^(1/3) + cc = 0
        // Using quadratic formula after substituting V^(1/3) = x
        // This gives equilibrium volume where growth = settling + diffusion
        double V = pow((-bb + sqrt(bb * bb - 4.0 * aa * cc)) / 2.0 / aa, 3.0 / 2.0);
        
        // Handle cases where updraft dominates (no real solution)
        if (isnan(V) || V <= 0.0) {
            // Updraft-dominated regime: use asymptotic solution
            //V_asym = 39.9 * [μ*u*exp(ln²σ)/(ρ_p*g*C_c)]^(3/2)
            V = 39.9 * pow(mu * u * exp(pow(log(sig), 2.0)) / (rho * g * Cc0), 3.0 / 2.0);
        }
        
        // PARTICLE DIAMETER from volume
        // d = (6V/π)^(1/3) * exp(-ln²(σ)) [m]
        double d1 = pow(6.0 * V / PI_local, 1.0 / 3.0) * exp(-pow(log(sig), 2.0));

        // UPDATE SLIP CORRECTION FACTORS
        // Knudsen number: Kn = λ/d (ratio of mean free path to particle size)
        double kn = lambda / d1;
        
        // CUNNINGHAM SLIP CORRECTION: Cc = 1 + Kn*(1.257 + 0.4*exp(-1.1/Kn))
        // Accounts for non-continuum effects when particles approach molecular size
        double Cc1 = 1.0 + kn * (1.257 + 0.4 * exp(-1.1 / kn));
        
        // VENTILATION CORRECTION: Accounts for enhanced mass transfer during settling
        // fa = (1 + Kn)/(1 + 2*Kn*(1+Kn)/α) where α is accommodation coefficient
        double fa1 = (1.0 + kn) / (1.0 + 2.0 * kn * (1.0 + kn) / acc);
        
        // CHECK CONVERGENCE: Continue iteration until slip factors stabilize
        if (fabs(Cc1 - Cc0) + fabs(fa1 - fa) < 0.001) {
            Vs = V; // Final equilibrium volume [m³]
            break;
        } else {
            Cc0= Cc1; // Update slip correction
            fa = fa1;  // Update ventilation factor
        }
    }
    
    // CALCULATE CHARACTERISTIC PARTICLE RADII from log-normal distribution
    // For log-normal distribution with geometric standard deviation σ:
    // r_i = r_mean * exp(n * ln²(σ)) where n depends on moment
    
    //r0, r1, r2, VP are in microns
    // r₀: Mean radius weighted by number (smallest particles) [μm]
    *r0 = pow((3.0 * Vs) / (4.0 * PI_local), 1.0 / 3.0) * exp(-1.5 * pow(log(sig), 2.0)) * 1.0E+6;
    
    // r₁: Mean radius weighted by surface area [μm] 
    *r1 = pow((3.0 * Vs) / (4.0 * PI_local), 1.0 / 3.0) * exp(-pow(log(sig), 2.0)) * 1.0E+6;
    
    // r₂: Mean radius weighted by volume (largest/most massive particles) [μm]
    // This is used for sedimentation calculations since fall velocity ∝ r²
    *r2 = pow((3.0 * Vs) / (4.0 * PI_local), 1.0 / 3.0) * exp(-0.5 * pow(log(sig), 2.0)) * 1.0E+6;
    
    // VP: Total particle volume [μm³]
    *VP = Vs * 1.0E+6;
    
    // CALCULATE ADDITIONAL PHYSICS PARAMETERS FOR SEDIMENTATION
    // Use r2 (volume-weighted radius) for sedimentation calculations
    double particle_radius_m = *r2 * 1.0e-6; // Convert μm to m
    
    // Calculate gravitational settling velocity using final particle size
    // v_fall = 2*r²*ρ*g*Cc/(9*μ) - Stokes law with Cunningham correction
    double gravitational_settling = (2.0 * particle_radius_m * particle_radius_m * rho * g * Cc0) / (9.0 * mu);
    
    // This is the key correction from Hu+2019: v_d = max[v_fall - u, 0]
    // When updraft dominates, particles are carried upward (v_d = 0)
    // When settling dominates, particles fall with reduced velocity (v_d = v_fall - u)
    // double v_d = fmax(0.0, gravitational_settling - u);
    double v_d = gravitational_settling;
    
    // Calculate sedimentation parameter using effective settling velocity
    // f_sed = v_d * H / K_zz where v_d is the effective settling velocity
    double f_sed = (v_d * H) / Kzz;
    
    // Calculate retention factor using Ackerman & Marley (2001) formula
    // α = 1 / (1 + f_sed) - higher f_sed means more material falls out
    double retention = 1.0 / (1.0 + f_sed);
    
    // Apply physical bounds to retention
    if (retention < 0.001) retention = 0.001; // Minimum 0.1% retention
    if (retention > 0.95) retention = 0.95;   // Maximum 95% retention
    
    // Return additional physics parameters
    *effective_settling_velocity = v_d; // [m/s]
    *scale_height = H; // [m]
    *retention_factor = retention; // [dimensionless]

    printf ("VALUES FROM PARTICLESIZEF_LOCAL\n");
    printf ("viscosity: %.2e\n", mu);
    printf ("retention: %.2e\n", retention);
    printf ("gravitational_settling: %.2e\n", gravitational_settling);
    printf ("f_sed: %.2e\n", f_sed);
    printf ("scale_height: %.2e\n", H);
}

void apply_enhanced_cloud_physics(int layer, double gravity)
{
    // ENHANCED CLOUD PHYSICS: ACKERMAN & MARLEY (2001) + HU+2019
    // Uses actual supersaturation (deltaP) for particle growth
    // If Graham+2021 says all vapor is condensed, deltaP = 0 (or very small)
    // Particle size is set by microphysics, but will not grow if no supersaturation
    // Cloud profile is set by exponential decay from the cloud base
    // Total cloud mass never exceeds Graham+2021 equilibrium value

    int enhanced_cloud_physics_enabled = 1;
    if (!enhanced_cloud_physics_enabled) return;

    double Kzz_cm = 1.0e8;  // Eddy diffusion coefficient [cm^2/s]
    double Kzz = Kzz_cm * 1.0e-4;  // Eddy diffusion coefficient [m^2/s]

    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        double cloud_vmr = clouds[layer][species_id] / MM[layer];
        double vapor_vmr = xx[layer][species_id] / MM[layer];
        double psat = 0.0;
        switch (species_id) {
            case 7:  psat = ms_psat_h2o(tl[layer]); break; // H₂O (Water)
            case 9:  psat = ms_psat_nh3(tl[layer]); break; // NH₃ (Ammonia)  
            case 21: psat = ms_psat_ch4(tl[layer]); break; // CH₄ (Methane)
            case 52: psat = ms_psat_co2(tl[layer]); break; // CO₂ (Carbon dioxide)
            case 20: psat = ms_psat_co(tl[layer]);  break; // CO (Carbon monoxide)
            case 45: psat = ms_psat_h2s(tl[layer]); break; // H₂S (Hydrogen sulfide)
            case 53: psat = ms_psat_h2(tl[layer]);  break; // H₂ (Hydrogen)
            case 54: psat = ms_psat_o2(tl[layer]);  break; // O₂ (Oxygen)
            case 55: psat = ms_psat_n2(tl[layer]);  break; // N₂ (Nitrogen)
            default: psat = 1000.0; break; // Default case - arbitrary reference
        }
        // Calculate actual supersaturation (deltaP)
        double partial_pressure = vapor_vmr * pl[layer];
        double deltaP = partial_pressure - psat;
        if (deltaP < 0.0) deltaP = 0.0; // No supersaturation if at or below equilibrium
        // For numerical stability, set a small minimum if exactly zero
        if (deltaP == 0.0) deltaP = 1.0e-12;

        // Call particlesizef_local to get particle size and physics parameters
        // This function now returns all necessary physics parameters to avoid redundancy
        double r0, r1, r2, VP;
        double effective_settling_velocity, scale_height, retention;
        particlesizef_local(gravity, tl[layer], pl[layer], meanmolecular[layer],
                           species_id, Kzz, deltaP, &r0, &r1, &r2, &VP, 
                           &effective_settling_velocity, &scale_height, &retention);

        //Use r2 for sedimentation (volume-weighted radius, most relevant for fall velocity)
        // double particle_size_um = r2;
        // double particle_radius_m = particle_size_um * 1.0e-6;
        // double particle_density, accommodation_coeff, molecular_mass;
    
        // get_particle_properties(species_id, tl[layer], &particle_density, &accommodation_coeff, &molecular_mass);
        // double viscosity = 8.76e-6 * pow(tl[layer] / 293.85, 1.5);
        // double mean_free_path = 2.0 * viscosity / pl[layer] / sqrt(8.0 * meanmolecular[layer] * 1.0e-3 / PI / R_GAS / tl[layer]);
        // double knudsen = mean_free_path / particle_radius_m;
        // knudsen = fmin(knudsen, 100.0);
        // double cunningham = 1.0 + knudsen * (1.257 + 0.4 * exp(-1.1 / knudsen));
        // double fall_velocity = (2.0 * particle_radius_m * particle_radius_m * particle_density * gravity * cunningham) / (9.0 * viscosity);
        // scale_height = (KBOLTZMANN * tl[layer]) / (meanmolecular[layer] * AMU * gravity);
        // double f_sed = (fall_velocity * scale_height) / Kzz;
        // retention = 1.0 / (1.0 + f_sed);
        // if (retention < 0.001) retention = 0.001;
        // if (retention > 0.95) retention = 0.95;        

        // printf ("==============================================\n");
        // printf ("values from apply_enhanced_cloud_physics\n");
        // printf ("==============================================\n"); 
        // printf ("viscosity: %.2e\n", viscosity);
        // printf ("retention: %.2e\n", retention);
        // printf ("fall_velocity: %.2e\n", fall_velocity);
        // printf ("f_sed: %.2e\n", f_sed);
        // printf ("scale_height: %.2e\n", scale_height);
        //Store results for use in sedimentation
        if (i < MAX_CONDENSIBLES) {
            particle_radius_um[layer][i] = r2; // Use r2 (volume-weighted radius) for sedimentation
            fall_velocity_ms[layer][i] = effective_settling_velocity; // Already includes max[... - u, 0] correction
            cloud_retention[layer][i] = retention; // Already calculated with proper physics
        }
        
        // Debug output for key layer and species
        if (layer == 60 && cloud_vmr > 1.0e-15) {
            printf("Enhanced cloud physics - Layer %d, Species %d:\n", layer, species_id);
            printf("  Condensate VMR: %.2e, deltaP: %.2e Pa\n", cloud_vmr, deltaP);
            printf("  Particle sizes: r0=%.1f, r1=%.1f, r2=%.1f μm, VP=%.2e μm³\n", r0, r1, r2, VP);
            printf("  Using r2 for sedimentation: %.1f μm\n", r2);
            printf("  Effective settling velocity (v_d): %.4e m/s (includes max[... - u, 0])\n", effective_settling_velocity);
            printf("  Scale height: %.1f m\n", scale_height);
            printf("  Retention: %.3f\n", retention);
            printf("  NOTE: All physics calculated consistently in particlesizef_local\n");
        }
    }
}


// CORRECT ACKERMAN & MARLEY (2001) IMPLEMENTATION
// Applies retention factors to Graham's equilibrium amounts
// Uses real particle sizes from particlesizef_local for sedimentation physics
void apply_equilibrium_cloud_distribution(double gravity)
{
    // CORRECT APPROACH: Apply retention factors to Graham's equilibrium amounts
    // Don't redistribute - just apply sedimentation physics to what thermodynamics says should condense
    
    printf("=== APPLYING SEDIMENTATION TO GRAHAM'S EQUILIBRIUM AMOUNTS ===\n");
    printf("LAYER NUMBERING: Layer 0 = Surface (high P), Layer %d = Top (low P)\n", zbin);
    
    // Calculate total cloud amounts before sedimentation for diagnostics
    double total_cloud_before = 0.0;
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < NCONDENSIBLES; i++) {
            int species_id = CONDENSIBLES[i];
            total_cloud_before += clouds[layer][species_id];
        }
    }
    
    // Store fallen material for each species
    double fallen_material[zbin+1][NCONDENSIBLES];
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < NCONDENSIBLES; i++) {
            fallen_material[layer][i] = 0.0;
        }
    }
    
    // Apply sedimentation layer by layer (TOP TO BOTTOM)
    // In EPACRIS: layer zbin = top, layer 0 = surface
    // Material falls from high layer numbers to low layer numbers
    for (int layer = zbin; layer >= 1; layer--) {
        for (int i = 0; i < NCONDENSIBLES; i++) {
            int species_id = CONDENSIBLES[i];
            
            // Get Graham's equilibrium amount (what thermodynamics says should condense)
            double graham_condensate = clouds[layer][species_id];
            
            if (graham_condensate < 1.0e-20) continue;
            
            // Get retention factor from enhanced cloud physics
            double retention = 1.0; // Default: keep everything
            if (i < MAX_CONDENSIBLES) {
                retention = cloud_retention[layer][i];
            } else {
                // Fallback: use ALPHA_RAINOUT if enhanced physics didn't run
                retention = ALPHA_RAINOUT;
            }
            
            // Apply physical bounds to retention
            if (retention < 0.001) retention = 0.001; // Minimum 0.1% retention
            if (retention > 0.95) retention = 0.95;   // Maximum 95% retention
            
            // Calculate retained vs. falling material
            double retained = graham_condensate * retention;
            double falling = graham_condensate * (1.0 - retention);
            
            // Update this layer (retained amount only)
            clouds[layer][species_id] = retained;
            
            // CORRECTED: Add falling material to layer below (lower layer number)
            // Material falls from layer N to layer N-1 (toward surface)
            if (layer > 1) {
                fallen_material[layer-1][i] += falling;
            }
            
            // Debug output for key layer and species
            if (layer == 60 && species_id == 7 && graham_condensate > 1.0e-15) {
                printf("Layer %d, Species %d: Graham=%.2e, retention=%.3f, retained=%.2e, falling=%.2e\n",
                       layer, species_id, graham_condensate, retention, retained, falling);
                printf("  Using particle size: %.1f μm, fall velocity: %.2f m/s\n",
                       (i < MAX_CONDENSIBLES) ? particle_radius_um[layer][i] : 0.0,
                       (i < MAX_CONDENSIBLES) ? fall_velocity_ms[layer][i] : 0.0);
                printf("  Material falls from layer %d to layer %d (toward surface)\n", layer, layer-1);
            }
            
            // DEBUG: Check for very low retention (potential issue)
            if (retention < 0.1 && graham_condensate > 1.0e-15) {
                printf("WARNING: Very low retention %.3f in layer %d, species %d\n", retention, layer, species_id);
            }
        }
    }
    
    // Add fallen material to receiving layers
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < NCONDENSIBLES; i++) {
            int species_id = CONDENSIBLES[i];
            
            if (fallen_material[layer][i] > 0.0) {
                clouds[layer][species_id] += fallen_material[layer][i];
                
                // Debug significant additions
                if (layer >= 55 && layer <= 65 && species_id == 7 && fallen_material[layer][i] > 1.0e-15) {
                    printf("Layer %d receives fallen material: %.2e\n", layer, fallen_material[layer][i]);
                }
            }
        }
    }
    
    // Calculate total cloud amounts after sedimentation for diagnostics
    double total_cloud_after = 0.0;
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < NCONDENSIBLES; i++) {
            int species_id = CONDENSIBLES[i];
            total_cloud_after += clouds[layer][species_id];
        }
    }
    
    printf("Total cloud mass: before=%.3e, after=%.3e, conservation ratio=%.6f\n", 
           total_cloud_before, total_cloud_after, total_cloud_after/total_cloud_before);
    
    // DEBUG: Show cloud profile for H2O (species 7) around key layers
    printf("CLOUD PROFILE DEBUG (H2O around layers 55-65):\n");
    for (int layer = 65; layer >= 55; layer--) {
        double h2o_vmr = clouds[layer][7] / MM[layer];
        printf("  Layer %d: H2O VMR = %.2e\n", layer, h2o_vmr);
    }
    
    printf("=== SEDIMENTATION COMPLETE ===\n");
}





void ms_conv_check(double tempb[],double P[],double lapse[],int isconv[],int* ncl,int* nrl)
{
// Check each pair of layers for convective stability
//ncl, nrl, and isconv are handed over as 0 initially

    //LOCAL parameters
    int i;
    for (i=0; i<=zbin; i++) isconv[i] = 0;
    
    //TESTING!!
    //isconv[0] = 1;

    //START function
    for (i=0; i<zbin; i++)
    {
//        printf("%s %d%s %f%s%f\n","Layer",i,": Tabove=",tempb[i+1]," but adiabatically ",( tempb[i] * pow(P[i+1]/P[i],lapse[i+1])));
//        if( tempb[i+1] < ( tempb[i] * pow(P[i+1]/P[i],lapse[i+1]) /0.98 ) )
        
        // STABILITY FIX: Use larger tolerance to prevent oscillations at cloud decks
        // where latent heating can cause rapid lapse rate changes
        double stability_tolerance = 0.5;  // Increased from 0.1 K to 0.5 K
        double adiabatic_temp = tempb[i] * pow(P[i+1]/P[i], lapse[i+1]);
        
        if( tempb[i+1] < (adiabatic_temp + stability_tolerance) ) //for numerical stability
        {
        //printf("%s %d\t%s %d\t%f %s %f\n","Layer",i+1,"isconv",isconv[i],tempb[i+1],"<",( tempb[i] * pow(P[i+1]/P[i],lapse[i+1])));
            isconv[i+1] = 1;
            //isconv[i] = 1;
        }
//    printf("%s %d %s %d\n","ms_i",i,"ms_isconv",isconv[i]);
    }
    for (i=1;i<=zbin;i++)
    {
        if(isconv[i]) *ncl+=1;
        else *nrl+=1;
    }
//    printf("%s %d %s %d\n","ms_ncl",*ncl,"ms_nrl",*nrl);
}// END: void ms_conv_check()
//*********************************************************
//*********************************************************
void ms_temp_adj(double tempb[],double P[],double lapse[],int isconv[], double cp[], int ncreg[], double pot_temp[])
{
// Adjust temperatures due to convection

    //LOCAL parameters
//    int ncreg[zbin+1] = {0}; //determine convective regions
    int i,j,cr,nccount = 0;
    double potT[zbin+1]={0};//, pot_temp[zbin+1];
    double enthalpy[zbin+1]={0}, potTint[zbin+1]={0};
    double ppk[zbin+1];
    double t_old[zbin+1];

    //START function
    for (i=0;i<=zbin;i++)
    {   
        t_old[i]=tempb[i]; //copy for later
        ppk[i]=1.0;
    }

//    printf("%s\t%s\n","i","ppk[i]");
    if(isconv[0]) {nccount+=1; ncreg[0] = nccount;}
    for (i=1; i<=zbin; i++) //number convective regions
    {
//        printf("%d\t%e\n",i,ppk[i]);
        if(isconv[i]==1) 
        {
            if(isconv[i-1]==0) nccount += 1;
            ncreg[i] = nccount;
        }
    }//END: number convective regions
//    print:f("%s\t%s\t%s\n","Layer","isconv","Conv_Region");
//    for (i=zbin; i>=1; i--) printf("%d\t%d\t%d\n",i,isconv[i],ncreg[i]);

    //CALC potential Temps
//    printf("%s\t%s\t%s\n","cr","enthalpy","potT enthalpy integral");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}

    for (i=1; i<=zbin; i++)
    {
//temps too high:        enthalpy[ncreg[i]] += cp[i] * (tempb[i]+tempb[i-1])/2.0 * (P[i-1] - P[i]);
        enthalpy[ncreg[i]] += cp[i] * tempb[i] * (P[i-1] - P[i]);
//    printf("%s\t%s\t%s\n","cr","enthalpy-2","potT enthalpy integral-2");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}

        j=i;
//        printf("%s\t%s\t%s\t%s\t%s\t%s\n","i","j","ppk","Pj","Pj-1","lapse");
        while (abs(ncreg[j]-ncreg[i]) <1.0 && isconv[i]==1 && j>0)
        {
            ppk[i] *= pow(P[j] / P[j-1],lapse[j]);
//            printf("%d\t%d\t%e\t%e\t%e\t%f\n",i,j,ppk[i],P[j],P[j-1],lapse[j]);
            j -= 1;
        }
        potTint[ncreg[i]] += cp[i] * ppk[i] * (P[i-1] - P[i]);
//    printf("%s\t%s\t%s\n","cr","enthalpy-3","potT enthalpy integral-3");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}
        
//        printf("%s\t%s\t%s\n","cr","enthalpy","potT enthalpy integral");
//        for (cr=0; cr<=nccount; cr++) printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);
    }
    
//    printf("%s\t%s\t%s\t%s\n","region","enthalpy","potT enthalpy integral","potT");

    for (cr=1; cr<=nccount; cr++)
    {
        potT[cr] = enthalpy[cr] / potTint[cr];
//        printf("%d\t%e\t%e\t%e\n",cr,enthalpy[cr],potTint[cr],potT[cr]);
    }

    //ADJUST temps
    for (i=1; i<=zbin; i++) //first mid-layer:
    {
        pot_temp[i] = potT[ncreg[i]];
        if(isconv[i]) 
        {
            double new_temp = potT[ncreg[i]] * ppk[i];
            
            // STABILITY FIX: Add temperature damping to prevent oscillations at cloud decks
            // Blend old and new temperatures to smooth rapid changes
            double damping_factor = 0.3;  // Use 30% of new temperature, 70% of old
            tempb[i] = damping_factor * new_temp + (1.0 - damping_factor) * t_old[i];
            
            if (isconv[i-1]==0) {
                // At convection boundary, also apply damping to avoid discontinuities
                double new_base_temp = potT[ncreg[i]];
                tempb[i-1] = damping_factor * new_base_temp + (1.0 - damping_factor) * t_old[i-1];
            }
        }
    }
/*    for (i=1; i<zbin; i++) //now layer boundaries
    {
        tempb[i] = 0.5 * ( tl[i] + tl[i+1]);
    }*/
    //TOA:
    //tempb[zbin] = tl[zbin] + (tl[zbin] - tempb[zbin-1]);
    //ADJUST surface
    //if(isconv[1]) tempb[0] = potT[ncreg[1]] * (1 + 0.5 * (pow(P[0]/(P[1]),lapse[1]) - 1.0));
    // covered already above if(isconv[1]) tempb[0] = potT[ncreg[1]];

    //mantas debug print for conv_test
    // printf("%s\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\n",
    //     "Layer",      // Layer number (from top of atmosphere down)
    //     "cp",         // Specific heat capacity
    //     "lapse",      // Lapse rate (rate of temperature change with height)
    //     "isconv",     // Is this layer convective? (1=yes, 0=no)
    //     "Conv_Region", // Which convective region this layer belongs to
    //     "Pot_Temp",   // Potential temperature
    //     "T_old",      // Temperature before adjustment
    //     "T_new"       // Temperature after adjustment
    // );
    // for (i=zbin; i>=0; i--) {
    //     printf("%d\t%f\t%f\t%d\t%d\t\t%e\t%f\t%f\n",
    //         i,            // Layer number
    //         cp[i],        // Heat capacity for this layer
    //         lapse[i],     // Lapse rate for this layer
    //         isconv[i],    // Convective flag (1/0)
    //         ncreg[i],     // Convective region number
    //         pot_temp[i],  // Potential temperature (in scientific notation)
    //         t_old[i],     // Old temperature
    //         tempb[i]      // New temperature
    //     );

//atexit(pexit);exit(0); //ms debugging mode
}//END: void ms_temp_adj()
//*********************************************************
//*********************************************************

//=========================================================
//=== Dynamic Condensation Detection Functions ===========

void initialize_condensibles_mode() {
    // Initialize condensibles based on CONDENSATION_MODE
    
    if (CONDENSATION_MODE == 0) {
        // Manual mode: use predefined lists from AlphaAb.h
        NCONDENSIBLES = NCONDENSIBLES_MANUAL;
        
        // Directly copy values from compound literal macros
        // This avoids the "array initialized from non-constant array expression" error
        for (int i = 0; i < NCONDENSIBLES; i++) {
            CONDENSIBLES[i] = CONDENSIBLES_MANUAL[i];
        }
        
        printf("CONDENSATION MODE: Manual - Using predefined %d condensible species\n", NCONDENSIBLES);
        printf("Manual condensibles: ");
        for (int i = 0; i < NCONDENSIBLES; i++) {
            printf("%d ", CONDENSIBLES[i]);
        }
        printf("\n");
        printf("Alpha rainout value: %.3f (%.1f%% retention)\n", ALPHA_RAINOUT, ALPHA_RAINOUT * 100.0);
        
    } else if (CONDENSATION_MODE == 1) {
        // Automatic mode: detect condensibles dynamically
        printf("CONDENSATION MODE: Automatic - Will detect condensibles dynamically\n");
        printf("Alpha rainout value: %.3f (%.1f%% retention)\n", ALPHA_RAINOUT, ALPHA_RAINOUT * 100.0);
        NCONDENSIBLES = 0; // Will be set by detect_condensibles_atmosphere()
        
    } else if (CONDENSATION_MODE == 2) {
        // Hybrid mode: start with manual list, validate and expand automatically
        NCONDENSIBLES = NCONDENSIBLES_MANUAL;
        
        // Directly copy values from compound literal macros
        // This avoids the "array initialized from non-constant array expression" error
        for (int i = 0; i < NCONDENSIBLES; i++) {
            CONDENSIBLES[i] = CONDENSIBLES_MANUAL[i];
        }
        
        printf("CONDENSATION MODE: Hybrid - Starting with manual list, will validate and expand\n");
        printf("Alpha rainout value: %.3f (%.1f%% retention)\n", ALPHA_RAINOUT, ALPHA_RAINOUT * 100.0);
    }
}

int check_species_condensible(int species_id, double temp, double pressure, double partial_pressure) {
    // Check if a species should be considered condensible
    double psat = 0.0;
    
    // Get saturation pressure for the species
    switch(species_id) {
        case 7:  // H2O
            psat = ms_psat_h2o(temp);
            break;
        case 9:  // NH3
            psat = ms_psat_nh3(temp);
            break;
        case 20: // CO
            psat = ms_psat_co(temp);
            break;
        case 21: // CH4
            psat = ms_psat_ch4(temp);
            break;
        case 52: // CO2
            psat = ms_psat_co2(temp);
            break;
        case 53: // H2
            psat = ms_psat_h2(temp);
            break;
        case 54: // O2
            psat = ms_psat_o2(temp);
            break;
        case 55: // N2
            psat = ms_psat_n2(temp);
            break;
        case 45: // H2S
            psat = ms_psat_h2s(temp);
            break;
        default:
            return 0; // Unknown species, treat as dry
    }
    
    // Check if saturation pressure is reasonable (not the 1e+20 "above critical" case)
    if (psat >= 1e+15) {
        return 0; // Above critical point or undefined - definitely dry
    }
    
    // Calculate saturation ratio
    if (psat <= 0.0) {
        return 0; // Invalid saturation pressure
    }
    
    double saturation_ratio = partial_pressure / psat;
    
    // If we're within SATURATION_THRESHOLD of saturation, treat as condensible
    if (saturation_ratio > SATURATION_THRESHOLD) {
        return 1; // Condensible (Xv)
    }
    
    return 0; // Dry (Xd)
}

void update_condensibles_list(int layer) {
    // Update condensibles list for a specific layer (for layer-specific detection)
    // This function can be called during each convective adjustment if needed
    
    if (CONDENSATION_MODE == 0) {
        return; // Manual mode - no updates needed
    }
    
    // For automatic or hybrid modes, check all species in this layer
    int temp_condensibles[MAX_CONDENSIBLES];
    int temp_ncondensibles = 0;
    
    // Check all species that have saturation pressure functions
    int species_to_check[] = {7, 9, 20, 21, 52, 53, 54, 55, 45}; // H2O, NH3, CO, CH4, CO2, H2, O2, N2, H2S
    int num_species_to_check = 9;
    
    for (int i = 0; i < num_species_to_check; i++) {
        int species_id = species_to_check[i];
        
        // Calculate partial pressure for this species
        double partial_pressure = (xx[layer][species_id] / MM[layer]) * pl[layer];
        
        if (check_species_condensible(species_id, tl[layer], pl[layer], partial_pressure)) {
            // Check if species is already in the list
            int already_in_list = 0;
            for (int j = 0; j < temp_ncondensibles; j++) {
                if (temp_condensibles[j] == species_id) {
                    already_in_list = 1;
                    break;
                }
            }
            
            if (!already_in_list && temp_ncondensibles < MAX_CONDENSIBLES) {
                temp_condensibles[temp_ncondensibles] = species_id;
                temp_ncondensibles++;
            }
        }
    }
    
    // For layer-specific updates, you could store this information
    // For now, we'll use the global detection approach
}

void detect_condensibles_atmosphere() {
    // Detect condensibles across the entire atmosphere
    
    if (CONDENSATION_MODE == 0) {
        return; // Manual mode - no detection needed
    }
    
    printf("DETECTING CONDENSIBLE SPECIES ACROSS ATMOSPHERE...\n");
    
    int detected_condensibles[MAX_CONDENSIBLES];
    int detected_count = 0;
    
    // Species that have saturation pressure functions
    int species_to_check[] = {7, 9, 20, 21, 52, 53, 54, 55, 45}; // H2O, NH3, CO, CH4, CO2, H2, O2, N2, H2S
    int num_species_to_check = 9;
    
    for (int i = 0; i < num_species_to_check; i++) {
        int species_id = species_to_check[i];
        int is_condensible_anywhere = 0;
        
        // Check across all atmospheric layers
        for (int layer = 1; layer <= zbin; layer++) {
            double partial_pressure = (xx[layer][species_id] / MM[layer]) * pl[layer];
            
            if (check_species_condensible(species_id, tl[layer], pl[layer], partial_pressure)) {
                is_condensible_anywhere = 1;
                break; // Found at least one layer where it's condensible
            }
        }
        
        if (is_condensible_anywhere && detected_count < MAX_CONDENSIBLES) {
            detected_condensibles[detected_count] = species_id;
            detected_count++;
            printf("  Species %d detected as condensible\n", species_id);
        }
    }
    
    if (CONDENSATION_MODE == 1) {
        // Pure automatic mode: use only detected species
        NCONDENSIBLES = detected_count;
        for (int i = 0; i < NCONDENSIBLES; i++) {
            CONDENSIBLES[i] = detected_condensibles[i];
        }
        
    } else if (CONDENSATION_MODE == 2) {
        // Hybrid mode: merge manual list with detected species
        int original_count = NCONDENSIBLES;
        
        // Add any newly detected species not in the manual list
        for (int i = 0; i < detected_count; i++) {
            int species_id = detected_condensibles[i];
            int already_in_manual = 0;
            
            for (int j = 0; j < original_count; j++) {
                if (CONDENSIBLES[j] == species_id) {
                    already_in_manual = 1;
                    break;
                }
            }
            
            if (!already_in_manual && NCONDENSIBLES < MAX_CONDENSIBLES) {
                CONDENSIBLES[NCONDENSIBLES] = species_id;
                NCONDENSIBLES++;
                printf("  Added species %d to manual list (detected automatically)\n", species_id);
            }
        }
    }
    
    printf("FINAL CONDENSIBLES LIST: %d species\n", NCONDENSIBLES);
    printf("Species IDs: ");
    for (int i = 0; i < NCONDENSIBLES; i++) {
        printf("%d ", CONDENSIBLES[i]);
    }
    printf("\n");
    printf("Using single alpha value: %.3f (%.1f%% retention) for all species\n", ALPHA_RAINOUT, ALPHA_RAINOUT * 100.0);
}

void report_condensibles_changes(int iteration) {
    // Report changes in condensibles list
    printf("Condensibles changes after iteration %d:\n", iteration);
    printf("Species IDs: ");
    for (int i = 0; i < NCONDENSIBLES; i++) {
        printf("%d ", CONDENSIBLES[i]);
    }
    printf("\n");
    printf("Using single alpha value: %.3f (%.1f%% retention) for all species\n", ALPHA_RAINOUT, ALPHA_RAINOUT * 100.0);
}

//=========================================================

double ms_psat_h2o(double temp) //calculate saturation pressure of water
{
    double a,P;
    
    if (temp<273.16) 
    {
        // over ice
        // Formulation from Murphy & Koop (2005)
        P = exp(9.550426-5723.265/temp+3.53068*log(temp)-0.00728332*temp); // [Pa] 
    }
    else
    {
        // over water
        // Formulation from Seinfeld & Pandis (2006)
        a = 1-373.15/temp;
        P = 101325*exp(13.3185*a-1.97*a*a-0.6445*a*a*a-0.1229*a*a*a*a); // [Pa]
    }
    
    return P;
}

//*********************************************************
//*********************************************************
double ms_psat_nh3(double temp) //calculate saturation pressure of ammonia
{
    double P;
    
    if (temp<195.40) 
    {
        // solid
        // Formulation from Lodders, Fegley (1998)
        P = pow(10.0,6.9-1588/temp) * 1.e+5; // [Pa] 
    }
    else if(temp<300.00) 
    {
        // liquid
	P = pow(10.0,5.201-1248/temp) * 1.e+5; // [Pa] 
    }
    else
    {
	// undefined - set for unsaturated case
	P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_co(double temp) //calculate saturation pressure of carbon monoxide
{
    double P;
    
    if (temp < 68.12) // below triple point
    {
        // solid CO
        // From Fray & Schmitt (2009)
        P = exp(26.97 - 764.2/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 134.45) // between triple point and critical point
    {
        // liquid CO
        // From Span & Wagner (1996)
        P = exp(24.63 - 697.4/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_ch4(double temp) //calculate saturation pressure of methane
{
    double P;
    
    if (temp < 90.67) // below triple point
    {
        // solid methane
        // From Fray & Schmitt (2009)
        P = exp(27.48 - 1190.0/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 190.44) // between triple point and critical point
    {
        // liquid methane
        // From Goodwin (1974)
        P = exp(27.1688 - 1022.9/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_co2(double temp) //calculate saturation pressure of carbon dioxide
{
    double P;
    
    if (temp < 216.54) // below triple point
    {
        // solid CO2 (dry ice)
        // Based on Fray & Schmitt (2009)
        P = exp(27.85 - 3182.4/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 304.2) // between triple point and critical point
    {
        // liquid CO2
        // From Span & Wagner (1996)
        P = exp(35.34 - 2648.0/temp - 2.74*log(temp)) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_h2(double temp) //calculate saturation pressure of hydrogen
{
    double P;
    
    if (temp < 13.95) // below triple point
    {
        // solid hydrogen
        // From Roder et al. (1973)
        P = exp(24.215 - 101.1/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 33.2) // between triple point and critical point
    {
        // liquid hydrogen
        // From Leachman et al. (2009)
        P = exp(22.5981 - 91.55/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_o2(double temp) //calculate saturation pressure of oxygen
{
    double P;
    
    if (temp < 54.3) // below triple point
    {
        // solid oxygen
        // Based on Fray & Schmitt (2009)
        P = exp(25.75 - 734.6/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 154.54) // between triple point and critical point
    {
        // liquid oxygen
        // From Jacobsen et al. (1997)
        P = exp(24.632 - 732.3/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_n2(double temp) //calculate saturation pressure of nitrogen
{
    double P;
    
    if (temp < 63.14) // below triple point
    {
        // solid nitrogen
        // Based on Fray & Schmitt (2009)
        P = exp(27.517 - 772.29/temp) * 1.0e+0; // [Pa]
    }
    else if (temp < 126.2) // between triple point and critical point
    {
        // liquid nitrogen
        // From Span et al. (2000)
        P = exp(24.31 - 704.55/temp) * 1.0e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_h2s(double temp) //calculate saturation pressure of hydrogen sulfide
{
    double P;
    
    if (temp < 187.7) // below triple point (187.7 K)
    {
        // solid H2S
        // NIST Antoine equation for solid H2S (138.8-212.8 K)
        // log10(P_bar) = A - B/(T + C)
        // A = 4.43681, B = 829.439, C = -25.412
        double log10_P_bar = 4.43681 - 829.439/(temp - 25.412);
        P = pow(10.0, log10_P_bar) * 1.0e+5; // Convert bar to Pa
    }
    else if (temp < 373.1) // between triple point and critical point (373.1 K)
    {
        // liquid H2S
        // NIST Antoine equation for liquid H2S (212.8-349.5 K)
        // log10(P_bar) = A - B/(T + C)
        // A = 4.52887, B = 958.587, C = -0.539
        double log10_P_bar = 4.52887 - 958.587/(temp - 0.539);
        P = pow(10.0, log10_P_bar) * 1.0e+5; // Convert bar to Pa
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_latent(int mol, double temp)
{
    double l;
    double t_crit,t_triple;  //Temps [K] of boiling point, triple point, and critical point (not used so far)
    double lvap_bp, lvap_tp; //latent heats [J/kg] for liquid-vapor phase change at boiling point and triple point
    double lfuse;            //latent heat [J/kg] for liquid-solid phase change
    double lsubl;            // latent heat [J/kg] for solid-vapor phase change
    double molweight;        // g/mol of each molecule to convert to J/mol
    
// Values from Pierrehumbert 2010 2.3.2, Table2.1, unless stated otherwise:
// NIST lvap and lsubl calculated from dH[kj/mol] / MolWeight [g/mol] *1e6
//H2O
    if(mol==7)  {t_crit= 647.1; t_triple=273.16; 
        lvap_bp=22.55e+05; lvap_tp=24.93e+05; 
        lfuse=3.34e+05; lsubl=28.4e+05;
        molweight = 18.02;
    }
//NH3
    if(mol==9)  {t_crit= 405.5; t_triple=195.4; 
        lvap_bp=13.71e+05; lvap_tp=16.58e+05; 
        lfuse=3.314e+05; lsubl=19.89e+05;
        molweight = 17.03;
    } 
//CO --- Temps: NIST Chemistry Webbook SRD 69
    if(mol==20) {t_crit= 134.45; t_triple=68.12; 
        lvap_bp=2.14e+05; lvap_tp=2.14e+05 /*NIST dH=6.0 for both*/; 
        lfuse=2.98e+04 /*engineeringtoolbox.com - calculated from BTU/lb*/; 
        lsubl=2.89e+05 /*NIST dH=8.1*/;
        molweight = 28.01;
            //since lsubl=lfuse+lvap_bp+cp_liquid*(Tboil=81.65 - Tmelt=68.12) -> cp_CO_liquid =~ 4006 [J/kg/K] = 0.1122 [kJ/mol/K] ?
    } 
//CH4
    if(mol==21) {t_crit= 190.44; t_triple=90.67; 
        lvap_bp=5.1e+05; lvap_tp=5.36e+05; 
        lfuse=5.868e+04; lsubl=5.95e+05;
        molweight = 16.04;
    } 
//CO2
    if(mol==52) {t_crit= 304.2; t_triple=216.54; 
        lvap_bp=3.79e+05 /*NIST from dH=16.7*/; lvap_tp=3.97e+05; 
        lfuse=1.96e+05; lsubl=5.93e+05;
        molweight = 44.01;
    } 
//H2
    if(mol==53) {t_crit= 33.2; t_triple=13.95; 
        lvap_bp=4.54e+05; lvap_tp=4.54e+05 /*from bp since no values avail.*/; 
        lfuse=5.82e+04; 
        lsubl=5.65e+05; //=lfuse+lvap_bp +dH_liquid [=cp_liquid*(Tboil=20.39K - Ttriple)] - since no data avail.
        molweight = 2.02;
            //cp[tp]=~6.57e3[J/kg/K], cp[bp]=~9.77e3[J/Kg/K] -> dH_liquid=0.5*(6.57e3+9.77e3)*(20.39-13.95) = 5.26e4
            //lsubl=lfuse+dH_liquid+lvap_bp
            //Liquid H2 data from KAERI/TR-2723/2004 (https://inis.iaea.org/collection/NCLCollectionStore/_Public/36/045/36045728.pdf)
    } 
//O2
    if(mol==54) {t_crit= 154.54; t_triple=54.3; 
        lvap_bp=2.13e+05; lvap_tp=2.42e+05; 
        lfuse=13.9e+04; lsubl=2.56e+05;
        molweight = 32.00;
    } 
//N2
    if(mol==55) {t_crit= 126.2; t_triple=63.14; 
        lvap_bp=1.98e+05; lvap_tp=2.18e+05; 
        lfuse=2.573e+04; lsubl=2.437e+05;
        molweight = 28.01;
    } 
//H2S --- NIST Chemistry Webbook SRD 69
    if(mol==45) {t_crit= 373.3; t_triple=187.7; 
        lvap_bp=18.6e+03 * 1.0e+3 / 34.081; // NIST 18.6 kJ/mol at 243K -> J/kg
        lvap_tp=21.9e+03 * 1.0e+3 / 34.081; // NIST 21.9 kJ/mol at 200K -> J/kg
        lfuse=(25.4e+03 - 21.9e+03) * 1.0e+3 / 34.081; // lfuse = lsubl - lvap_tp
        lsubl=25.4e+03 * 1.0e+3 / 34.081; // NIST 25.4 kJ/mol at 175K -> J/kg
        molweight = 34.081;
    } 

//Now calculate the appropriate l:
    if (temp<=t_triple) {l = lsubl;}
    else if (temp<t_crit) {l = (t_crit-temp)/(t_crit-t_triple)*lvap_tp + (temp-t_triple)/(t_crit-t_triple)*lvap_bp;}
    else if (temp>=t_crit) {l = lvap_bp;}

    l = l * molweight / 1.0e+3; // convert J/kg to J/mol
    //printf("%s%f\n","Psat = ",*P);
    return l;
}// END: couble ms_latent()
//*********************************************************
//*********************************************************


