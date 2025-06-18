#include <math.h>
#include "constant.h"
#include "condensed_heat.h"  // Include new condensed heat capacity functions

// Note: Global variables for dynamic condensibles management are defined in the main file
// NCONDENSIBLES and CONDENSIBLES[] are declared in global_temp.h
// ALPHA_RAINOUT is a single constant value defined in AlphaAb.h

//=========================================================
//=== Function identifiers ================================

void ms_adiabat(int lay, double lapse[], double xxHe, double* cp, double saturation_ratios[]); // calculate lapse rate (NO rainout)
void ms_rainout(int lay, double* mass_loss_ratio); // apply rainout physics
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
    double Xd,cp_d; //mole fraction of non-condensible gas, heat capacity of non-condensible gas
    double Xv[NCONDENSIBLES],Xvold,cp_v[NCONDENSIBLES]; //condensibles mol fractions in vapor form and heat cap. [Xvold is the old value of Xv]
    double Xc[NCONDENSIBLES],cp_c[NCONDENSIBLES]; //condensibles fractions in condensed form (Xc) and heat cap. (cp_c)
    double alpha[NCONDENSIBLES], beta[NCONDENSIBLES], latent[NCONDENSIBLES]; // Latent heat and B=L/RT for each condensible

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


    // STEP 1: Calculate total heat capacity of the ENTIRE atmosphere (all species)
    // cpxx = Σ(number_density × molar_heat_capacity) for ALL species  [molecules/cm³] × [J/(mol·K)]
    cpxx = xxHe*HeHeat(tl[lay])+xx[lay][7]*H2OHeat(tl[lay])+xx[lay][9]*NH3Heat(tl[lay])+
        xx[lay][20]*COHeat(tl[lay])+xx[lay][21]*CH4Heat(tl[lay])+xx[lay][45]*H2SHeat(tl[lay])+
        xx[lay][52]*CO2Heat(tl[lay])+xx[lay][53]*H2Heat(tl[lay])+xx[lay][54]*O2Heat(tl[lay])+xx[lay][55]*N2Heat(tl[lay]);
        // potentially add more heat capacities of possible atmospheric species here

    // STEP 2: Calculate total number density of ALL species 
    MM_cp = xxHe+xx[lay][7]+xx[lay][9]+xx[lay][20]+xx[lay][21]+xx[lay][45]+xx[lay][52]+xx[lay][53]+xx[lay][54]+xx[lay][55];
    
    // STEP 3: Calculate average molar heat capacity of ENTIRE atmosphere
    cp_d = cpxx / MM_cp;
    
    // if(lay==1) printf("STEP 1-3: cpxx=%.3e, MM_cp=%.3e, cp_d_initial=%.3e\n", cpxx, MM_cp, cp_d);
    
    // STEP 4: Loop through each condensible species to subtract them out from total heat capacity and number density
    for (i=0; i<NCONDENSIBLES; i++)
    {
        // if(lay==1)  printf("STEP 4.%d: Processing condensible species %d\n", i+1, CONDENSIBLES[i]);
        // if(lay==1)  printf("  Before: cpxx=%.3e, MM_cp=%.3e\n", cpxx, MM_cp);
        
        if(CONDENSIBLES[i]==7)  {
            cp_v[i] = H2OHeat(tl[lay]); //heat capacity of water vapor
            cp_c[i] = H2O_liquid_heat_capacity(tl[lay]); //heat capacity of condensed water - NOW USING NEW LIBRARY
            //if(lay==1) printf("H2O heat capacities - Vapor: %.3e, Liquid: %.3e\n", cp_v[i], cp_c[i]);
            MM_cp-=xx[lay][7];  // SUBTRACT H2O number density from total number density
            cpxx-=xx[lay][7]*H2OHeat(tl[lay]); // SUBTRACT H2O heat capacity contribution from total heat capacity
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        if(CONDENSIBLES[i]==9)  {
            cp_v[i] = NH3Heat(tl[lay]); 
            cp_c[i] = NH3_liquid_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][9];  // SUBTRACT NH3 number density  
            cpxx-=xx[lay][9]*NH3Heat(tl[lay]); // SUBTRACT NH3 heat capacity contribution
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        if(CONDENSIBLES[i]==20) {
            cp_v[i] = COHeat(tl[lay]);  
            cp_c[i] = CO_condensed_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][20]; 
            cpxx-=xx[lay][20]*COHeat(tl[lay]); 
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        
        if(CONDENSIBLES[i]==21) {
            cp_v[i] = CH4Heat(tl[lay]); 
            cp_c[i] = CH4_condensed_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][21]; 
            cpxx-=xx[lay][21]*CH4Heat(tl[lay]); 
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        
        if(CONDENSIBLES[i]==52) {
            cp_v[i] = CO2Heat(tl[lay]); 
            cp_c[i] = CO2_solid_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][52]; 
            cpxx-=xx[lay][52]*CO2Heat(tl[lay]); 
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        
        if(CONDENSIBLES[i]==53) {
            cp_v[i] = H2Heat(tl[lay]);  
            cp_c[i] = H2_condensed_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][53]; 
            cpxx-=xx[lay][53]*H2Heat(tl[lay]); 
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        
        if(CONDENSIBLES[i]==54) {
            cp_v[i] = O2Heat(tl[lay]);  
            cp_c[i] = O2_condensed_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][54]; 
            cpxx-=xx[lay][54]*O2Heat(tl[lay]); 
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        
        if(CONDENSIBLES[i]==55) {
            cp_v[i] = N2Heat(tl[lay]);  
            cp_c[i] = N2_condensed_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][55];  
            cpxx-=xx[lay][55]*N2Heat(tl[lay]); 
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
        
        if(CONDENSIBLES[i]==45) {
            cp_v[i] = H2SHeat(tl[lay]);  
            cp_c[i] = H2S_condensed_heat_capacity(tl[lay]); // NOW USING NEW LIBRARY
            MM_cp-=xx[lay][45];  
            cpxx-=xx[lay][45]*H2SHeat(tl[lay]); 
            alpha[i] = ALPHA_RAINOUT; // Use config-defined alpha value
        }
    
        // Calculate latent heat of the condensible species
        latent[i] = ms_latent(CONDENSIBLES[i],tl[lay]);
        beta[i] =  latent[i] / R_GAS / tl[lay];
        

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
    
    // After subtracting the condensable species, calculate cp_d ONCE:
    cp_d = cpxx / MM_cp;  // Now cp_d = average heat capacity of remaining (dry) species

    //if(lay==1)  printf("FINAL: cpxx=%.3e, MM_cp=%.3e, cp_d_final=%.3e\n", cpxx, MM_cp, cp_d);

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
    
}// END: void ms_rainout()
//*********************************************************
//*********************************************************
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

//atexit(pexit);exit(0); //ms debugging mode



