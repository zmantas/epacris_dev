'''
Function for generating a 1D temperature structure using the Parmentier & Guillot (2014/15)
analytical non-grey picket fence radiative-equilibrium approach. 
A convective correction to the radiative-equilibrium solution can also be optionally applied.

NOTES: 
'''

import numpy as np
from atm_module import k_Ross_Freedman

def Bond_Parmentier(Teff0, grav):

  # Calculates the Bond Albedo according to Parmentier et al. (2015) expression

  # Input:
  # Teff0 - Atmospheric profile effective temperature [K] with zero albedo
  # grav - Surface gravity of planet [cm s-2]

  # Output:
  # AB - Bond albedo

  grav = grav/100.0 # Convert to m s-2

  # a and b cofficents dependent on T_eff and grav
  if (Teff0 <= 250.0):
    a = -0.335 * grav**(0.070)
    b = 0.0
  elif ((Teff0 > 250.0) and (Teff0 <= 750.0)):
    a = -0.335 * grav**(0.070) + 2.149 * grav**(0.135)
    b = -0.896 * grav**(0.135)
  elif ((Teff0 > 750.0) and (Teff0 < 1250.0)):
    a = -0.335 * grav**(0.070) -  0.428 * grav**(0.135)
    b = 0.0
  elif (Teff0 >= 1250.0):
    a = 16.947 - 3.174 * grav**(0.070) - 4.051 * grav**(0.135)
    b = -5.472 + 0.917 * grav**(0.070) + 1.170 * grav**(0.135)

  # Final Bond Albedo expression
  AB = 10.0**(a + b * np.log10(Teff0))

  return AB

def gam_Parmentier(Teff, table_num):

  # Calculates 3 band grey visual gamma values and 2 picket fence IR gamma values
  # according to the coefficients and equations in:
  # Parmentier & Menou (2014) and Parmentier et al. (2015)

  # Input:
  # Teff - Effective temperature [K] (See Parmentier papers for various ways to calculate this)
  # for non-irradiated atmosphere Teff = Tint
  # table_num - Table selection from Parmentier et al. (2015): 1 = w. TiO/VO, 2 = w.o. TiO/VO

  # Output:
  # gam_V(3) - gamma ratio for 3 visual bands (gam_V = kV_Ross/kIR_Ross)
  # beta_V(3) - fraction of total incident stellar flux in band (1/3 for Parmentier values)
  # Beta - equivalent bandwidth for picket fence IR model
  # gam_1 - gamma ratio for IR band 1 (gam_1 = kIR_1/kIR_Ross)
  # gam_2 - gamma ratio for IR band 2 (gam_2 = kIR_2/kIR_Ross)
  # gam_P - gamma ratio for Planck mean (gam_P = kIR_Planck/kIR_Ross)
  # tau_lim - tau limit variable (usually for IC system)

  # Log 10 T_eff variables
  l10T = np.log10(Teff)
  l10T2 = l10T**2

  if (table_num == 1):
    # First table in Parmentier et al. (2015) w. TiO/VO
    # Start large if statements with visual band and Beta coefficients
    if (Teff <= 200.0):
      aV1 = -5.51 ; bV1 = 2.48
      aV2 = -7.37 ; bV2 = 2.53
      aV3 = -3.03 ; bV3 = -0.20
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 200.0) and (Teff <= 300.0)):
      aV1 = 1.23 ; bV1 = -0.45
      aV2 = 13.99 ; bV2 = -6.75
      aV3 = -13.87 ; bV3 = 4.51
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 300.0) and (Teff <= 600.0)):
      aV1 = 8.65 ; bV1 = -3.45
      aV2 = -15.18 ; bV2 = 5.02
      aV3 = -11.95 ; bV3 = 3.74
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 600.0) and (Teff <= 1400.0)):
      aV1 = -12.96 ; bV1 = 4.33
      aV2 = -10.41 ; bV2 = 3.31
      aV3 = -6.97 ; bV3 = 1.94
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 1400.0) and (Teff < 2000.0)):
      aV1 = -23.75 ; bV1 = 7.76
      aV2 = -19.95 ; bV2 = 6.34
      aV3 = -3.65 ; bV3 = 0.89
      aB = 0.84  ; bB = 0.0
    elif (Teff >= 2000.0):
      aV1 = 12.65 ; bV1 = -3.27
      aV2 = 13.56 ; bV2 = -3.81
      aV3 = -6.02 ; bV3 = 1.61
      aB = 6.21  ; bB = -1.63

    # gam_P coefficients
    aP = -2.36
    bP = 13.92
    cP = -19.38

  elif (table_num == 2):
    # Appendix table from Parmentier et al. (2015) - without TiO and VO
    if (Teff <= 200.0):
      aV1 = -5.51 ; bV1 = 2.48
      aV2 = -7.37 ; bV2 = 2.53
      aV3 = -3.03 ; bV3 = -0.20
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 200.0) and (Teff <= 300.0)):
      aV1 = 1.23 ; bV1 = -0.45
      aV2 = 13.99 ; bV2 = -6.75
      aV3 = -13.87 ; bV3 = 4.51
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 300.0) and (Teff <= 600.0)):
      aV1 = 8.65 ; bV1 = -3.45
      aV2 = -15.18 ; bV2 = 5.02
      aV3 = -11.95 ; bV3 = 3.74
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 600.0) and (Teff <= 1400.0)):
      aV1 = -12.96 ; bV1 = 4.33
      aV2 = -10.41 ; bV2 = 3.31
      aV3 = -6.97 ; bV3 = 1.94
      aB = 0.84  ; bB = 0.0
    elif ((Teff > 1400.0) and (Teff < 2000.0)):
      aV1 = -1.68 ; bV1 = 0.75
      aV2 = 6.96 ; bV2 = -2.21
      aV3 = 0.02 ; bV3 = -0.28
      aB = 3.0  ; bB = -0.69
    elif (Teff >= 2000.0):
      aV1 = 10.37 ; bV1 = -2.91
      aV2 = -2.4 ; bV2 = 0.62
      aV3 = -16.54 ; bV3 = 4.74
      aB = 3.0  ; bB = -0.69

    # gam_P coefficients
    if (Teff <= 1400.0):
      aP = -2.36
      bP = 13.92
      cP = -19.38
    else:
      aP = -12.45
      bP = 82.25
      cP = -134.42

  # Calculation of all values

  gam_V = np.zeros(3)
  Beta_V = np.zeros(3)
  Beta = np.zeros(2)

  # Visual band gamma
  gam_V[0] = 10.0**(aV1 + bV1 * l10T)
  gam_V[1] = 10.0**(aV2 + bV2 * l10T)
  gam_V[2] = 10.0**(aV3 + bV3 * l10T)

  #Visual band fractions
  Beta_V[:] = 1.0/3.0

  # gamma_Planck - if < 1 then make it grey approximation (k_Planck = k_Ross, gam_P = 1)
  gam_P = 10.0**(aP * l10T2 + bP * l10T + cP)
  if (gam_P < 1.0000001):
    gam_P = 1.0000001

  # equivalent bandwidth value
  Beta[0] = aB + bB * l10T
  Beta[1] = 1.0 - Beta[0]

  # IR band kappa1/kappa2 ratio - Eq. 96 from Parmentier & Menou (2014)
  RT = (gam_P - 1.0)/(2.0*Beta[0]*Beta[1])
  R = 1.0 + RT + np.sqrt(RT**2 + RT)

  # gam_1 and gam_2 values - Eq. 92, 93 from Parmentier & Menou (2014)
  gam_1 = Beta[0] + R - Beta[0]*R
  gam_2 = gam_1 / R

  # Calculate tau_lim parameter
  tau_lim = 1.0/(gam_1*gam_2) * np.sqrt(gam_P/3.0)

  return gam_V, Beta_V, Beta, gam_1, gam_2, gam_P, tau_lim

def Parmentier_T_p(nlay, pl, Tint, mu_in, Tirr, grav, met, table_num):

  mu = np.maximum(1e-6,mu_in)
    
  # Effective temperature parameter
  Tmu = (mu * Tirr**4)**(1.0/4.0)

  # Find Bond albedo of planet - Bond albedo is given by mu = 1/sqrt(3)
  Teff0 = (Tint**4 + (1.0/np.sqrt(3.0)) * Tirr**4)**(1.0/4.0)

  Bond = Bond_Parmentier(Teff0, grav)

  Teff = (Tint**4 + (1.0 - Bond) * mu * Tirr**4)**(1.0/4.0)

  a2 = np.zeros(3)
  a3 = np.zeros(3)
  b1 = np.zeros(3)
  b2 = np.zeros(3)
  b3 = np.zeros(3)
  Av1 = np.zeros(3)
  Av2 = np.zeros(3)
  C = np.zeros(3)
  D = np.zeros(3)
  E = np.zeros(3)

  gam_V = np.zeros(3)
  Beta_V = np.zeros(3)
  Beta = np.zeros(2)

  tau = np.zeros(nlay + 1)
  kRoss = np.zeros(nlay)
  Tl = np.zeros(nlay)

  # Find the V band gamma, beta and IR gamma and beta ratios for this profile
  # Passed mu, so make lat = acos(mu) and lon = 0
  gam_V, Beta_V, Beta, gam_1, gam_2, gam_P, tau_lim = gam_Parmentier(Teff, table_num)

  gam_V[:] = gam_V[:] / mu

  # Hard work starts here - first calculate all the required coefficients
  At1 = gam_1**2*np.log(1.0 + 1.0/(tau_lim*gam_1))
  At2 = gam_2**2*np.log(1.0 + 1.0/(tau_lim*gam_2))
  Av1[:] = gam_1**2*np.log(1.0 + gam_V[:]/gam_1)
  Av2[:] = gam_2**2*np.log(1.0 + gam_V[:]/gam_2)

  a0 = 1.0/gam_1 + 1.0/gam_2

  a1 = -1.0/(3.0 * tau_lim**2) * (gam_P/(1.0-gam_P) * (gam_1 + gam_2 - 2.0)/(gam_1 + gam_2) \
    + (gam_1 + gam_2)*tau_lim - (At1 + At2)*tau_lim**2)

  a2[:] = tau_lim**2/(gam_P*gam_V[:]**2) \
    * ((3.0*gam_1**2-gam_V[:]**2)*(3.0*gam_2**2-gam_V[:]**2)*(gam_1+gam_2) \
    - 3.0*gam_V[:]*(6.0*gam_1**2*gam_2**2-gam_V[:]**2*(gam_1**2+gam_2**2))) \
    / (1.0-gam_V[:]**2 * tau_lim**2)

  a3[:] = -tau_lim**2*(3.0*gam_1**2-gam_V[:]**2)*(3.0*gam_2**2-gam_V[:]**2)*(Av2[:]+Av1[:]) \
    / (gam_P*gam_V[:]**3*(1.0-gam_V[:]**2*tau_lim**2))

  b0 = 1.0/(gam_1*gam_2/(gam_1-gam_2)*(At1-At2)/3.0-(gam_1*gam_2)**2/np.sqrt(3.0*gam_P)-(gam_1*gam_2)**3 \
    / ((1.0-gam_1)*(1.0-gam_2)*(gam_1+gam_2)))

  b1[:] = gam_1*gam_2*(3.0*gam_1**2-gam_V[:]**2)*(3.0*gam_2**2-gam_V[:]**2)*tau_lim**2 \
    / (gam_P*gam_V[:]**2*(gam_V[:]**2*tau_lim**2-1.0))

  b2[:] = 3.0*(gam_1+gam_2)*gam_V[:]**3/((3.0*gam_1**2-gam_V[:]**2)*(3.0*gam_2**2-gam_V[:]**2))

  b3[:] = (Av2[:]-Av1[:])/(gam_V[:]*(gam_1-gam_2))

  A = 1.0/3.0*(a0+a1*b0)
  B = -1.0/3.0*(gam_1*gam_2)**2/gam_P*b0
  C[:] = -1.0/3.0*(b0*b1[:]*(1.0+b2[:]+b3[:])*a1+a2[:]+a3[:])
  D[:] = 1.0/3.0*(gam_1*gam_2)**2/gam_P*b0*b1[:]*(1.0+b2[:]+b3[:])
  E[:] = (3.0-(gam_V[:]/gam_1)**2)*(3.0-(gam_V[:]/gam_2)**2)/(9.0*gam_V[:]*((gam_V[:]*tau_lim)**2-1.0))

  # T-p structure calculation - we follow exactly V. Parmentier's method
  # Estimate the skin temperature by setting tau = 0
  tau[0] = 0.0
  Tskin = 3.0*Tint**4/4.0*(tau[0]+A+B*np.exp(-tau[0]/tau_lim)) + sum(3.0*Beta_V[:] \
    * Tmu**4/4.0*(C[:]+D[:]*np.exp(-tau[0]/tau_lim)+E[:]*np.exp(-gam_V[:]*tau[0])))
  Tskin = Tskin**(1.0/4.0)

  # Estimate the opacity TOA at the skin temperature - assume this is = first layer opacity
  kRoss[0] =  k_Ross_Freedman(Tskin, pl[0], met)

  # Recalculate the upmost tau with new kappa
  tau[0] = kRoss[0]/grav * pl[0]
  # More accurate layer T at uppermost layer
  Tl[0] = 3.0*Tint**4/4.0*(tau[0]+A+B*np.exp(-tau[0]/tau_lim)) + sum(3.0*Beta_V[:] \
    * Tmu**4/4.0*(C[:]+D[:]*np.exp(-tau[0]/tau_lim)+E[:]*np.exp(-gam_V[:]*tau[0])))
  Tl[0] = Tl[0]**(1.0/4.0)

  # Now we can loop in optical depth space to find the T-p profile
  for i in range(1, nlay):
    # Initial guess for layer
    kRoss[i] = k_Ross_Freedman(Tl[i-1], np.sqrt(pl[i-1]*pl[i]), met)
    tau[i] = tau[i-1] + kRoss[i]/grav * (pl[i] - pl[i-1])
    Tl[i] = 3.0*Tint**4/4.0*(tau[i]+A+B*np.exp(-tau[i]/tau_lim)) + sum(3.0*Beta_V[:] \
      * Tmu**4/4.0*(C[:]+D[:]*np.exp(-tau[i]/tau_lim)+E[:]*np.exp(-gam_V[:]*tau[i])))
    Tl[i] = Tl[i]**(1.0/4.0)
    # Convergence loop
    for j in range(5):
      kRoss[i] = k_Ross_Freedman(np.sqrt(Tl[i-1]*Tl[i]), np.sqrt(pl[i-1]*pl[i]), met)
      tau[i] = tau[i-1] + kRoss[i]/grav * (pl[i] - pl[i-1])
      Tl[i] = 3.0*Tint**4/4.0*(tau[i]+A+B*np.exp(-tau[i]/tau_lim)) + sum(3.0*Beta_V[:] \
        * Tmu**4/4.0*(C[:]+D[:]*np.exp(-tau[i]/tau_lim)+E[:]*np.exp(-gam_V[:]*tau[i])))
      Tl[i] = Tl[i]**(1.0/4.0)

  return Tl