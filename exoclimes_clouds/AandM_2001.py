'''
A simple interpretation and coding of the Ackerman & Marley (2001) scheme.
To find the particle settling radius, we use the same method as VIRGA, but include potential for Epstein regime particle sizes.


Main functions: 
  AandM_2001 - atmosphere stepper function
Aux functions:
  v_f - function for the settling radius root finder
'''

import numpy as np
from atm_module import vapour_pressure
from scipy.optimize import root_scalar

Scloud = 0.0
R = 8.31446261815324e7

def v_f(r_c, grav, rho_d, rho, eta, mfp, cT, w):

  # Settling velocity function for the root finding algorithm

  # Knudsen number
  Kn = mfp/r_c
  Kn_b = np.minimum(Kn, 100.0)

  # Cunningham slip factor (Jung et al. 2012)
  beta = 1.0 + Kn_b*(1.165 + 0.480 * np.exp(-1.001/Kn_b))

  # Stokes regime (Kn << 1) settling velocity (Ohno & Okuzumi 2017)
  v_f_St = (2.0 * beta * grav * r_c**2 * (rho_d - rho))/(9.0 * eta) \
    * (1.0 + ((0.45*grav*r_c**3*rho*rho_d)/(54.0*eta**2))**(0.4))**(-1.25)

  # Epstein regime (Kn >> 1) regime settling velocity (Woitke & Helling 2003)
  v_f_Ep = (np.sqrt(np.pi)*grav*rho_d*r_c)/(2.0*cT*rho)

  # tanh interpolation function for Kn ~ 1
  fx = 0.5 * (1.0 - np.tanh(2.0*np.log10(Kn)))

  # Interpolation for settling velocity
  v_f = fx*v_f_St + (1.0 - fx)*v_f_Ep

  # Settling velocity to fit is fall velocity minus vertical velocity
  w_diff = v_f - w

  return w_diff


def AandM_2001(nlay, vap_VMR, vap_mw, cld_sp, fsed, sigma, alpha, rho_d, cld_mw, grav, altl, Tl, pl, met, al, Hp, Kzz, mu, eta, rho, cT, mfp):


  # Step through atmosphere to calculate condensate fraction at each layer using A&M 2001 Eq. (8)

  q_v = np.zeros(nlay)
  q_c = np.zeros(nlay)
  q_t = np.zeros(nlay)
  q_s = np.zeros(nlay)

  # Set lowest boundary values - convert vapour VMR to MMR
  q_c[-1] = 0.0
  q_v[-1] = vap_VMR * vap_mw/mu[-1]
  q_t[-1] = q_v[-1]
  p_vap, q_s[-1] = vapour_pressure(vap_mw, cld_sp, Tl[-1], pl[-1], rho[-1], met)

  # Set through atmosphere from lower to upper boundary and perform calculation
  for k in range(nlay-2,-1,-1):
    # First calculate saturation mmr of species (qs) at layer
    p_vap, q_s[k] = vapour_pressure(vap_mw, cld_sp, Tl[k], pl[k], rho[k], met)
    # If the condensate mmr vapour is > than saturation mmr then it 
    # condenses in this layer, else it does not
    if (q_v[k+1] <= q_s[k]):
      # No condensation - solution is not changed from layer below
      q_v[k] = q_v[k+1]
      q_c[k] = 0.0
      q_t[k] = q_v[k] + q_c[k]
      continue
    else:
      # Condensation is triggered - adjust the condensate mixing ratio
      # according to equilibrium solution
      
      # Use A&M Eq. (7) to get total mixing ratio in layer
      dz = altl[k] - altl[k+1]
      L = al * (Hp[k] + Hp[k+1])/2.0
      #q_t[k] = q_t[k+1] * np.exp(-fsed * dz / L)

      q_t[k] = q_s[k] + (q_t[k+1] - q_s[k]) * np.exp(-fsed * dz / L)

      # Use A&M Eq. (8) to get condensate fraction
      q_c[k] = np.maximum(0.0,q_t[k] - (Scloud + 1.0)*q_s[k])

      # Account for any differences between vapour and condensate molecular weight
      q_c[k] = q_c[k]

      q_v[k] = q_t[k] - q_c[k] # (should be = qs[k] by definition)

  # We now have the condensation profile, next we use
  # the balance equation to estimate the cloud properties
  # The upward diffusion of total condensate must equal the downward velocity of the condensate
 
  # Find the vertical convective velocity using Kzz = w * alpha*Hp (Marley & Robinson 2015)
  w = np.zeros(nlay)
  w[:] = Kzz[:]/(al*Hp[:])

  r_w = np.zeros(nlay)
  r_m = np.zeros(nlay)
  r_eff = np.zeros(nlay)
  N_c = np.zeros(nlay)
  for k in range(nlay):
    if (q_c[k] < 1e-10):
      # If low condensate mass fraction, assume zero
      r_w[k] = 0.0
      r_m[k] = 0.0
      N_c[k] = 0.0
    else:

      # Follow the exact same `og' method VIRGA uses, VIRGA has other, probably better, methods included.

      # Range of particle sizes to optimise for [cm]
      r_low = 1.0e-7
      r_up = 10.0

      # Find the optimal radius value that minimises the v_f function
      opt_val = root_scalar(v_f, bracket=[r_low, r_up], method='brentq', 
        args=(grav, rho_d, rho[k], eta[k], mfp[k], cT[k], w[k]))

      # Settling particle radius
      r_w[k] = opt_val.root

      # Median particle radius given log-normal distribution
      r_m[k] = r_w[k] * fsed**(1.0/alpha) * np.exp(-(alpha+6.0)/2.0 * np.log(sigma)**2)

      # Effective particle radius given log-normal distribution
      r_eff[k] = r_w[k] * fsed**(1.0/alpha) * np.exp(-(alpha+1.0)/2.0 * np.log(sigma)**2)

      # Total number density given log-normal distribution
      N_c[k] = (3.0 * q_c[k] * rho[k])/(4.0*np.pi*rho_d*r_m[k]**3) \
        * np.exp(-9.0/2.0 * np.log(sigma)**2)
 
  return q_v, q_c, q_t, q_s, r_w, r_m, r_eff, N_c