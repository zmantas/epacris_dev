'''
  Here we use the implicit Radau 5th order method (implicit Runge-Kutta) for time integration, other implicit schemes could be used.
  For example, BDF are typically more efficient than Radau, requiring less RHS calls.

  Cloud microphysics is generally a stiff problem and requires much more computational power compared to the tracer saturation scheme.

  Main functions: 
    two_moment_mixed - prepares input and integrator
    dqdt - finds the RHS of the microphysics scheme for the time-stepping scheme (dq/dt term)
  Aux functions: 
    calc_cond - calculates the change in mass from condensation/evaporation (dmdt [g s^-1])
    calc_hom_nuc - calculate the nucleation rate, using modified classical (homogenous) nucleation theory J_* [cm^-3 s^-1]
    calc_seed_evap - calculate the evaporation rate of seed particle, J_evap [cm^-3 s^-1]
    calc_coag - calculate the Brownian diffusion coagulation rate [cm^3 s^-1]
    calc_coal - calculate the gravitational coalescence rate [cm^3 s^-1]
'''


import numpy as np
from atm_module import vapour_pressure, surface_tension
from scipy.integrate import solve_ivp

kb = 1.380649e-16
R = 8.31446261815324e7
amu = 1.66053906892e-24
Avo = 6.02214076e23

r_seed = 1e-7
V_seed = 4.0/3.0 * np.pi * r_seed**3

def calc_cond(ncld, r_c, Kn, n_v, D, vth, sat, m0, V_mix):

  # Calculates the condensation or evaporation rate of mass from the bulk (dm/dt [g s^-1])

  dmdt = np.zeros(ncld)
  for n in range(ncld):

    if (sat[n] == 1.0):
      dmdt[n] = 0.0
      continue

    # Diffusive limited regime (Kn << 1) [g s-1]
    dmdt_low = 4.0 * np.pi * r_c * D[n] * m0[n] * n_v[n] \
      *  (1.0 - 1.0/sat[n])

    # Free molecular flow regime (Kn >> 1) [g s-1]
    dmdt_high = 4.0 * np.pi * r_c**2 * vth[n] * m0[n] * n_v[n] * 1.0 \
      * (1.0 - 1.0/sat[n])

    # If evaporation, weight rate by current condensed volume ratio (Woitke et al. 2020)
    if (sat[n] < 1.0):
      dmdt_high = dmdt_high * V_mix[n]
      dmdt_low = dmdt_low * V_mix[n]

    # Critical Knudsen number
    Kn_crit = Kn * (dmdt_high/dmdt_low)

    # Kn' (Woitke & Helling 2003)
    Knd = Kn/Kn_crit

    # tanh interpolation function
    fx = 0.5 * (1.0 - np.tanh(2.0*np.log10(Knd)))

    # Mass change rate (First moment)
    dmdt[n] = dmdt_low * fx + dmdt_high * (1.0 - fx)

  return dmdt

def calc_hom_nuc(ncld, n_v, T, p, mw_cld, rho_d, cld_nuc, sig_inf, sat, r0):

  # Uses Modified Classical (Homogenous) Nucleation Theory (MCRT) to calculate the rate of formation of seed particles (J_* [cm^-3 s^-1])

  J_hom = np.zeros(ncld)
  for n in range(ncld):

    # Only calculate J_* for pre-determined nucleation species
    if (cld_nuc[n] == False):
      J_hom[n] = 0.0
      continue

    # If unsaturated, nucleation cannot take place
    if (sat[n] <= 1.0):
      J_hom[n] = 0.0
    else:

      # Nucleation parameters
      alpha = 1.0
      Nf = 0.0
      third = 1.0/3.0
      twothird = 2.0/3.0

      # Efficiency Variables
      ln_ss = np.log(sat[n]) 
      f0 = 4.0 * np.pi * r0[n]**2 
      kbT = kb * T       

      theta_inf = (f0 * sig_inf[n])/(kbT)  
      N_inf = (((2.0/3.0) * theta_inf) / ln_ss)**3

      # N_* variable
      N_star = 1.0 + (N_inf / 8.0) \
        * (1.0 + np.sqrt(1.0 + 2.0*(Nf/N_inf)**third) \
        - 2.0*(Nf/N_inf)**third)**3
      N_star = np.maximum(1.00001, N_star)
      N_star_1 = N_star - 1.0

      dg_rt = theta_inf * (N_star_1 / (N_star_1**third + Nf**third))

      Zel = np.sqrt((theta_inf / (9.0 * np.pi * (N_star_1)**(4.0/3.0))) \
        * ((1.0 + 2.0*(Nf/N_star_1)**1.0/3.0)/(1.0 + (Nf/N_star_1)**third)**3))

      tau_gr = (f0 * N_star**twothird) * alpha * np.sqrt(kbT \
        / (2.0 * np.pi * mw_cld[n] * amu)) * n_v[n]

      J_hom[n] = n_v[n] * tau_gr * Zel * np.exp(np.maximum(-300.0, N_star_1*ln_ss - dg_rt))

  return J_hom

def calc_seed_evap(ncld, N_c, m_c, m_seed, sat, cld_nuc):

  # Calculates the evaporation of the seed particle core, given conditions are `reasonable' (J_evap [cm^-3 s^-1])

  J_evap = np.zeros(ncld)
  for n in range(ncld):

    # Only evaporate species designated as nucleation speces
    if (cld_nuc[n] == False):
      J_evap[n] = 0.0
      continue

    # If saturated then evaporation can't take place
    if (sat[n] >= 1.0):
      J_evap[n] = 0.0

    else:
      # Check if average mass is around 0.1% the seed particle mass
      # This means the core is (probably) exposed to the air and can evaporate freely
      if (m_c <= (1.001 * m_seed[n])):
        # Assume a reasonable evaporation rate (here 0.1 seconds)
        tau_evap = 0.1
        # Seed particle evaporation rate [cm-3 s-1]
        J_evap[n] = -N_c/tau_evap
      else:
        # There is still some mantle to evaporate from
        J_evap[n] = 0.0

  return J_evap

def calc_coag(T, m_c, r_c, beta, eta):

  # Calculates the monodisperse Brownian diffusion coagulation rate (f_coag [cm^3 s^-1])

  # Particle diffusion limit rate (Kn << 1)
  D_r = (kb*T*beta)/(6.0*np.pi*eta*r_c)

  # Thermal velocity limit rate (Kn >> 1)
  V_r = np.sqrt((8.0*kb*T)/(np.pi*m_c))

  # Moran (2022) Kn ~ 1 interpolation method using diffusive Knudsen number
  Knd = (8.0*D_r)/(np.pi*V_r*r_c)
  phi = 1.0/np.sqrt(1.0 + np.pi**2/8.0 * Knd**2)

  # Coagulation rate (Zeroth moment) [cm^3 s^-1]
  f_coag = (4.0*kb*T*beta)/(3.0*eta) * phi

  return f_coag

def calc_coal(grav, r_c, Kn, vf):

  # Calculates the monodisperse gravitational coalescence rate (f_coal [cm^3 s^-1])

  # Epsilon value parameter - approximate the relative settling velocity as a spread of monodisperse settling velocity
  # (Ohno & Okuzumi 2017,Sato et al. 2016)
  eps = 0.5

  # Estimate differential velocity
  d_vf = eps * vf

  # Calculate E - collisional efficiency factor
  if (Kn >= 1.0):
    # E = 1 when Kn >= 1
    E = 1.0
  else:
    # Calculate Stokes number
    Stk = (vf * d_vf)/(grav * r_c)
    # Guillot et al. (2014) expression
    E = np.maximum(0.0,1.0 - 0.4*Stk**(-0.75))

  # Coalescence rate (Zeroth moment) [cm3 s-1]
  f_coal = 2.0*np.pi*r_c**2*d_vf*E

  return f_coal

def dqdt(t, y, ncld, p_vap, vth, D, sig_inf, Rd_v, cld_nuc, nd_atm, rho, mfp, T, p, m_seed, rho_d, m0, cld_mw, r0, eta, grav, cT):

  # Limit y values
  y[:] = np.maximum(y[:],1e-30)

  f = np.zeros(len(y))

  rho_c = np.zeros(ncld)
  rho_v = np.zeros(ncld)
  p_v = np.zeros(ncld)
  n_v = np.zeros(ncld)
  sat = np.zeros(ncld)
  V_mix = np.zeros(ncld)

  # Convert y to real physical numbers to calculate f
  N_c = y[0]*nd_atm
  rho_c[:] = y[1:2+ncld-1]*rho 
  rho_v[:] = y[2+ncld-1:]*rho 

  # Find the partial pressure of the vapour and number density
  p_v[:] = rho_v[:] * Rd_v[:] * T
  n_v[:] = p_v[:]/(kb*T) 

  # Find supersaturation ratio
  sat[:] = p_v[:]/p_vap[:]

  # Total condensed mass [g cm^-3]
  rho_c_t = np.sum(rho_c[:])

  # Mean mass of particle [g]
  m_c = np.maximum(rho_c_t/N_c, m_seed[0])

  # Bulk density of particle mixture [g cm^-3]
  rho_d_m = 0.0
  for n in range(ncld):
    rho_d_m += (rho_c[n]/rho_c_t) * rho_d[n]

  # Mass weighted mean radius of particle [cm]
  r_c = np.maximum(((3.0*m_c)/(4.0*np.pi*rho_d_m))**(1.0/3.0), r_seed)

  # Bulk material volume mixing ratio of mixture
  V_tot = np.sum(rho_c[:]/rho_d[:]) # Total condensed volume
  V_mix[:] = (rho_c[:]/rho_d[:])/V_tot # Condensed volume mixing ratio

  # Knudsen number
  Kn = mfp/r_c
  Kn_b = np.minimum(Kn, 100.0)

  # Cunningham slip factor (Jung et al. 2012)
  beta = 1.0 + Kn_b*(1.165 + 0.480 * np.exp(-1.001/Kn_b))

  # Settling velocity (Stokes regime)
  vf_s = (2.0 * beta * grav * r_c**2 * (rho_d_m - rho))/(9.0 * eta) \
    * (1.0 + ((0.45*grav*r_c**3*rho*rho_d_m)/(54.0*eta**2))**(0.4))**(-1.25)

  # Settling velocity (Epstein regime)
  vf_e = (np.sqrt(np.pi)*grav*rho_d_m*r_c)/(2.0*cT*rho)

  # tanh interpolation function
  fx = 0.5 * (1.0 - np.tanh(2.0*np.log10(Kn)))

  # Interpolation for settling velocity
  vf = fx*vf_s + (1.0 - fx)*vf_e
  vf = np.maximum(vf, 1e-30)

  # Calculate condensation/evaporation rate
  f_cond = calc_cond(ncld, r_c, Kn, n_v, D, vth, sat, m0, V_mix)

  # Calculate homogenous nucleation rate
  f_nuc_hom = calc_hom_nuc(ncld, n_v, T, p, cld_mw, rho_d, cld_nuc, sig_inf, sat, r0)

  # Calculate seed particle evaporation rate
  f_seed_evap = calc_seed_evap(ncld, N_c, m_c, m_seed, sat, cld_nuc)

  # Calculate Brownian coagulation rate
  f_coag = calc_coag(T, m_c, r_c, beta, eta)

  # Calculate gravitational coalescence rate
  f_coal = calc_coal(grav, r_c, Kn, vf)

  # Calculate final net flux rate for each moment and vapour
  f[0] = np.sum(f_nuc_hom[:] + f_seed_evap[:]) - (f_coag + f_coal)*N_c**2
  f[1:2+ncld-1] = m_seed[:]*(f_nuc_hom[:]  + f_seed_evap[:]) + f_cond[:]*N_c
  f[2+ncld-1:] = -f[1:2+ncld-1]

  # Convert f to ratios
  f[0] = f[0]/nd_atm
  f[1:] = f[1:]/rho

  # Check if condensation from vapour is viable
  for n in range(ncld):
    if (rho_v[n]/rho <= 1e-28):
      if (f[1+n] > 0.0):
        f[1+n] = 0.0
        f[1+ncld+n] = 0.0

  # Check if evaporation from condensate is viable
  for n in range(ncld):
    if (rho_c[n]/rho <= 1e-28):
      if (f[1+ncld+n] > 0.0):
        f[1+n] = 0.0
        f[1+ncld+n] = 0.0

  return f

def two_moment_mixed(nlay, ncld, t_step, vap_VMR, vap_mw, cld_sp, rho_d, cld_mw, Tl, pl, nd_atm, rho, mfp, eta, cT, mu, grav, met, cld_nuc, q_v, q_0, q_1):

  # work arrays
  p_vap = np.zeros(ncld)
  q_s = np.zeros((nlay, ncld))
  vth = np.zeros(ncld)
  D = np.zeros(ncld)
  sig_inf = np.zeros(ncld)
  Rd_v = np.zeros(ncld)

  # Calculate come fundamental properties, such as unit masses, volumes etc
  m0 = np.zeros(ncld)
  V0 = np.zeros(ncld)
  r0 = np.zeros(ncld)
  d0 = np.zeros(ncld)
  m_seed = np.zeros(ncld)
  for n in range(ncld):
    m0[n] = cld_mw[n] * amu
    V0[n] = m0[n] / rho_d[n]
    r0[n] = ((3.0*V0[n])/(4.0*np.pi))**(1.0/3.0)
    d0[n] = 2.0 * r0[n]
    Rd_v[n] = R/vap_mw[n]
    if (cld_nuc[n] == False):
      m_seed[n] = 0.0
    else:
      m_seed[n] = V_seed * rho_d[n]

  # Prepare integration solver
  y0 = np.zeros(1+ncld+ncld) # one q_0, ncld q_1 and ncld q_v

  # Tolerances and time-stepping
  rtol = 1e-3
  atol = 1e-99
  max_step = np.inf
  t_span = [0.0, t_step]

  for k in range(nlay):

    # Calculate constant variables for each species
    for n in range(ncld):

      # Saturation vapour mass mixing ratio for this layer
      p_vap[n], q_s[k,n] = vapour_pressure(vap_mw[n], cld_sp[n], Tl[k], pl[k], rho[k], met)

      # Thermal velocity of vapour [cm s-1]
      vth[n] = np.sqrt((kb*Tl[k])/(2.0*np.pi*vap_mw[n]*amu))

      # Gaseous diffusion constant of vapour [cm2 s-1] 
      D[n] = 5.0/(16.0*Avo*d0[n]**2*rho[k]) * np.sqrt((R*Tl[k]*mu[k])/(2.0*np.pi) * (vap_mw[n] + mu[k])/vap_mw[n])

      # Surface tension of species [erg cm-2]
      sig_inf[n] = surface_tension(cld_sp[n], Tl[k])

    # Give tracer values to y array
    y0[0] = q_0[k]
    y0[1:2+ncld-1] = q_1[k,:]
    y0[2+ncld-1:] = q_v[k,:]

    # Use implicit stiff method to integrate the tracer values in time
    sol = solve_ivp(dqdt, t_span, y0, method='Radau', rtol=rtol, atol=atol, \
      args=(ncld, p_vap, vth, D, sig_inf, Rd_v, cld_nuc, nd_atm[k], rho[k], mfp[k], Tl[k], pl[k], m_seed, rho_d, m0, cld_mw, r0, eta[k], grav, cT[k]))

    # Give back results to the vapour and condensate array
    q_0[k] = sol.y[0,-1]
    q_1[k,:] = sol.y[1:2+ncld-1,-1]
    q_v[k,:] = sol.y[2+ncld-1:,-1]

    q_0[k] = np.maximum(q_0[k], 1e-30)
    for n in range(ncld):
      q_1[k,n] = np.maximum(q_1[k,n], 1e-30)
      q_v[k,n] = np.maximum(q_v[k,n], 1e-30)
      q_s[k,n] = np.maximum(q_s[k,n], 1e-30)

  return q_v, q_0, q_1, q_s