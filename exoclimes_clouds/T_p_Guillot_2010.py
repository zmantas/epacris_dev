'''
Function for generating a 1D temperature structure using the Guillot (2010) 
semi-grey analytical radiative-equilibrium 

NOTES: 

We use Eq. 
'''
import numpy as np

def Guillot_T_p(nlay, p0, pl, k_v, k_ir, Tint, mu_z, Tirr, grav):

  gam = k_v/k_ir # visual to IR band gamma ratio
  tau0 = k_ir/grav * p0 # IR optical depth to lowest level

  # Since constant opacities, IR band is just scaled to the lowest optical depth
  tau_irl = np.zeros(nlay) 
  tau_irl[:] = (pl[:]/p0 * tau0)

  Tl = np.zeros(nlay)

  if (mu_z <= 1.0e-6):
    # Zenith angle is very small, assume T-p profile is Eddington relation
    Tl[:] = ((3.0/4.0) * Tint**4 * (tau_irl[:] + 2.0/3.0))

  elif (mu_z > 1.0):
    # Zenith angle cannot be > 1
    print('Invalid zenith angle: ', mu_z)
    quit() 
  else:
    # Use Guillot (2010) Eq. (27)
    Tl[:] = ((3.0/4.0) * Tint**4 * (tau_irl[:] + 2.0/3.0))
    Tl[:] = Tl[:] + (mu_z * 3.0 * Tirr**4)/4.0 *  \
      (2.0/3.0 + mu_z/gam + ((gam/(3.0*mu_z)) - mu_z/gam) * np.exp(-gam*tau_irl[:]/mu_z))
      
  Tl[:] = Tl[:]**(1.0/4.0)

  return Tl