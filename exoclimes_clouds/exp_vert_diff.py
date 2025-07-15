'''

'''

import numpy as np
from numba import jit

@jit
def exp_vert_diff(t_step, nlay, Kzz, rho, alte, altl, q, q0):

  CFL = 0.95
  q_min = 1e-99

  nlev = nlay + 1

  # Find Kzz at levels (layer interfaces)
  Kzz_e = np.zeros(nlev)
  Kzz_e[0] = Kzz[0]
  for k in range(nlay-1):
    Kzz_e[k+1] = (Kzz[k] + Kzz[k+1])/2.0
  Kzz_e[-1] = Kzz[-1]

  # Find rho at levels (layer interfaces)
  rho_e = np.zeros(nlev)
  rho_e[0] = rho[0]
  for k in range(nlay-1):
    rho_e[k+1] = (rho[k] + rho[k+1])/2.0
  rho_e[-1] = rho[-1]

  # Find difference in altitude between levels
  dz = np.zeros(nlay)
  for k in range(nlay):
    dz[k] = alte[k] - alte[k+1]

  # Find difference in altitude between midpoint of layers
  dz_m = np.zeros(nlay-1)
  for k in range(nlay-1):
    dz_m[k] = altl[k] - altl[k+1]

  # Prepare time-stepping routine
  # Find minimum timestep that satisfies the CFL condition
  dt_max = t_step
  for k in range(nlay):
    dt_max = np.minimum(dt_max,CFL*(dz[k]**2/(2.0*Kzz[k])))
  dt = dt_max

  # Copy tracers to work array
  qc = np.zeros(nlay)
  qc[:] = q[:]

  # Work arrays
  phit = np.zeros(nlay)
  phil = np.zeros(nlay)
  flux = np.zeros(nlay)


  t_now = 0.0
  while t_now < t_step:

    # Avoid overshooting the timestep value
    if ((t_now + dt) > t_step):
      dt = t_step - t_now


    # Apply tracer upper and lower boundary conditions
    qc[0] = qc[1]
    qc[-1] = q0
    qc[:] = np.maximum(qc[:],q_min)

    # Compute fluxes at each layer
    for k in range(1,nlay-1):
      phit[k] = rho_e[k+1]*Kzz_e[k+1]*(qc[k+1] - qc[k])/dz_m[k]
      phil[k] = rho_e[k]*Kzz_e[k]*(qc[k] - qc[k-1])/dz_m[k-1]
      flux[k] = ((phit[k] - phil[k])/dz[k])/rho[k]

    qc[:] = qc[:] + dt * flux[:]

    t_now += dt

  # Apply tracer upper and lower boundary conditions
  q[0] = qc[1]
  q[-1] = q0
  q[:] = np.maximum(qc[:],q_min)

  return q

