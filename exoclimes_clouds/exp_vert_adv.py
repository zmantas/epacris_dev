'''

'''

import numpy as np
from numba import jit

@jit
def exp_vert_adv(t_step, nlay, v_f, rho, alte, altl, q_in, q0):


  CFL = 0.95
  q_min = 1e-99

  nlev = nlay + 1

  # Find v_f at levels (layer interfaces)
  v_f_e = np.zeros(nlev)
  v_f_e[0] = v_f[0]
  for k in range(nlay-1):
    v_f_e[k+1] = (v_f[k] + v_f[k+1])/2.0
  v_f_e[-1] = v_f[-1]


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
    dt_max = np.minimum(dt_max,CFL*(dz[k]/(v_f[k])))
  dt = dt_max

  # Copy tracers to work array
  q = np.zeros(nlay)
  q[:] = q_in[:] * rho[:]

  # Work arrays
  c = np.zeros(nlay-1)
  qc = np.zeros(nlay)
  phi = np.zeros(nlay)

  eps = 1e-10

  t_now = 0.0
  while t_now < t_step:

    # Avoid overshooting the timestep value
    if ((t_now + dt) > t_step):
      dt = t_step - t_now

     # Find the courant number
    c[:] = (np.absolute(v_f_e[1:nlay]) * dt) / dz_m[:]

    # Calculate flux limiter using koren method
    phi[0] = 0.0
    for k in range(1, nlay-1):
      dq_minus = (q[k] - q[k-1]) / dz_m[k-1]
      dq_plus  = (q[k+1] - q[k]) / dz_m[k]
      r = dq_minus / (dq_plus + eps)
      phi[k] = np.maximum(0.0, np.minimum(np.minimum(2.0 * r, (1.0 + 2.0 * r) / 3.0), 2.0))
    phi[-1] = 0.0

    # Perform McCormack step - do not evolve index nlay as we have fixed boundary condition
    # Predictor step (forward in space)
    qc[:] = q[:]
    qc[:-1] = q[:-1] - phi[:-1] * c[:] * (q[1:] - q[:-1])

    # Corrector step (backward in space)
    q[1:nlay-1] = 0.5 * (q[1:nlay-1] + qc[1:nlay-1] - c[1:nlay-1] * (qc[1:nlay-1] - qc[0:nlay-2]))

    q[:] = np.maximum(q[:],q_min)
    q[0] = q[1]
    q[-1] = q0

    t_now += dt

  # Apply tracer upper and lower boundary conditions
  q[:] = np.maximum(qc[:]/rho[:],q_min)
  q[0] = q[1]
  q[-1] = q0

  return q