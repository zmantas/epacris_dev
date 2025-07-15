'''
  Here we use the Dormand-Price 5(4) Runge-Kutta method for time integration, which is accurate, but probably overkill.
  Simpler methods can be used such as Euler (see commented out section), or lower/higher order Runge-Kutta methods.

  Tracer saturation adjustment should generally not be stiff if the tau_cond value is reasonable, 
  if stiffness is encountered (for some reason) then you would need to switch to an implicit ODE solver, which would be much slower.

  Main functions: 
    tracer_sat_adj - prepares input and integrator
    dqdt - finds the RHS of the tracer saturation adjustment scheme for the time-stepping scheme (dq/dt term)
  Aux functions: 
    None
'''

import numpy as np
from atm_module import vapour_pressure
from scipy.integrate import solve_ivp

def dqdt(t, y, q_s_k, tau_cond):

  # y is length two - index 0 is vapour, index 1 is condensate
  #y[0] = np.maximum(y[0],1e-30)
  #y[1] = np.maximum(y[1],1e-30)

  # RHS array
  f = np.zeros(len(y))

  if (y[0] < q_s_k):
    # Vapour is undersaturated - adjust vapour by evaporating from q_c
    f[0] = np.minimum(q_s_k - y[0], y[1])/tau_cond
  elif (y[0] > q_s_k):
    # Vapour is supersaturated - adjust vapour by condensing from q_v
    f[0] = -(y[0] - q_s_k)/tau_cond
  else:
    # Do nothing as q_v = q_s
    f[0] = 0.0

  # RHS of condensate is negative vapour RHS
  f[1] = -f[0]

  return f

def tracer_sat_adj(nlay, t_step, vap_VMR, vap_mw, cld_sp, rho_d, cld_mw, Tl, pl, rho, met, tau_cond, q_v, q_c):

  # Saturation vapour mass mixing ratio at each layer
  q_s = np.zeros(nlay)
    
  # Prepare integration solver
  y0 = np.zeros(2)

  # Tolerances and time-stepping
  rtol = 1e-3
  atol = 1e-30
  max_step = np.inf
  t_span = [0.0, t_step]

  for k in range(nlay):

    # Saturation vapour mass mixing ratio for this layer
    p_vap, q_s[k] = vapour_pressure(vap_mw, cld_sp, Tl[k], pl[k], rho[k], met)

    # Initial y value array
    y0[0] = q_v[k]
    y0[1] = q_c[k]

    # Perform sub-time-stepping integration  (alternative handmade simple Euler method)
    # t_now = 0.0
    # t_sub = tau_cond # Usually a suitable timestep is the condensation timescale
    # while t_now < t_step:
    #   # Avoid overshooting the timestep value
    #   if ((t_now + t_sub) > t_step):
    #     t_sub = t_step - t_now
    #   f = dqdt(t_now, y0, q_s[k], tau_cond)
    #   y0[:] = y0[:] + t_sub * f[:]
    #   t_now += t_sub
    # q_v[k] = y0[0]
    # q_c[k] = y0[1]

    # Use explicit Runge-Kutta method to integrate the tracer values in time
    sol = solve_ivp(dqdt, t_span, y0, method='RK45', rtol=rtol, atol=atol, args=(q_s[k],tau_cond))

    # Give back results to the vapour and condensate array
    q_v[k] = sol.y[0,-1]
    q_c[k] = sol.y[1,-1]

    q_v[k] = np.maximum(q_v[k], 1e-30)
    q_c[k] = np.maximum(q_c[k], 1e-30)
    q_s[k] = np.maximum(q_s[k], 1e-30)

  return q_v, q_c, q_s