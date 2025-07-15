'''
Example 3: Mass moment method, including mixed materials

Steps:

1. Generate a T-p profile using semi-grey or picket-fence semi-analytical solutions.
2. Timestep the vertical settling and diffusion scheme with mass moment method, given parameters in parameters.yaml.
3. Generate output plots and png.
'''

import numpy as np # for efficient array numerical operations 
import matplotlib.pylab as plt # for plotting
import seaborn as sns # for good colourblind colour scheme
import yaml # for reading in the input parameter YAML file

from atm_module import hypsometric, visc_mixture, adiabat_correction, v_f_two_moment

from T_p_Guillot_2010 import Guillot_T_p # Import function for Guillot 2010 semi-grey profile
from T_p_Parmentier_2015 import Parmentier_T_p # Import function for Parmentier & Guillot 2015 picket-fence profile

from two_moment_mixed import two_moment_mixed
from exp_vert_adv import exp_vert_adv
from exp_vert_diff import exp_vert_diff

R = 8.31446261815324e7
kb = 1.380649e-16
amu = 1.66053906892e-24

# Open parameter YAML file and read parameters for A&M profile
with open('parameters.yaml', 'r') as file:
  param = yaml.safe_load(file)['two_moment_mixed']

# Now extract the parameters from the YAML file into local variables

# Give number of layers in 1D atmosphere - number of levels is nlay + 1
nlay = param['nlay']
nlev = nlay + 1

mu_z = param['mu_z']
Tirr = param['Tirr']
Tint = param['Tint']
k_v =  param['k_v']
k_ir = param['k_ir']
met = param['met']
grav = param['grav']
kappa = param['kappa']

# Get top and bottom pressure layers in bar - convert to dyne
ptop = param['ptop'] * 1e6
pbot = param['pbot'] * 1e6

# Get pressure at levels (edges) - assume log-spaced
pe = np.logspace(np.log10(ptop),np.log10(pbot),nlev)

# Get pressure at layers using level spacing
pl = np.zeros(nlay)
pl[:] = (pe[1:] - pe[0:-1])/np.log(pe[1:]/pe[0:-1])
 

# Kzz profile - here we assume constant but a 1D array for easy change to
# non-constant or some function with some editing
# A&M (2001) use a convective heat flux mixing length theory prescription
Kzz = np.zeros(nlay)
Kzz[:] = param['Kzz']
# Example change - use expression from Parmentier et al. (2013):
#Kzz[:] = 5e8/np.sqrt(pl[:]*1e-6)

# We do the same for molecular weight of the atmosphere 
# Can be changed to varying with height/pressure
mu = np.zeros(nlay)
mu[:] = param['mu']

# Get 1D analytical temperature pressure profile at layers
Tl = np.zeros(nlay)
if (param['Guillot'] == True):
  Tl[:] = Guillot_T_p(nlay, pbot, pl, k_v, k_ir, Tint, mu_z, Tirr, grav)
elif (param['Parmentier'] == True):
  Tl[:] = Parmentier_T_p(nlay, pl, Tint, mu_z, Tirr, grav, met, 1)
else:
  print('Invalid T structure selection')
  quit()

if (param['adibat_corr'] == True):
  Tl[:] = adiabat_correction(nlay, Tl, pl, kappa)

# Atmosphere mass density
rho = np.zeros(nlay)
rho[:] = (pl[:]*mu[:]*amu)/(kb * Tl[:])

# Atmosphere number density
nd_atm = np.zeros(nlay)
nd_atm[:] = (pl[:])/(kb * Tl[:])

# Atmosphere thermal velocity
cT = np.zeros(nlay)
cT[:] = np.sqrt((2.0 * kb * Tl[:]) / (mu[:] * amu))

# Find the altitude grid and scale heights at levels and layers using the hypsometric equation
alte = np.zeros(nlev)
alte, Hp = hypsometric(nlev, Tl, pe, mu, grav)
altl = (alte[0:-1] + alte[1:])/2.0

# Volume mixing ratio of condensate vapour at lower boundary
# Deep abyssal mixing ratio (assume infinite supply at constant value)
bg_sp = param['bg_sp']
bg_VMR = param['bg_VMR']
bg_mw = param['bg_mw']
bg_d = param['bg_d']
bg_LJ = param['bg_LJ']
nbg = len(bg_sp)

# Find dynamical viscosity of each layer given a background gas mixture
eta = np.zeros(nlay)
for k in range(nlay):
  eta[k] = visc_mixture(Tl[k], nbg, bg_VMR, bg_mw, bg_d, bg_LJ)

# Find mean free path of each layer
mfp = np.zeros(nlay)
mfp[:] = (2.0*eta[:]/rho[:]) * np.sqrt((np.pi * mu[:])/(8.0*R*Tl[:]))

cld_sp = param['cld_sp']
rho_d = param['rho_d']
cld_mw = param['cld_mw']
cld_nuc = param['nuc_sp']

vap_mw =  param['vap_mw']
vap_VMR = param['vap_VMR']

ncld = len(cld_sp)

# Lower boundary conditions for vapour and condensate
q_v_bot = np.zeros(ncld)
for n in range(ncld):
  q_v_bot[n] = vap_VMR[n] * vap_mw[n]/mu[-1]

q_0_bot = 1e-30
q_1_bot = np.zeros(ncld)
q_1_bot[:] = 1e-30

n_step = param['n_step']
n_step = int(n_step)
t_step = param['t_step']

q_v = np.zeros((nlay,ncld))
q_0 = np.zeros(nlay) 
q_1 = np.zeros((nlay,ncld))
q_s = np.zeros((nlay,ncld))
v_f = np.zeros(nlay)

# Being timestepping loop for n_step
t_now = 0.0
for t in range(n_step):

  # Perform cloud microphysics in each layer
  q_v[:,:], q_0[:], q_1[:,:], q_s[:,:] = \
    two_moment_mixed(nlay, ncld, t_step, vap_VMR, vap_mw, cld_sp, rho_d, cld_mw, Tl, pl, nd_atm, rho, mfp, eta, cT, mu, grav, met, cld_nuc, q_v, q_0, q_1)

  # Calculate vertical settling velocity for each moment in each layer
  v_f[:] = v_f_two_moment(nlay, ncld, q_0, q_1, grav, rho_d, nd_atm, rho, eta, mfp, cT)

  # Advect condensate tracer (downwards) - give condensate boundary conditions
  q_0[:] = exp_vert_adv(t_step, nlay, v_f, rho, alte, altl, q_0, q_0_bot)
  for n in range(ncld):
    q_1[:,n] = exp_vert_adv(t_step, nlay, v_f, rho, alte, altl, q_1[:,n], q_1_bot[n])

  # Diffuse vapour and condensate tracer (upward and downward) - give vapour and condensate boundary conditions
  q_0[:] = exp_vert_diff(t_step, nlay, Kzz, rho, alte, altl, q_0, q_0_bot)
  for n in range(ncld):
    q_1[:,n] = exp_vert_diff(t_step, nlay, Kzz, rho, alte, altl, q_1[:,n], q_1_bot[n])
    q_v[:,n] = exp_vert_diff(t_step, nlay, Kzz, rho, alte, altl, q_v[:,n], q_v_bot[n])

  # Limit values to a maximum
  q_0[:] = np.maximum(q_0[:], 1e-30)
  for n in range(ncld):
    q_1[:,n] = np.maximum(q_1[:,n], 1e-30)
    q_v[:,n] = np.maximum(q_v[:,n], 1e-30)

  t_now += t_step

  print(t, t_now)

# After time-stepping, find N_c and r_c at each layer
N_c = np.zeros(nlay)
N_c[:] = q_0[:] * nd_atm[:]

rho_c = np.zeros((nlay,ncld))
rho_c_t = np.zeros(nlay)
for n in range(ncld):
  rho_c[:,n] = q_1[:,n]*rho[:]
  rho_c_t[:] += rho_c[:,n]

m_c = np.zeros(nlay)
m_c[:] = rho_c_t[:]/N_c[:]

rho_d_m = np.zeros(nlay)
for n in range(ncld):
  rho_d_m[:] += (rho_c[:,n]/rho_c_t[:]) * rho_d[n]

r_c = np.zeros(nlay)
r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d_m[:]))**(1.0/3.0), 1e-7)

V_mix = np.zeros((nlay,ncld))
V_mix_t = np.zeros(nlay)
for k in range(k):
  V_mix_t = 0.0
  for n in range(ncld):
    V_mix_t +=  rho_c[k,n]/rho_d[n]
  V_mix[k,:] = (rho_c[k,:]/rho_d[:])/V_mix_t

fig = plt.figure() # Start figure 
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

Tp = ax2.plot(Tl,pl/1e6,c='black',label=r'T-p',ls='dotted')

colour = sns.color_palette('colorblind') # Decent colourblind wheel (10 colours)

for n in range(ncld):
  if (n == 0):
    qv = ax1.plot(q_v[:,n],pl/1e6,c=colour[n],ls='dashed',lw=2,label=r'$q_{\rm v}$')
    q0 = ax1.plot(q_0[:],pl/1e6,c=colour[n],ls='dashdot',lw=2,label=r'$q_{\rm 0}$')
    q1 = ax1.plot(q_1[:,n],pl/1e6,c=colour[n],ls='solid',lw=2,label=r'$q_{\rm 1}$')
    qs = ax1.plot(q_s[:,n],pl/1e6,c=colour[n],ls='dotted',lw=1,label=r'$q_{\rm s}$')
  else:
    ax1.plot(q_v[:,n],pl/1e6,c=colour[n],ls='dashed',lw=2)
    ax1.plot(q_1[:,n],pl/1e6,c=colour[n],ls='solid',lw=2)
    ax1.plot(q_s[:,n],pl/1e6,c=colour[n],ls='dotted',lw=1)    

ax1.set_xlabel(r'$q$ [g g$^{-1}$]',fontsize=16)
ax2.set_xlabel(r'$T$ [K]',fontsize=16)
ax1.set_ylabel(r'$p$ [bar]',fontsize=16)
ax1.tick_params(axis='both',which='major',labelsize=14)
ax2.tick_params(axis='both',which='major',labelsize=14)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1e-8,1e-2)
ax2.set_zorder(1)
lns = Tp + qv + q0 + q1 + qs
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs,fontsize=10,loc='upper right')
plt.gca().invert_yaxis()
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('example_3_q.png',dpi=300,bbox_inches='tight')

fig = plt.figure() # Start figure 

colour = sns.color_palette('colorblind') # Decent colourblind wheel (10 colours)

plt.plot(N_c[:],pl/1e6,c=colour[n],ls='solid',lw=2,label=r'$N_{\rm c}$')
plt.xlabel(r'$N_{\rm c}$ [cm$^{-3}$]',fontsize=16)
plt.ylabel(r'$p$ [bar]',fontsize=16)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('example_3_N_c.png',dpi=300,bbox_inches='tight')


fig = plt.figure() # Start figure 

colour = sns.color_palette('colorblind') # Decent colourblind wheel (10 colours)

plt.plot(r_c[:]*1e4,pl/1e6,c=colour[n],ls='solid',lw=2,label=r'$r_{\rm c}$')
plt.xlabel(r'$r$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$p$ [bar]',fontsize=16)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('example_3_r_c.png',dpi=300,bbox_inches='tight')

fig = plt.figure() # Start figure 

colour = sns.color_palette('colorblind') # Decent colourblind wheel (10 colours)

plt.plot(v_f[:],pl/1e6,c=colour[n],ls='solid',lw=2,label=r'$v{\rm f}$')
plt.xlabel(r'$v_{\rm f}$ [cm s$^{-1}$]',fontsize=16)
plt.ylabel(r'$p$ [bar]',fontsize=16)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('example_3_v_f.png',dpi=300,bbox_inches='tight')


fig = plt.figure() # Start figure 

colour = sns.color_palette('colorblind') # Decent colourblind wheel (10 colours)

for n in range(ncld):
  if (n == 0):
    plt.plot(V_mix[:,n],pl/1e6,c=colour[n],ls='solid',lw=2,label=r'$V_{\rm mix}$')
  else:
    plt.plot(V_mix[:,n],pl/1e6,c=colour[n],ls='solid',lw=2)
plt.xlabel(r'$V_{\rm mix}$',fontsize=16)
plt.ylabel(r'$p$ [bar]',fontsize=16)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('example_3_V_mix.png',dpi=300,bbox_inches='tight')


plt.show()