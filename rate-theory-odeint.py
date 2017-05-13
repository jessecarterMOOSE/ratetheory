#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# set some parameters
T = 700.0      # temperature in Kelvin
flux = 1.0e14  # flux in n/cm^2/s
rho = 1.0e8    # dislocation density in cm^-2

# set some material parameters
Di0 = 3.5e-4         # interstitial diffusion pre-exponential factor in cm^2/s
Dv0 = 2.2e-2         # vacancy diffusion pre-exponential factor in cm^2/s
Eim = 0.06           # interstitial migration energy in eV
Evm = 0.93           # vacancy migration energy in eV
Evf = 1.79           # vacancy formation energy in eV
Ziv = 500.0          # combinatorial factor for recombination (unitless)
a = 3.234e-8         # lattice constant in cm
Rd = 7.0e-8          # dislocation core radius in cm
N = 4.29e22          # atomic density in atoms/cm^3
dpa_norm = 1.6e-21   # specific dpa rate to fluence ratio
casc_eff = 1.0/3.0

# constants
k_boltz = 8.617e-5   # boltzmann constant in eV/K

# calculate rate-theory parameters
Di = Di0*np.exp(-Eim/k_boltz/T)  # interstitial diffusion coefficient in cm^2/s
Dv = Dv0*np.exp(-Evm/k_boltz/T)  # vacancy diffusion coefficient in cm^2/s
# K_recom = Ziv*Di/a/a/N     # recombination coefficient in cm^3/s
K_recom = 1.0e-8*Ziv*Di        # recombination coefficient in cm^3/s
Q = 2.0*np.pi*rho/np.log(1.0/Rd/np.sqrt(np.pi*rho))  # sink strength for dislocation line in cm^-2
P = dpa_norm*N*flux*casc_eff

# initial conditions
y0 = np.zeros(2)

# times to solve at
t = np.logspace(-8,2,num=1000)
t = np.insert(t,0,np.zeros(1))

# set up RHS of ODE
def ode_system(y, t, K0, Kiv, Q, Di, Dv):
    Ci, Cv = y
    return np.array([K0 - Kiv*Ci*Cv - Q*Di*Ci, K0 - Kiv*Ci*Cv - Q*Dv*Cv])

# solve!
sol = odeint(ode_system, y0, t, args=(P, K_recom, Q, Di, Dv))
Ci = sol[:,0]
Cv = sol[:,1]

# print other info to screen
print 'T: {0:.2f}'.format(T)
print 'flux: {0:.2e}'.format(flux)
print 'rho: {0:.2e}'.format(rho)
print 'P: {0:.2e}'.format(P)
print 'Di: {0:.2e}'.format(Di)
print 'Dv: {0:.2e}'.format(Dv)
print 'K_recom: {0:.2e}'.format(K_recom)
print 'Q: {0:.2e}'.format(Q)
print 'final Ci: {0:.2e}'.format(Ci[-1])
print 'final Cv: {0:.2e}'.format(Cv[-1])

# plot!
plt.loglog(t, Ci, label='Ci')
plt.loglog(t, Cv, label='Cv')
plt.legend()
plt.grid()
plt.xlabel('time (sec)')
plt.ylabel('defect concentration (defects/$\mathdefault{cm^3}$)')
plt.show()
