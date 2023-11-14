"""
AdiabaticTests.py

Compute the lowest order modes for adiabatic discs with gamma=5/3 and gamma=2. For the latter compare the two gamma=2 solvers.
These differ in the evaluation of the geometric part of the Hamiltonian. EccentricDiscAdiabatic2DSolver impliments a series
method (valid for real gamma), while EccentricDiscGamma2_2DSolver makes use of an exact evaluation of the orbit averages possible
when gamma is an integer.

"""
import numpy as np
import matplotlib.pyplot as plt
import eccentric_disc.ModeSolvers as msolvers
import eccentric_disc.geometry as egeo
import eccentric_disc.plottools as eplt

# sound speed
cs0=0.05

def rigid_inner(self,a):
  return [0.0,ea0]

def rigid_outer(self,a,y):
  return y[0]

# put powerlaw into begin
def forcing_freq(a,e):
  return 0.0

def csound2(a):
  return cs0*cs0

def dlnHcirc_dlna(a):
  return 1.0  #(-2.0*ab*((acav/a)**zeta)*l0*zeta + a*np.sqrt(ab/a)*(1.0 + 2.0*zeta*((acav/a)**zeta)))/(2.0*a*np.sqrt(ab/a) - 2.0*ab*l0)

# setup general (2D) adiabatic solver
solver=msolvers.EccentricDiscAdiabatic2DSolver()
solver.forcing_frequency = forcing_freq
solver.csound2 = csound2
solver.dlnHcirc_dlna = dlnHcirc_dlna

solver.gamma=5.0/3.0

solver.inner_bc = rigid_inner
solver.outer_bc = rigid_outer

amin = 1.0
amax = 2.0

print('start')

solver.amin=amin
solver.amax=amax

a_s = 10.0**np.linspace(np.log10(amin),np.log10(amax),200000)

# inner bounday eccentricity gradients and initial guesses to omega to compute
ea0s=[0.1,0.2,0.5,0.7,0.8]
omegas=[-0.01,-0.01,-0.01,-0.2141592274082687,-0.2141592274082687]


#ea0=0.5
#ea0=0.8
for i in range(len(ea0s)):

  ea0=ea0s[i]
  solver.omega0=omegas[i]

  #solver.omega0=-0.2692

  #solver.omega0=-0.01
  #solver.omega0=-0.2141592274082687

  res = solver.solve()

  print('solved')
  print(res)

  soly = solver(a_s)

  print('soly')

  plt.plot(a_s,soly[0],'k-',linewidth=2)


plt.ylabel('e(a)')
plt.xlabel('a')

plt.show()

# switch to gamma=2 as we can compare the implementation
# of gamma=2 in EccentricDiscAdiabatic2DSolver and EccentricDiscGamma2_2DSolver
solver.gamma=2.0

# alternate implentation of (2D) gamma=2 adiabatic disc
# based on exact evaluation of the orbit averages.
solver_gamma2=msolvers.EccentricDiscGamma2_2DSolver()
solver_gamma2.forcing_frequency = forcing_freq
solver_gamma2.csound2 = csound2
solver_gamma2.dlnHcirc_dlna = dlnHcirc_dlna

solver_gamma2.inner_bc = rigid_inner
solver_gamma2.outer_bc = rigid_outer

# inner bounday eccentricity gradients and initial guesses to omega to compute
ea0s=[0.1,0.2,0.5,0.7,0.8]
omegas=[-0.01,-0.01,-0.01,-0.2141592274082687,-0.2141592274082687]

omega_actual=[]
omega_actual_gamma2=[]

#ea0=0.5
#ea0=0.8
for i in range(len(ea0s)):

  ea0=ea0s[i]
  solver.omega0=omegas[i]
  solver_gamma2.omega0=omegas[i]

  res = solver.solve()
  res_gamma2 = solver_gamma2.solve()

  omega_actual+=[res]
  omega_actual_gamma2+=[res_gamma2]

  soly = solver(a_s)
  soly_gamma2 = solver_gamma2(a_s)

  plt.plot(a_s,soly[0],'k-',linewidth=2)
  plt.plot(a_s,soly_gamma2[0],'kx--',linewidth=2)

plt.ylabel('e(a)')
plt.xlabel('a')

plt.show()

# show the computed precession frequencies are the same.
plt.plot(omega_actual,omega_actual_gamma2,'k-')
plt.xlabel('EccentricDiscAdiabatic2DSolver')
plt.ylabel('EccentricDiscGamma2_2DSolver')

plt.show()








