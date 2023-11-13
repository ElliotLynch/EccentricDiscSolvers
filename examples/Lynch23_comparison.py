"""
Lynch23_comparison.py

Computes the fiducial modes of Lynch and Dewberry 2023

"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import eccentric_disc.ModeSolvers as msolvers
import eccentric_disc.geometry as egeo
import eccentric_disc.plottools as eplt
import eccentric_disc.eos as eos

# sound speed
cs0=0.05

# referance beta0 (without taper)
beta0=100
IMRI=cs0

bz=1.0
bt=0

# grid size
N1=200
N2=200
#N1,N2 = 100,100

NTOT = N1*N2

amin=1.0
amax=5.0
wtransition=0.1

ain=amin+0.5
aout=amax-0.5

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

def Ma(a):
  return 2.0*np.pi*a

def entropy(X):
  return np.ones_like(X[0])


  # functions with taper for magnetic field
 
def Vz(a):
  Bz0 = bz*IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  ss=0.5*(1.0+np.tanh((a-ain)/wtransition))*(1.0-np.tanh((a-aout)/wtransition))

  OmegaK=(a**-1.5)

  return (Bz0/cs0)*ss*OmegaK

# being lazy
def sech(x): return 1.0/np.cosh(x)

 # to alow use of mathematica CForm:
def Power(x,n): return x**n

def dVz_da(a):
  
  Bz0 = bz*IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  w=wtransition

  return (Bz0/cs0)*(Power(sech((a - aout)/wtransition),2)*Power(sech((a - ain)/wtransition),2)*(np.cosh((aout - ain)/wtransition) + np.sinh((aout - ain)/wtransition))*(-3*wtransition*np.cosh((aout - ain)/wtransition) - 3*wtransition*np.cosh((-2*a + aout + ain)/wtransition) + 4*a*np.sinh((-2*a + aout + ain)/wtransition)))/(8.*Power(a,2.5)*wtransition)


def Vt(a):
  Bt0 = bt*IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  ss=0.5*(1.0+np.tanh((a-ain)/wtransition))*(1.0-np.tanh((a-aout)/wtransition))

  OmegaK=(a**-1.5)

  return (Bt0/cs0)*ss*OmegaK

def dVt_da(a):

  Bt0 = bt*IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  w=wtransition

  return (Bt0/cs0)*(Power(sech((a - aout)/wtransition),2)*Power(sech((a - ain)/wtransition),2)*(np.cosh((aout - ain)/wtransition) + np.sinh((aout - ain)/wtransition))*(-3*wtransition*np.cosh((aout - ain)/wtransition) - 3*wtransition*np.cosh((-2*a + aout + ain)/wtransition) + 4*a*np.sinh((-2*a + aout + ain)/wtransition)))/(8.*Power(a,2.5)*wtransition)


N_A_POINTS = 200000

def maxe_search(solver,maxe,omega0,ea00):

  a_s = 10.0**np.linspace(np.log10(amin),np.log10(amax),N_A_POINTS)

  def objective(x):
    def rigid_inner_obj(self,a):
      return [0.0,x]

    solver.inner_bc = rigid_inner_obj
    solver.omega0 = omega0

    res = solver.solve()
    soly = solver(a_s)

    return np.max(soly[0]) - maxe

  res = opt.newton(objective,ea00)
  
  return res

def limiting_slope(a):
  abar = 0.5*(amax + amin)

  return np.where(a<abar,a-amin,-a+amax)/a


solver=msolvers.EccentricDiscIsoMHD_2DSolver()
solver.forcing_frequency = forcing_freq
solver.csound2 = csound2
solver.dlnHcirc_dlna = dlnHcirc_dlna

solver.Vz=Vz
solver.Vt=Vt
solver.dVz_da=dVz_da
solver.dVt_da=dVt_da

solver.inner_bc = rigid_inner
solver.outer_bc = rigid_outer

solver.is3D = False

print('start')

solver.amin=amin
solver.amax=amax


# Generate Figure 1 of Lynch and Dewberry 2023
a_s = 10.0**np.linspace(np.log10(amin),np.log10(amax),200000)


ea0s=[0.6,0.85,0.9622957544001679]
omega0s=[-0.001,-0.002,-0.01]
emaxs=[0.2,0.35,0.5]

for i in range(len(ea0s)):

  ea0 = maxe_search(solver,emaxs[i],omega0s[i],ea0s[i])


  solver.omega0 = omega0s[i]
  res = solver.solve()

  print('omega = ',res)


  soly = solver(a_s)

  plt.plot(a_s,soly[0],'k-',linewidth=2)


plt.plot(a_s, limiting_slope(a_s),'k:',linewidth=2)

plt.ylabel(r'$e(a)$',fontsize=18)
plt.xlabel(r'$a$',fontsize=18)

plt.xlim([1.0,5.0])
plt.ylim([0.0,0.7])

plt.show()






