"""
initial_conditions_for_PLUTO2p5D.py

Computes initial conditions for a 3D unstratified PLUTO simulation

Based on code provided by Loren E. Held

"""
import numpy as np
import matplotlib.pyplot as plt
import eccentric_disc.ModeSolvers as msolvers
import eccentric_disc.geometry as egeo
import eccentric_disc.plottools as eplt
import eccentric_disc.eos as eos

# sound speed
cs0=0.05

# referance beta0 (without taper)
beta0=100
IMRI=cs0

bt_to_bz=0.0

# grid size
N1=128
N2=128

NTOT = N1*N2

amin=1.0
amax=5.0
wtransition=0.05

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
  Bz0 = IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  ss=0.5*(1.0+np.tanh((a-ain)/wtransition))*(1.0-np.tanh((a-aout)/wtransition))

  OmegaK=(a**-1.5)

  return (Bz0/cs0)*ss*OmegaK

# being lazy
def sech(x): return 1.0/np.cosh(x)

 # to alow use of mathematica CForm:
def Power(x,n): return x**n

def dVz_da(a):
  
  Bz0 = IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  w=wtransition

  return (Bz0/cs0)*(Power(sech((a - aout)/wtransition),2)*Power(sech((a - ain)/wtransition),2)*(np.cosh((aout - ain)/wtransition) + np.sinh((aout - ain)/wtransition))*(-3*wtransition*np.cosh((aout - ain)/wtransition) - 3*wtransition*np.cosh((-2*a + aout + ain)/wtransition) + 4*a*np.sinh((-2*a + aout + ain)/wtransition)))/(8.*Power(a,2.5)*wtransition)


def Vt(a):
  Bt0 = bt_to_bz*IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  ss=0.5*(1.0+np.tanh((a-ain)/wtransition))*(1.0-np.tanh((a-aout)/wtransition))

  OmegaK=(a**-1.5)

  return (Bt0/cs0)*ss*OmegaK

def dVt_da(a):

  Bt0 = bt_to_bz*IMRI/(2.0*np.pi*np.sqrt(16.0/15.0))

  w=wtransition

  return (Bt0/cs0)*(Power(sech((a - aout)/wtransition),2)*Power(sech((a - ain)/wtransition),2)*(np.cosh((aout - ain)/wtransition) + np.sinh((aout - ain)/wtransition))*(-3*wtransition*np.cosh((aout - ain)/wtransition) - 3*wtransition*np.cosh((-2*a + aout + ain)/wtransition) + 4*a*np.sinh((-2*a + aout + ain)/wtransition)))/(8.*Power(a,2.5)*wtransition)


solver=msolvers.EccentricDiscIsoMHD_2DSolver()
solver.forcing_frequency = forcing_freq
solver.csound2 = csound2
solver.dlnHcirc_dlna = dlnHcirc_dlna
#solver.beta0=100.0

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

a_s = 10.0**np.linspace(np.log10(amin),np.log10(amax),200000)

fig = plt.figure()

ax = []

ax += [fig.add_subplot(2,1,1)]

ax += [fig.add_subplot(2,1,2)]
ax[1] = eplt.ScalorDiscAxes(fig,ax[1].get_position())

ea0 = 0.3

solver.omega0 = -0.001
res = solver.solve()

soly = solver(a_s)

ax[0].plot(a_s,soly[0],'k-',linewidth=2)

ax[0].set_ylabel('E')
ax[0].set_xlabel('a')

disc = egeo.EccentricDisc()
disc.solver = solver
disc.eos = eos.EOSIsothermal()
disc.Ma = Ma
disc.entropy = entropy
disc.is3D = False

grid = np.meshgrid(np.linspace(amin,amax,N1),np.linspace(0,2.0*np.pi,N2), indexing = 'ij')

ax[1].to_cart_grid = disc.X
ax[1].crdgrid = grid
ax[1].subsample=20

ax[1].plot_orbits()

ax[1].set_ylabel('y')
ax[1].set_xlabel('x')

plt.show()

disc.coord_system = "CYL"

Rmin = amin
Rmax = amax

cylgrid = np.meshgrid(np.linspace(Rmin,Rmax,N1),np.linspace(0,2.0*np.pi,N2), indexing = 'ij')

agrid=grid[0]

#den = interp_to_cylgrid(disc,'DEN',agrid,cylgrid[0],cylgrid[1])


#plt.pcolormesh(cylgrid[0],cylgrid[1],den)

#plt.show()


#print np.max(disc.e(agrid))
#print np.max(disc.q(agrid))

#print np.max(disc.ey(agrid))

#print np.max(disc.jacobian(*grid)), np.min(disc.jacobian(*grid))

#raise

rho, vx1, vx2, vx3, bx1, bx2, bx3 = disc.save(grid,cylgrid,Nz=16,isMHD=True,filename="data.dbl")

# plots showing slices throught the initial conditions

plt.figure(1)
plt.imshow(rho[:,0,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('z')
plt.colorbar()
plt.title('density (xz-plane)')

plt.figure(2)
plt.imshow(vx1[:,0,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('z')
plt.colorbar()
plt.title('vx1 (xz-plane)')

plt.figure(3)
plt.imshow(vx2[:,0,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('z')
plt.colorbar()
plt.title('vx2 (xz-plane)')

plt.figure(4)
plt.xlabel('x')
plt.ylabel('z')
plt.imshow(vx3[:,0,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.colorbar()
plt.title('vx3 (xz-plane)')

plt.figure(5)
plt.imshow(bx1[:,0,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('z')
plt.colorbar()
plt.title('bx1 (xz-plane)')

plt.figure(6)
plt.imshow(bx2[:,0,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('z')
plt.colorbar()
plt.title('bx2 (xz-plane)')

plt.figure(7)
plt.xlabel('x')
plt.ylabel('z')
plt.imshow(bx3[:,0,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.colorbar()
plt.title('bx3 (xz-plane)')


plt.show()


#PLOT 2D SLICES IN xy-PLANE:
plt.figure(8)
plt.xlabel('x')
plt.ylabel('y')
plt.imshow(rho[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.colorbar()
plt.title('density (xy-plane)')

plt.figure(9)
plt.imshow(vx1[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('vx1 (xy-plane)')

plt.figure(10)
plt.imshow(vx2[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('vx2 (xy-plane)')

plt.figure(11)
plt.imshow(vx3[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('vx3 (xy-plane)')

plt.figure(12)
plt.imshow(bx1[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('bx1 (xy-plane)')

plt.figure(13)
plt.imshow(bx2[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('bx2 (xy-plane)')

plt.figure(14)
plt.imshow(bx3[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('bx3 (xy-plane)')

plt.show()





