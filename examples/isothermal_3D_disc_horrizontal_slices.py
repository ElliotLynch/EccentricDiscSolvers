"""
isothermal_3D_disc_horrizontal_slices.py

plots horrizontal slice through fluid properties of a 3D isothermal disc

"""
import numpy as np
import matplotlib.pyplot as plt
import eccentric_disc.ModeSolvers as msolvers
import eccentric_disc.geometry as egeo
import eccentric_disc.plottools as eplt
import eccentric_disc.eos as eos
import eccentric_disc.VerticalStructure as evert

# sound speed
cs0=0.05

H0 = 0.05

# grid size
N1=200
N2=200
N3=100

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
  return np.ones_like(a)

def entropy(X):
  return np.ones_like(X[0])

def Hcirc(a):

  return H0*np.ones_like(a)

solver=msolvers.EccentricDiscIso3DSolver()
solver.forcing_frequency = forcing_freq
solver.csound2 = csound2
solver.dlnHcirc_dlna = dlnHcirc_dlna

solver.inner_bc = rigid_inner
solver.outer_bc = rigid_outer

solver.is3D = True

amin = 1.0
amax = 2.0

#print 'start'

solver.amin=amin
solver.amax=amax

ea0 = 0.75

solver.omega0 = -0.01
res = solver.solve()

vsolver = evert.VerticalStructureIsothermal()
vsolver.scale_heightcirc = Hcirc
vsolver.coord_system = 'AE' #I think

disc = egeo.EccentricDisc()
disc.solver = solver
disc.coord_system = "CART"
disc.eos = eos.EOSIsothermal()
disc.Ma = Ma
disc.entropy = entropy

vsolver.orbits = disc

disc.vertical_structure = vsolver

Zmax = 0.25

alin = np.linspace(amin,amax,N1)
Elin = np.linspace(0.0,2.0*np.pi,N2)
Zlin = np.linspace(-Zmax,Zmax,N3)

z1 = Zlin[N3//2]/H0
z2 = Zlin[N3//2+2]/H0
z3 = Zlin[N3//2+20]/H0


#raise

grid = np.meshgrid(alin,Elin,Zlin, indexing = 'ij') 
plot_grid = np.meshgrid(alin,Elin, indexing = 'ij')

fig = plt.figure()

ax = []

ax += [fig.add_subplot(2,3,1)]

ax[0] = eplt.ScalorDiscAxes(fig,ax[0].get_position())
ax[0].crdgrid = plot_grid
ax[0].subsample=15
ax[0].to_cart_grid = disc.X

#fig, ax = eplt.scalordisc_subplot(1, 3, grid=plot_grid,to_cart_grid=disc.X,subsample=15)

pressure = disc["PRS"](grid)



#print 'AX: ', ax[0]

 # density

ax[0].plot_orbits()
ax[0].pcolormesh(Z = pressure[:,:,N3//2])
ax[0].set_ylabel('y')
ax[0].set_xlabel('x')

ax[0].set_title(r'$z/H^{\circ} = $'+str(z1))
ax[0].set_colorbar(orientation='horizontal')

  # density



ax += [fig.add_subplot(2,3,2)]

ax[1] = eplt.ScalorDiscAxes(fig,ax[1].get_position())
ax[1].crdgrid = plot_grid
ax[1].subsample=15
ax[1].to_cart_grid = disc.X

ax[1].plot_orbits()
ax[1].pcolormesh(Z = pressure[:,:,N3//2+2])
ax[1].set_ylabel('y')
ax[1].set_xlabel('x')

ax[1].set_title(r'$z/H^{\circ} = $'+str(z2))
ax[1].set_colorbar(orientation='horizontal')


ax += [fig.add_subplot(2,3,3)]

ax[2] = eplt.ScalorDiscAxes(fig,ax[2].get_position())
ax[2].crdgrid = plot_grid
ax[2].subsample=15
ax[2].to_cart_grid = disc.X

 # density

ax[2].plot_orbits()
ax[2].pcolormesh(Z = pressure[:,:,N3//2+20])
ax[2].set_ylabel('y')
ax[2].set_xlabel('x')

ax[2].set_title(r'$z/H^{\circ} = $'+str(z3))
ax[2].set_colorbar(orientation='horizontal')


ax += [fig.add_subplot(2,3,5)]


disc.coord_system='CART'
X, Y, Z = disc.X(grid)

p3 = ax[3].pcolormesh(X[:,N2//2,:],Z[:,N2//2,:],pressure[:,N2//2,:])

ax[3].set_xlabel('X')
ax[3].set_ylabel('Z')

fig.colorbar(p3, ax = ax[3])

#ax[3].colorbar()

plt.show()


#plt.contourf(X[:,N2//2,:],Z[:,N2//2,:],pressure[:,N2//2,:])

#plt.xlabel('X')
#plt.ylabel('Z')

#plt.colorbar()

#plt.show()


