"""
Barker16_comparison.py

Compute the lowest order modes from Barker and Ogilvie 2016 and plot orbits, eccentricity and surface density

Again using this for agile development

TODO:

 - correct weird issue with axis ticks
 - allow use of lambda phi orbital coordinate systems

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

solver=msolvers.EccentricDiscIso3DSolver()
solver.forcing_frequency = forcing_freq
solver.csound2 = csound2
solver.dlnHcirc_dlna = dlnHcirc_dlna

solver.inner_bc = rigid_inner
solver.outer_bc = rigid_outer

solver.is3D = False

amin = 1.0
amax = 2.0

print('start')

solver.amin=amin
solver.amax=amax

a_s = 10.0**np.linspace(np.log10(amin),np.log10(amax),200000)

ea0=0.99755

solver.omega0=-0.2692

res = solver.solve()

print('solved')

soly = solver(a_s)

print('soly')

plt.plot(a_s,soly[0],'k-',linewidth=2)

plt.ylabel('E')
plt.xlabel('a')

plt.show()


print('omega :', res)

fig = plt.figure()

ax = []

ax += [fig.add_subplot(2,1,1)]

ax += [fig.add_subplot(2,1,2)]
ax[1] = eplt.ScalorDiscAxes(fig,ax[1].get_position())


ea0s = [0.1,0.25,0.5,0.75,0.9,0.99755]
omega0s = [-0.01,-0.01,-0.01,-0.01,-0.05,-0.2692]

for i in range(len(ea0s)):

  ea0=ea0s[i]

  solver.omega0=omega0s[i]

  res = solver.solve()

  soly = solver(a_s)

  ax[0].plot(a_s,soly[0],'k-',linewidth=2)

ax[0].set_ylabel('E')
ax[0].set_xlabel('a')

ea0=0.99755

solver.omega0=-0.2692

res = solver.solve()

orbits = egeo.EccentricOrbits()
orbits.solver = solver
orbits.coord_system = "CART"

grid = np.meshgrid(np.linspace(amin,amax,200),np.linspace(0,2.0*np.pi,200), indexing = 'ij')


ax[1].to_cart_grid = orbits.X
ax[1].crdgrid = grid
ax[1].subsample=20

#fig, ax = eplt.scalordisc_subplot(1, 1, grid=grid,to_cart_grid=orbits.X,subsample=20)

ax[1].plot_orbits()

#ax[0].get_xaxis().set_visible(False)
#ax[0].get_yaxis().set_visible(False)
#ax[0].set_frame_on(False)

#ax[1].set_title(r'Orbits')
ax[1].set_ylabel('y')
ax[1].set_ylabel('x')

plt.show()


