"""
initial_conditions_for_PLUTO2p5D.py

Computes initial conditions for a hydro, 3D unstratified, PLUTO simulation

Based on code provided by Loren E. Held

"""
import numpy as np
import matplotlib.pyplot as plt
import eccentric_disc.ModeSolvers as msolvers
import eccentric_disc.geometry as egeo
import eccentric_disc.plottools as eplt
import eccentric_disc.eos as eos

from scipy.interpolate import CubicSpline

# sound speed
cs0=0.05

# grid size
N1=128
N2=128

NTOT = N1*N2


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

def ecc_to_cyl(aa,ee,aea,Nr=128,Np=128):
    ''' 
        takes in a profile of eccentricity e=e(a) as a function of semimajor
        axis, and returns the associated density and velocity fields on 
        uniformly spaced meshes in cylindrical radius and phi.
    '''
    rr,ph = np.linspace(aa[0],aa[-1],Nr), np.linspace(0,2*np.pi,Np)

    # semilatus rectum:
    lam = aa*(1. - ee**2)

    # gradient (quick and dirty):
    lep = lam*np.gradient(ee,lam)

    #lep = aea*(1 - ee*ee)/(1.0 - ee*ee - 2.0*ee*aea)


    # 2D meshes in lam,phi coordinates:
    Lm,Ss = np.meshgrid(lam,np.sin(ph),indexing='ij')
    Ee,Cc = np.meshgrid(ee,np.cos(ph),indexing='ij')
    Le,Ph = np.meshgrid(lep,ph,indexing='ij')

    Aea,Ph=np.meshgrid(aea,ph,indexing='ij')

    #cosE=(Cc + Ee)/(1.0+Ee*Cc)

    #j = (1.0 - Ee*(Ee + Aea))/np.sqrt(1.0 - Ee**2) - Aea*cosE/np.sqrt(1.0 - Ee**2)


    # cylindrical radius R=R(lam,phi) 
    RR = Lm/(1. + Ee*Cc)
    RR = np.minimum(RR,aa[-1]) # can get floating pt R>a[-1] otherwise
    RR = np.maximum(RR,aa[0])

    # primitive variables (from eqns in Barker & Ogilvie, 2016), on a grid 
    # in which R=R(a,phi) depends both on lam and phi:
    Sg = (1. - Ee**2)**(3./2.)*(1. + Ee*Cc)/(1. + (Ee - Le)*Cc)
    #Sg = (1. + Ee*Cc)/(1. + (Ee - Le)*Cc)
    #Sg = 1.0/j
    Ur = Ee*Ss/np.sqrt(Lm)
    Up = (1. + Ee*Cc)/np.sqrt(Lm)

    # interpolate to put primitive vars on a grid mesh with R independent from phi:
    R_cyl,P_cyl = np.meshgrid(rr,ph,indexing='ij')
    Sg_cyl = np.array([CubicSpline(x=RR[:,j],y=Sg[:,j])(rr) for j in range(Np)]).T
    Ur_cyl = np.array([CubicSpline(x=RR[:,j],y=Ur[:,j])(rr) for j in range(Np)]).T
    Up_cyl = np.array([CubicSpline(x=RR[:,j],y=Up[:,j])(rr) for j in range(Np)]).T

    return R_cyl,P_cyl,Sg_cyl,Ur_cyl,Up_cyl



solver=msolvers.EccentricDiscIso3DSolver()
solver.forcing_frequency = forcing_freq
solver.csound2 = csound2
solver.dlnHcirc_dlna = dlnHcirc_dlna

solver.inner_bc = rigid_inner
solver.outer_bc = rigid_outer

solver.is3D = False

amin = 1.0
amax = 5.0

print('start')

solver.amin=amin
solver.amax=amax

a_s = 10.0**np.linspace(np.log10(amin),np.log10(amax),200000)

fig = plt.figure()

ax = []

ax += [fig.add_subplot(2,1,1)]

ax += [fig.add_subplot(2,1,2)]
ax[1] = eplt.ScalorDiscAxes(fig,ax[1].get_position())

#ea0 = 0.8459459371420034

ea0=0.35

solver.omega0 = -0.001
res = solver.solve()

print('omega = ', res)

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


R=cylgrid[0]
phi=cylgrid[1]

print(R[:,0])

#es=disc.e(R[:,0])
#qs=disc.q(R[:,0])

#emesh, dummy=np.meshgrid(es,np.linspace(0,2.0*np.pi,N2), indexing = 'ij')
#qmesh, dummy=np.meshgrid(qs,np.linspace(0,2.0*np.pi,N2), indexing = 'ij')

#meanmotion=R**(-1.5)

#vx1=R*meanmotion*emesh*np.sin(phi)
#vx2=R*meanmotion*(1.0 + 2.0*emesh*np.cos(phi))
#vx3=np.zeros_like(vx1)

#rho=1.0 + qmesh*np.cos(phi)

rho, vx1, vx2, vx3 = disc.save(grid,cylgrid,Nz=16,filename="data.dbl")


r_cyl, phi_cyl, rho_cyl, vx1_xyl, vx2_cyl = ecc_to_cyl(aa=a_s,ee=soly[0],aea=a_s*soly[1],Nr=128,Np=128)


plt.imshow(rho[0,:,:])
plt.colorbar()

plt.show()

plt.imshow(rho_cyl.T)
plt.colorbar()

plt.show()

plt.imshow((rho[0,:,:]-rho_cyl.T)/rho[0,:,:])
plt.colorbar()

plt.show()


plt.plot(R[:,0],vx2[0,0,:]/R[:,0],'k-')
plt.plot(R[:,0],R[:,0]**(-1.5),'k--')
plt.show()

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

plt.show()


#PLOT 2D SLICES IN xy-PLANE:
plt.figure(5)
plt.xlabel('x')
plt.ylabel('y')
plt.imshow(rho[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.colorbar()
plt.title('density (xy-plane)')

plt.figure(6)
plt.imshow(vx1[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('vx1 (xy-plane)')

plt.figure(7)
plt.imshow(vx2[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('vx2 (xy-plane)')

plt.figure(8)
plt.imshow(vx3[0,:,:],origin='lower',cmap='RdBu_r')
plt.axes().set_aspect('auto')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.title('vx3 (xy-plane)')

plt.show()





