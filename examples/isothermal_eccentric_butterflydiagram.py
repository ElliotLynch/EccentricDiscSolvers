"""
isothermal_eccentric_butterflydiagram

plots butterfly diagrams (iso-velocity surfaces) for (3D) isothermal eccentric disc


"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
from scipy.integrate import simps
import eccentric_disc.ModeSolvers as msolvers
import eccentric_disc.geometry as egeo
import eccentric_disc.plottools as eplt
import eccentric_disc.eos as eos
import eccentric_disc.VerticalStructure as evert

# sound speed
cs0=0.05

H0 = 0.05
h0=H0

GM=1.0
#h0=0.1

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

def EquilibriumEccentricFlux(x,y,z,nu,nustar):

  R = np.sqrt(x*x+y*y+z*z)

  flux_normal = nu*nu*nu/(np.exp(nu*np.sqrt(R)/nustar) - 1.0)

  return flux_normal

solver=msolvers.EccentricDiscIso3DSolver()
solver.forcing_frequency = forcing_freq
solver.csound2 = csound2
solver.dlnHcirc_dlna = dlnHcirc_dlna

solver.inner_bc = rigid_inner
solver.outer_bc = rigid_outer

solver.is3D = True

amin = 1.0
amax = 10.0

#print 'start'

solver.amin=amin
solver.amax=amax

ea0 = 0.75

solver.omega0 = -0.1
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

# grid size
N1=1000
N2=1000
N3=3


Zmax = 0.25

alin = 10.0**np.linspace(np.log10(amin),np.log10(amax),N1)
Elin = np.linspace(0.0,2.0*np.pi,N2)
Zlin = np.linspace(-Zmax,Zmax,N3)

#z1 = Zlin[N3//2]/H0
#z2 = Zlin[N3//2+2]/H0
#z3 = Zlin[N3//2+20]/H0

#print "max e = ", np.max(np.abs(disc.e(alin)))

grid = np.meshgrid(alin,Elin,Zlin, indexing = 'ij')
plot_grid = np.meshgrid(alin,Elin, indexing = 'ij')

agrid=grid[0]

print("max e = ", np.max(np.abs(disc.e(agrid))))

#den = disc["DEN"](grid)
#vx, vy, vz = disc.V(grid)


disc.coord_system='CART'
X, Y, Z = disc.X(grid)

 # hopefully works?

den = disc["DEN"](grid)
vx, vy, vz = disc.V(grid)

ht = disc["H"](grid)

#p3 = plt.pcolormesh(X[:,N2//2,:],Z[:,N2//2,:],den[:,N2//2,:])

print(den[:,:,N3//2])

plt.pcolormesh(X[:,:,N3//2],Y[:,:,N3//2],den[:,:,N3//2])  #,vmin=0.1,vmax=10.0)

plt.xlabel('X')
plt.ylabel('Y')

#fig.colorbar(p3, ax = ax[3])

#ax[3].colorbar()

plt.show()



plt.pcolormesh(X[:,:,N3//2],Y[:,:,N3//2],ht[:,:,N3//2])  #,vmin=0.1,vmax=10.0)
plt.colorbar()


plt.xlabel('X')
plt.ylabel('Y')

#fig.colorbar(p3, ax = ax[3])

#ax[3].colorbar()

plt.show()

#inc = np.pi/2.0-0.2
inc=np.pi/4.0
ztld=1.0

Z = ztld*ht
#Zl = -ztld*ht

 # not clear this is the easiest way of doing this,

Vz = disc["dH"](grid)*ztld 
#Vzl = -disc["dH"](grid)*ztld

x=X*np.cos(inc)-Z*np.sin(inc)
#xl=X*np.cos(inc)-Zl*np.sin(inc)
y=Y

VN = vx*np.cos(inc) + Vz*np.sin(inc)
#VNl = vx*np.cos(inc) + Vzl*np.sin(inc)

plt.pcolormesh(x[:,:,N3//2],y[:,:,N3//2],VN[:,:,N3//2])  #,vmin=0.1,vmax=10.0)

plt.xlabel('X')
plt.ylabel('Y')

plt.show()

#plt.pcolormesh(rgrid,phigrid,den)

#plt.show()

nustar=1.0

nu0=0.01
dnu=0.001

nus=np.linspace(nu0-dnu,nu0+dnu,50)

flux=np.zeros((agrid.shape[0],agrid.shape[1]))

for i in range(agrid.shape[0]):
  for j in range(agrid.shape[1]):
    flux[i,j]=simps(EquilibriumEccentricFlux(X[i,j,N3//2],Y[i,j,N3//2],Z[i,j,N3//2],nus,nustar),nus)

V0s=np.linspace(-0.3,0.3,6)

fig, axs = plt.subplots(2, 3)

for i in range(2):
  for j in range(3):

    V0=V0s[i*3+j]
#V0=0.1
    dV=0.01

    butterfly_basic=flux*(np.where(np.abs(VN-V0)<dV,1,0)[:,:,N3//2])

    axs[i,j].pcolormesh(x[:,:,N3//2],y[:,:,N3//2],np.nan_to_num(butterfly_basic))

#print flux

#plt.pcolormesh(x[:,:,N3//2],y[:,:,N3//2],np.nan_to_num(butterfly_basic)) #,vmin=0.0,vmax=1.0e-7)

#plt.xlabel('X')
#plt.ylabel('Y')

plt.show()




def hcirc(a):
  return h0*a

def dH_isolin(a,eccanom,e):

  n = np.sqrt(GM/(a*a*a))

  return 3.0*hcirc(a)*e*n*np.sin(eccanom)/(1 - e*np.cos(eccanom))
  #return 0.0

def H_isolin(a,eccanom,e):

  return hcirc(a)*(1.0 - 3.0*e*np.cos(eccanom))

  #return hcirc(a)

def getVN(a,eccanom,e,omega,inc,dH,ztld):

  n = np.sqrt(GM/(a*a*a))

  return -a*n*np.cos(inc)*(np.sin(eccanom)*np.cos(omega) + np.sqrt(1 - e*e)*np.cos(eccanom)*np.sin(omega))/(1 - e*np.cos(eccanom)) + dH*ztld*np.sin(inc)

  # assume H is given as above

def geta_xy(x,y,e,omega,inc,scale_height,ztld):

  #x = X*np.sin(inc) - Z*np.cos(inc)

  #X=x/np.sin(inc)
  #Y=y

  #phi = np.arctan2(Y,X)

  def objective(cvec):

    a=cvec[0]
    Z=cvec[1]

    X = (x + Z*np.cos(inc))/np.sin(inc)
    Y=y

    phi = np.arctan2(Y,X)

    #H = scale_height(a,eccanom,e(a))

    f = phi - omega(a)
 
    cosE=(np.cos(f) + e(a))/(1 + e(a)*np.cos(f))

    #H = scale_height(a,eccanom,e(a))

    #dH = dH_isolin(a,e(a),eccanom)
    #return VN-getVN(a,e(a),omega(a),eccanom,inc,dH,ztld)
    r=np.sqrt(X*X+Y*Y)     
 
    xrot = r*np.cos(f)
    yrot = r*np.sin(f)

    # check this makes sence
    eccanom = np.arctan2(yrot/np.sqrt(1.0 - e(a)*e(a)),xrot + a*e(a))

    H = scale_height(a,eccanom,e(a))

    return r-a*(1.0 - e(a)*cosE), Z-ztld*H 

  X0=x/np.sin(inc)
  Y0=y

  r0=np.sqrt(X0*X0+Y0*Y0)

  sol = opt.root(objective, [r0,0.0])
  a = sol.x[0]
  Z = sol.x[1] # do we care about Z? - yes
  
  X = (x + Z*np.cos(inc))/np.sin(inc)
  Y=y

  phi = np.arctan2(Y,X)

  f = phi - omega(a)

  # this is incorrect
  #eccanom=np.arccos((np.cos(f) + e(a))/(1 + e(a)*np.cos(f)))
 
  # copied from the eccentric disc module - combine at some point

  xrot = r0*np.cos(f)
  yrot = r0*np.sin(f)

  # check this makes sence
  eccanom = np.arctan2(yrot/np.sqrt(1.0 - e(a)*e(a)),xrot + a*e(a))

  return a, eccanom, Z

def getaEgrid(x,y,e,omega,inc,scale_height,ztld):

  a_s=np.zeros_like(x)
  Es=np.zeros_like(x) 
  Zs = np.zeros_like(x)

  for i in range(x.shape[0]):
    for j in range(x.shape[1]):
      a_s[i,j], Es[i,j], Zs[i,j] =geta_xy(x[i,j],y[i,j],e,omega,inc,scale_height,ztld)


  return a_s, Es, Zs

  # this should have an orientation/solid angle dependance
def EquilibriumEccentricFlux(a,eccanom,e,nu,nustar):

  flux_normal = nu*nu*nu/(np.exp(nu*np.sqrt(a*(1 - e*np.cos(eccanom)))/nustar) - 1.0)

  return flux_normal


omega0=1.0
inc= np.pi/4.0
e0=0.1 #2
#e0=0.0


def omega_aligned(a):
  return omega0

def efall(a):
  #return e0/a

  return e0*np.sin(np.pi*(a-1.0)/10.0)

ztld=1.0

xs=np.linspace(-10,10,200)
ys=np.linspace(-10,10,200)

xmesh, ymesh = np.meshgrid(xs,ys,indexing='ij')

amesh_u, Emesh_u, Zmesh_u = getaEgrid(xmesh,ymesh,efall,omega_aligned,inc,H_isolin,ztld)
amesh_l, Emesh_l, Zmesh_l = getaEgrid(xmesh,ymesh,efall,omega_aligned,inc,H_isolin,-ztld)

plt.pcolormesh(xmesh,ymesh,Emesh_u,vmin=-np.pi,vmax=np.pi)
plt.show()




#print np.min(Emesh), np.max(Emesh)

#raise

dH_u=dH_isolin(amesh_u,Emesh_u,efall(amesh_u))
dH_l=dH_isolin(amesh_l,Emesh_l,efall(amesh_l))

VNUpper=getVN(amesh_u,Emesh_u,efall(amesh_u),omega_aligned(amesh_u),inc,dH_u,ztld)
VNLower=getVN(amesh_l,Emesh_l,efall(amesh_l),omega_aligned(amesh_l),inc,dH_l,-ztld)

print(VNUpper)


plt.pcolormesh(xmesh,ymesh,VNUpper,vmin=-1.0,vmax=1.0)
plt.show()

nustar=1.0

nu0=0.01
dnu=0.001

nus=np.linspace(nu0-dnu,nu0+dnu,50)

flux_u=np.zeros_like(amesh_u)
flux_l=np.zeros_like(amesh_l)

for i in range(amesh_u.shape[0]):
  for j in range(amesh_u.shape[1]):
    flux_u[i,j]=simps(EquilibriumEccentricFlux(amesh_u[i,j],Emesh_u[i,j],efall(amesh_u[i,j]),nus,nustar),nus)
    flux_l[i,j]=simps(EquilibriumEccentricFlux(amesh_l[i,j],Emesh_l[i,j],efall(amesh_l[i,j]),nus,nustar),nus)

#print flux

plt.pcolormesh(xmesh,ymesh,flux_u,vmin=0.0,vmax=1.0e-7)
plt.show()

V0s=np.linspace(-0.6,0.6,12)

fig, axs = plt.subplots(3, )

for i in range(3):
  for j in range(4):

    V0=V0s[i*4+j]
#V0=0.1
    dV=0.01

    butterfly_basic=flux_u*np.where(np.abs(VNUpper-V0)<dV,1,0) + flux_l*np.where(np.abs(VNLower-V0)<dV,1,0)

    axs[i,j].pcolormesh(xmesh,ymesh,np.nan_to_num(butterfly_basic))

plt.show()


