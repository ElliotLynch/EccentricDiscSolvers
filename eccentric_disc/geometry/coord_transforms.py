import numpy as np
import scipy.optimize as opt

####
# contain function/class to convert between caresian,cylindrical and orbital coordinates
#
####



def cart2cyl(X):
  x = X[0]
  y = X[1]

  is3D = False

  if (len(X)==3):
    is3D = True
    z = X[2]

  r = np.sqrt(x*x+y*y)
  phi = np.arctan2(y,x)

  if is3D:
    return r, phi, z
  else:
    return r, phi


def cyl2cart(X):

  r = X[0]
  phi = X[1]

  is3D = False

  if (len(X)==3):
    is3D = True
    z = X[2]
  
  x = r*np.cos(phi)
  y = r*np.sin(phi)

  if is3D:
    return x, y, z
  else:
    return x, y

 # check which vertical coordinate to use

def aE2cyl(orbits,X):

  a = X[0] 
  eccanom = X[1]

  is3D = False

  if (len(X)==3):
    is3D = True
    z = X[2]

  eccentricity = orbits.e(a)

  r = a*(1.0 - eccentricity*np.cos(eccanom))
  
  xrot = a*(np.cos(eccanom) - eccentricity)

  yrot = a*np.sqrt(1.0 - eccentricity*eccentricity)*np.sin(eccanom)

  # probably not the sensible way of doing this
  if type(orbits.omega(a))==type(None):
    phi = np.arctan2(yrot,xrot)
  else:
    phi = np.arctan2(yrot,xrot) + orbits.omega(a)

  if is3D:
    return r, phi, z
  else:
    return r, phi


def aE2cart(orbits,X):
  return cyl2cart(aE2cyl(orbits,X))

  #rewrite to allow for arrays

def cyl2aE(X,orbits):

  # assume flattened for not, implement the general case
  # after

  r = X[0]
  phi = X[1]

  is3D = False

  if (len(X)==3):
    is3D = True
    z = X[2]

  a = []

   # assume r is an array
   # deal with array like later 
  for i in xrange(r.size):

    def objective(a):
      #r = a*(1.0 - orbits.e(a)*np.cos(eccanom))
  
      f = phi[i] - orbits.omega(a)
 
      return r[i] - a*(1.0 - (orbits.e(a)**2))/(1.0 + orbits.e(a)*np.cos(f))

  # check which solver we want to use
   
    a += [opt.root(objective,r[i])]

  a = np.array(a)

  f = phi - orbits.omega(a)

  xrot = r*np.cos(f)
  yrot = r*np.sin(f)

  #xrot + a*orbits.e(a) = a*np.cos(eccanom)
  #yrot/np.sqrt(1.0 - (orbits.e(a)**2)) = a*np.sin(eccanom)

  eccanom = np.arctan2(yrot/np.sqrt(1.0 - (orbits.e(a)**2)),xrot + a*orbits.e(a))

  if is3D:
    return a, eccanom, z
  else:
    return a, eccanom

def cart2aE(orbits,X):
  return aE2cyl(orbits,*cart2cyl(X))

 # helper function as we use this alot
 # doing this the lazy way
def phi2EccentricAnomaly(a,phi,e,omega=0,increasingEccanom=True):

  if omega==None:
    f=phi
  else:
    f=phi-omega

  r = a*(1 - e*e)/(1 + e*np.cos(f))

  xrot = r*np.cos(f)
  yrot = r*np.sin(f)

  # check this makes sence
  eccanom = np.arctan2(yrot/np.sqrt(1.0 - e*e),xrot + a*e)

  if increasingEccanom:
    eccanom=np.where(eccanom>=0,eccanom,eccanom+2.0*np.pi)

  return eccanom

def aphi2aE(orbits,X,increasingEccanom=True):
 
  a=X[0]
  phi=X[1]

  is3D=False
  if (len(X)==3):
    is3D = True
    z=X[2]

  eccanom=np.zeros_like(a)

  # this isn't going to work with the new way we deal with EccentricOrbits

  for j in xrange(a.shape[1]):
    if is3D:
      for k in xrange(a.shape[2]):
        eccanom[:,j,k] = phi2EccentricAnomaly(a[:,j,k],phi[:,j,k],orbits.e(a[:,j,k]),orbits.omega(a[:,j,k]),increasingEccanom=increasingEccanom)
    else:
      eccanom[:,j] = phi2EccentricAnomaly(a[:,j],phi[:,j],orbits.e(a[:,j]),orbits.omega(a[:,j]),increasingEccanom=increasingEccanom)

  if is3D:
    return a, eccanom, z
  else:
    return a, eccanom


if __name__ == '__main__':

  pass
