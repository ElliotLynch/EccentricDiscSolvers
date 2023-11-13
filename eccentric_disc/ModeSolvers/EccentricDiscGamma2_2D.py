import numpy as np


# defining these so we can make use
# of CForm
def sec(x): return 1.0/np.cos(x)

def csc(x): return 1.0/np.sin(x)

def Power(x,n): return x**n

# Using Trigonometric regulerisation from Ogilvie and Lynch (2019)
# to deal with the aparent singularity at e=f


 # geometric part of the Hamiltonian
def geoFGamma2_2D(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #this is probably wrong
  #return np.log(np.cos(2.0*alpha)) - 2.0*np.log(np.cos(alpha+beta)) + np.sin(2.0*alpha)*np.sin(alpha-beta)*sec(alpha+beta)

  # should be
  return ((-np.cos(3*alpha - beta) + 3*np.cos(alpha + beta))*sec(2*beta)*sec(alpha + beta))/2.0

# various derivatives of the geometric part of the Hamiltonian

def dgeoFGamma2_2Dde(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return ((-1 + np.cos(4*alpha) + 2*np.cos(2*(alpha - beta)))*sec(2*alpha)*Power(sec(alpha + beta),2)*np.tan(2*alpha))/4.0

  return (sec(2*alpha)*sec(2*beta)*Power(sec(alpha + beta),2)*(np.sin(4*alpha) + 2*np.sin(2*alpha - 2*beta)))/4.

def dgeoFGamma2_2Ddf(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0
 
  #return -(sec(2*beta)*Power(sec(alpha + beta),2)*(np.sin(4*alpha) - 2*np.sin(2*(alpha + beta))))/4.0

  return (Power(sec(2*beta),3)*Power(sec(alpha + beta),2)*(-4*np.sin(2*alpha) + np.sin(2*alpha - 4*beta) + 6*np.sin(2*beta) - 2*np.sin(4*alpha + 2*beta) + 3*np.sin(2*alpha + 4*beta)))/8.

def ddgeoFGamma2_2Ddf2(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return ((3*np.cos(alpha - beta) - np.cos(3*(alpha - beta)) + np.cos(3*alpha + beta) + 
  #     np.cos(alpha + 3*beta) + np.cos(5*alpha + 3*beta) - np.cos(3*alpha + 5*beta))*Power(sec(2*beta),3)*Power(sec(alpha + beta),3))/8.0

  return ((np.cos(alpha - 7*beta) - 13*np.cos(alpha - 3*beta) - 12*np.cos(3*alpha - beta) + 37*np.cos(alpha + beta) + 15*np.cos(3*(alpha + beta)) + 3*np.cos(5*(alpha + beta)) - 3*np.cos(5*alpha + beta) - 9*np.cos(alpha + 5*beta) - 3*np.cos(3*alpha + 7*beta))*Power(sec(2*beta),5)*Power(sec(alpha + beta),3))/16.

def ddgeoFGamma2_2Ddfde(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return (sec(2*beta)*Power(sec(alpha + beta),3)*(3*np.sin(alpha - beta) + np.sin(3*alpha + beta))*np.tan(2*alpha))/4.0

  return ((np.cos(alpha - 5*beta) - np.cos(alpha - beta) - 3*np.cos(3*alpha + beta) - np.cos(5*alpha + 3*beta))*sec(2*alpha)*Power(sec(2*beta),3)*Power(sec(alpha + beta),3))/8.

# the following are anti-aligned version of the above functions for use with toroidal magnetic fields

def geoFGamma2_2D_AntiAligned(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return sec(2.0*beta)*sec(alpha+beta)*np.cos(2.0*beta)*np.cos(alpha-beta)

  return np.cos(2*alpha)*np.cos(alpha - beta)*sec(2*beta)*sec(alpha + beta)

def dgeoFGamma2_2Dde_AntiAligned(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return -0.25*sec(2.0*alpha)*sec(2.0*beta)*(sec(alpha+beta)**2)*(np.sin(4.0*alpha) + 2.0*np.sin(2.0*(alpha-beta)))

  #doing from cform
  return -0.25*(sec(2*alpha)*sec(2*beta)*Power(sec(alpha + beta),2)*(np.sin(4*alpha) + 2*np.sin(2*alpha - 2*beta)))

def dgeoFGamma2_2Ddf_AntiAligned(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return 0.25*np.cos(2.0*alpha)*(sec(2.0*beta)**3)*(sec(alpha+beta)*2)*(np.sin(4.0*alpha) + 2.0*np.sin(2.0*(alpha+beta)))

  return (np.cos(2*alpha)*Power(sec(2*beta),3)*Power(sec(alpha + beta),2)*(np.sin(4*beta) + 2*np.sin(2*(alpha + beta))))/4.

def ddgeoFGamma2_2Ddf2_AntiAligned(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return ((-np.cos(alpha - 7*beta) + np.cos(alpha - 3*beta) + 8*np.cos(3*alpha - beta) + 11*np.cos(alpha + beta) + np.cos(3*(alpha + beta)) - 3*np.cos(5*(alpha + beta)) + 3*np.cos(5*alpha + beta) - 3*np.cos(alpha + 5*beta) - np.cos(3*alpha + 7*beta))*Power(sec(2*beta),5)*Power(sec(alpha + beta),3))/16.
 
  return (np.cos(2*alpha)*Power(sec(2*beta),5)*Power(sec(alpha + beta),3)*(4*np.cos(alpha - beta) + np.sin(2*beta)*(3*np.sin(3*(alpha + beta)) + np.sin(alpha + 5*beta))))/4.
  
def ddgeoFGamma2_2Ddfde_AntiAligned(e,f):

  alpha = np.arcsin(e)/2.0
  beta = np.arcsin(f)/2.0

  #return -0.125*((np.cos(alpha - 5*beta) - np.cos(alpha - beta) - 3*np.cos(3*alpha + beta) - np.cos(5*alpha + 3*beta))*sec(2*alpha)*Power(sec(2*beta),3)*Power(sec(alpha + beta),3))

  return -0.125*((np.cos(alpha - 5*beta) - np.cos(alpha - beta) - 3*np.cos(3*alpha + beta) - np.cos(5*alpha + 3*beta))*sec(2*alpha)*Power(sec(2*beta),3)*Power(sec(alpha + beta),3))


# class containing the eccentric mode equation in standard form for a shooting method
# self.forcing_frequency, self.csound2, self.dlnHcirc_dlna should be functions
# which specify the force precession frequency, circular counds speed squared and background
# pressure gradient (dlnHcirc_dlna = d ln (a P)/d ln a)
class EccentricDiscGamma2_2D( object ):

  def __init__(self):

    self.forcing_frequency = None # function pointer, function of (a,e)
    self.csound2 = None

    # function descibing the slope
    # of the internal energy in the disc
    # i.e. d ln(H_a^circ)/d ln a
    self.dlnHcirc_dlna = None
 
    self.GM = 1.0
 
    # add ability to set max e and f
    self.maxe = None
    self.maxf = None

    self._valid_solution=True

  # equation of motion, accepts precession frequency as an
  # argument
  def equation_of_motion(self,a,X,omega):

    e = X[0]
    ea = X[1]
    f = e+a*ea

    if not self._valid_solution:
      return [0.0,0.0]

    if self.maxe!=None:
      if np.abs(e)>self.maxe:
        self._valid_solution=False
        return [0.0,0.0]

    if self.maxf!=None:
      if np.abs(f)>self.maxe:
        self._valid_solution=False
        return [0.0,0.0]

    dFde = dgeoFGamma2_2Dde(e,f)

    # these are the same in 3d and 2d theory
    dFdf = dgeoFGamma2_2Ddf(e,f)
    ddFdf2 = ddgeoFGamma2_2Ddf2(e,f)
    ddFdfde = ddgeoFGamma2_2Ddfde(e,f)

    n = np.sqrt(self.GM/(a**3))
    
    omegaf = self.forcing_frequency(a,e)


    eaa = -(2.0/a)*ea + ((omega - omegaf)*n*a*a*e/(self.csound2(a)*np.sqrt(1 - e*e)) + dFde - a*ea*ddFdfde  - self.dlnHcirc_dlna(a)*dFdf)/(a*a*ddFdf2)

    return [ea,eaa]



# plot geometric part of the Hamiltonian as test
# probably good to test the 2D theory from the above
# against an existing implementation

if __name__ == '__main__':

  import matplotlib.pyplot as plt 

  elin = np.linspace(0.0,1.0,500)
  flin = np.linspace(-1.0,1.0,500)

  emesh, fmesh = np.meshgrid(elin,flin,indexing='ij')

  geoF = geoFGamma2_2D(emesh,fmesh)

  #plt.pcolormesh(emesh,fmesh,geoF)  

  plt.contourf(emesh,fmesh,geoF,levels=(-1.0,-0.5,-0.1,-0.01,0.0,0.01,0.1,0.5))
  plt.colorbar()


  plt.xlabel('e',fontsize=18)
  plt.ylabel('f',fontsize=18)

  plt.show()




