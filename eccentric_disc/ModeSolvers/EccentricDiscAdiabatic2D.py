import numpy as np
from . import SpecialIp, SpecialIp1, DSpecialIp, DSpecialIp1, DDSpecialIp, DDSpecialIp1

# defining these so we can make use
# of CForm
def sec(x): return 1.0/np.cos(x)

def csc(x): return 1.0/np.sin(x)

def Power(x,n): return x**n

 #check gamma vs gamma-1 throughout

 # geometric part of the Hamiltonian
def geoF2D(e,f,gamma):

  return ((Power((1 - e*f)/np.sqrt(1 - Power(e,2)),1 - gamma)*
     (SpecialIp((e - f)/(-1 + e*f),gamma-1.0) - e*SpecialIp1((e - f)/(-1 + e*f),gamma-1.0)))/(-1 + gamma))

# various derivatives of the geometric part of the Hamiltonian

def dgeoF2Dde(e,f,gamma):

  return ((Power(1 - Power(e,2),(-3 + gamma)/2.)*
     ((e - f)*(-1 + gamma)*(-SpecialIp((e - f)/(-1 + e*f),gamma-1.0) + e*SpecialIp1((e - f)/(-1 + e*f),gamma-1.0)) - 
       ((-1 + Power(e,2))*(-((-1 + Power(f,2))*DSpecialIp((e - f)/(-1 + e*f),gamma-1.0)) + 
            e*(-1 + Power(f,2))*DSpecialIp1((e - f)/(-1 + e*f),gamma-1.0) + 
            Power(-1 + e*f,2)*SpecialIp1((e - f)/(-1 + e*f),gamma-1.0)))/(-1 + e*f)))/(Power(1 - e*f,gamma)*(-1 + gamma)))

def dgeoF2Ddf(e,f,gamma):

  return -((Power(1 - Power(e,2),(-1 + gamma)/2.)*Power(1 - e*f,-1 - gamma)*
       ((-1 + Power(e,2))*DSpecialIp((e - f)/(-1 + e*f),gamma-1.0) - 
         e*(-1 + Power(e,2))*DSpecialIp1((e - f)/(-1 + e*f),gamma-1.0) + 
         e*(-1 + e*f)*(-1 + gamma)*(SpecialIp((e - f)/(-1 + e*f),gamma-1.0) - e*SpecialIp1((e - f)/(-1 + e*f),gamma-1.0))))/(-1 + gamma))

def ddgeoF2Ddf2(e,f,gamma):

  return -((Power(1 - Power(e,2),(-1 + gamma)/2.)*Power(1 - e*f,-3 - gamma)*
       (-(Power(-1 + Power(e,2),2)*DDSpecialIp((e - f)/(-1 + e*f),gamma-1.0)) + 
         e*Power(-1 + Power(e,2),2)*DDSpecialIp1((e - f)/(-1 + e*f),gamma-1.0) + 
         e*(-1 + e*f)*gamma*(-2*(-1 + Power(e,2))*DSpecialIp((e - f)/(-1 + e*f),gamma-1.0) + 
            2*e*(-1 + Power(e,2))*DSpecialIp1((e - f)/(-1 + e*f),gamma-1.0) + 
            e*(-1 + e*f)*(-1 + gamma)*(-SpecialIp((e - f)/(-1 + e*f),gamma-1.0) + e*SpecialIp1((e - f)/(-1 + e*f),gamma-1.0)))))/(-1 + gamma))


def ddgeoF2Ddfde(e,f,gamma):

  return ((Power(1 - Power(e,2),(-3 + gamma)/2.)*Power(1 - e*f,-3 - gamma)*
     (Power(-1 + Power(e,2),2)*(-1 + Power(f,2))*DDSpecialIp((e - f)/(-1 + e*f),gamma-1.0) - 
       e*Power(-1 + Power(e,2),2)*(-1 + Power(f,2))*DDSpecialIp1((e - f)/(-1 + e*f),gamma-1.0) + 
       (-1 + e*f)*((-1 + Power(e,2))*
           (e*Power(f,2)*(-1 + gamma) - 2*e*gamma + f*(1 + gamma))*
           DSpecialIp((e - f)/(-1 + e*f),gamma-1.0) - 
          (-1 + Power(e,2))*(1 + Power(e,3)*f + 
             Power(e,2)*(-1 + Power(f,2)*(-1 + gamma) - 2*gamma) + e*f*gamma)*
           DSpecialIp1((e - f)/(-1 + e*f),gamma-1.0) - 
          (-1 + e*f)*(-1 + gamma)*((-1 + Power(e,2)*gamma + e*(f - f*gamma))*
              SpecialIp((e - f)/(-1 + e*f),gamma-1.0) + 
             e*(2 + Power(e,3)*f + e*f*(-2 + gamma) - Power(e,2)*(1 + gamma))*SpecialIp1((e - f)/(-1 + e*f),gamma-1.0)))))/(-1 + gamma))



# class containing the eccentric mode equation in standard form for a shooting method
# self.forcing_frequency, self.csound2, self.dlnHcirc_dlna should be functions
# which specify the force precession frequency, circular counds speed squared and background
# pressure gradient (dlnHcirc_dlna = d ln (a P)/d ln a)
class EccentricDiscAdiabatic2D( object ):

  def __init__(self):

    self.forcing_frequency = None # function pointer, function of (a,e)
    self.csound2 = None

    # function descibing the slope
    # of the internal energy in the disc
    # i.e. d ln(H_a^circ)/d ln a
    self.dlnHcirc_dlna = None
 
    self.GM = 1.0

    self.is3D = True

    self.gamma=2.0

    # add ability to set max e and f
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

    dFde = dgeoF2Dde(e,f,self.gamma)

    # these are the same in 3d and 2d theory
    dFdf = dgeoF2Ddf(e,f,self.gamma)
    ddFdf2 = ddgeoF2Ddf2(e,f,self.gamma)
    ddFdfde = ddgeoF2Ddfde(e,f,self.gamma)

    n = np.sqrt(self.GM/(a**3))
    
    omegaf = self.forcing_frequency(a,e)

     #check there aren't factors of gamma or similar missing?
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

  geoF = geoF2D(emesh,fmesh,5.0/3.0)

  #plt.pcolormesh(emesh,fmesh,geoF)  

  plt.contourf(emesh,fmesh,geoF,levels=(-1.0,-0.5,-0.1,-0.01,0.0,0.01,0.1,0.5))
  plt.colorbar()


  plt.xlabel('e',fontsize=18)
  plt.ylabel('f',fontsize=18)

  plt.show()






