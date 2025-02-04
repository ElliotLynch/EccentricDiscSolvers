import numpy as np
from . import EccentricDiscIso3D #as ied
from scipy.integrate import solve_ivp
import scipy.optimize as opt

## To begin with we want to solve for 
## the closed boundary case to check
## the method works

class EccentricDiscIso3D_AABinarySolver( EccentricDiscIso3D ):

  def __init__( self ):

   self.amax = 2.0 
   self.amin = 1.0
   
   # initial guess
   self.Amplitude0 = 0.1

   EccentricDiscIso3D.__init__(self) 

   self.inner_bc = None
   self.outer_bc = None

   # function for setting an amplitude dependant cavity size
   # resulting in a variable amin
   self.cavity = None


   # something large to ensure it's rejected
   self._fail = 1.0e6

   self.params=None

   self._N_A_COLPOINTS = 200

   #for bisection type methods
   self.search_width = 0.1
  
   # stores for the amplitde dependant
   # versions of these functions
   #self._csound2 = None
   #self._dlnHcirc_dlna = None
   #self._Ma = None

   self.acav=None

  #TODO: correct acav vs amin


  def __call__(self, a, Amplitude=None):

    if Amplitude==None:
      Amplitude=self.Amplitude0
 

    if self.cavity!=None:
      amin=self.cavity(Amplitude)
      self.acav=amin #note amin!=acav -> need to deal with that conversion internal to the function
    else:
      amin=self.amin
      self.acav=amin

    if self.params==None:
      yinit = self.inner_bc(self,amin,Amplitude) # need Amplitude as well
    else:
      yinit = self.inner_bc(self,amin,Amplitude,self.params)
 
    #aspan = [self.amin,self.amax]

    # integrate to existing position
   
    if np.isscalar(a):

      aspan = [amin,a]

      # weird this is the way that suggested to implent this
      sol = solve_ivp(lambda ax, y: self.equation_of_motion(ax,y,0.0),aspan,yinit) 

      return sol.y[:,-1]

    else:

      if amin!=a[0]:

        aspaninit = [amin,a[0]]

      # weird this is the way that suggested to implent this
        solinit = solve_ivp(lambda ax, y: self.equation_of_motion(ax,y,0.0),aspaninit,yinit) #,t_eval=a)

        y0 = solinit.y[:,-1]

      else:
        y0=yinit

      aspan = [a[0],a[-1]]

      sol = solve_ivp(lambda ax, y: self.equation_of_motion(ax,y,0.0),aspan,y0,t_eval=a)
 
      return sol.y


  def shoot_once(self,Amplitude):

    if self.cavity!=None:
      amin=self.cavity(Amplitude)
      self.acav=amin
    else:
      amin=self.amin
      self.acav=amin

    if self.params==None:
      y0 = self.inner_bc(self,amin,Amplitude)
    else:
      y0 = self.inner_bc(self,amin,Amplitude,self.params)
 
    aspan = [amin,self.amax]
 
      # why isn't there tolerances set?
      # looking for stationary solutions
    sol = solve_ivp(lambda a, y: self.equation_of_motion(a,y,0.0),aspan,y0)
    #sol = solve_ivp(lambda a, y: self.equation_of_motion(a,y,0.0),aspan,y0,rtol=1.0e-10,atol=1.0e-10)
    #sol = solve_ivp(lambda a, y: self.equation_of_motion(a,y,0.0),aspan,y0,method='LSODA')


    # if inner and outer bc are the same then
    # an unsuccesful integration can yield a success
    # as this will return the inner bc values
    if (sol.success==True) and (self._valid_solution):
      return self.outer_bc(self,self.amax,sol.y[:,-1])
    else:
      return self._fail

  # want to solve for given ea0 on inner bc?
  #def solve(self,setSol=True,emax=None,method='newton'):
  def solve(self,setSol=True,method='newton'):

    if self.cavity!=None:
      amin=self.cavity(self.Amplitude0)
      self.acav=amin
    else:
      amin=self.amin
      self.acav=amin


     #possibly check if scipy root (or similar) takes this as an argument 
    if method=='newton':
      #res = opt.newton(self.shoot_once,self.Amplitude0)  #,full_output=True)
      res = opt.newton(self.shoot_once,self.Amplitude0)
    elif method=='bisect':
      res = opt.bisect(self.shoot_once,self.Amplitude0*(1.0 - self.search_width),self.Amplitude0*(1.0 + self.search_width))


    if setSol:
      self.Amplitude0=res
        #self.ea0=ea0   

    return res





