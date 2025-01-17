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


  # we must have to modify the inner boundary as that's where the mode amplitude is set
  # and we have replaced omega with the mode amplitude as a parameter

  #@property
  #def csound2(self,a):
  #  if self.cavity == None:
  #    return self._csound2(a)
  #  else:
  #    self.cavity = self._csound2(a,self.cavity(self.Amplitude0)) # will this give the correct behaviour?

  #@csound2.setter
  #def csound2(self, csound2):
  #  self._csound2 = csound2
  
  #@property
  #def dlnHcirc_dlna(self,a):
  #  if self.cavity == None:
  #    return self._dlnHcirc_dlna(a)
  #  else:
  #    self.dlnHcirc_dlna = self._dlnHcirc_dlna(a,self.cavity(self.Amplitude0)) # will this give the correct behaviour?

  #@dlnHcirc_dlna.setter
  #def dlnHcirc_dlna(self, dlnHcirc_dlna):
  #  self._dlnHcirc_dlna = dlnHcirc_dlna

  #@property
  #def Ma(self,a):
  #  if self.Ma == None:
  #    return self._Ma(a)
  #  else:
  #    self.Ma = self._Ma(a,self.cavity(self.Amplitude0)) # will this give the correct behaviour?

  #@Ma.setter
  #def Ma(self, Ma):
  #  self._Ma = Ma

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

     # want to solve for the mode with a given emax
     # does this makes sense for this solution - I suspect not as there
     # can only be one
    #if emax!=None:

    #  a_s=np.exp(np.linspace(np.log(amin),np.log(self.amax),self._N_A_COLPOINTS))

    #  Amplitude00=self.Amplitude0

    #  def _maxe_objective(inner_bc_param):

    #    self.params=inner_bc_param
    #    _res = opt.newton(self.shoot_once,Amplitude00) # not clear which we want to do

    #    self.Amplitude0=_res

         # does this assue set sol?
    #    e, ey = self(a_s)

    #    return np.max(e)-emax

        # not clear how to generalise - also is 0 singular?
    #  resEmax = opt.newton(_maxe_objective,0.0)

      # this doesn't make sense generally but trying to get this to work
      #resEmax = opt.bisect(_maxe_objective,0.0,1.0)

    #  self.params=resEmax

    #  res = opt.newton(self.shoot_once,Amplitude00)  #,full_output=True)

    #  self.Amplitude0=Amplitude00 # then switch back if needed

    #else:

      #need some form of switch
      
      #try: # bad practice need to do proper error catch
        #could try widening the search width until it crosses the root?    

      #res = opt.bisect(self.shoot_once,self.omega0*(1.0 - self.search_width),self.omega0*(1.0 + self.search_width))
      #except:
      #  res=self.omega0
      #  self.valid_solution=False

     #possibly check if scipy root (or similar) takes this as an argument 
    if method=='newton':
      res = opt.newton(self.shoot_once,self.Amplitude0)  #,full_output=True)
    elif method=='bisect':
      res = opt.bisect(self.shoot_once,self.Amplitude0*(1.0 - self.search_width),self.Amplitude0*(1.0 + self.search_width))

    # atempting with a different scheme
    #res_result = opt.root_scalar(self.shoot_once,x0=self.Amplitude0,x1=0.0)
    #res=res_result.root

    if setSol:
      self.Amplitude0=res
        #self.ea0=ea0   

    #print res

    return res



#if __name__ == '__main__':




