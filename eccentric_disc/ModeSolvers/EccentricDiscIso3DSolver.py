import numpy as np
from . import EccentricDiscIso3D #as ied
from scipy.integrate import solve_ivp
import scipy.optimize as opt

## To begin with we want to solve for 
## the closed boundary case to check
## the method works

class EccentricDiscIso3DSolver( EccentricDiscIso3D ):

  def __init__( self ):

   self.amax = 2.0 
   self.amin = 1.0
   
   # initial guess
   self.omega0 = 0.1
   self.ea0 = 0.0

   EccentricDiscIso3D.__init__(self) 

   self.inner_bc = None
   self.outer_bc = None

   # something large to ensure it's rejected
   self._fail = 1.0e6

   self.params=None

   self._N_A_COLPOINTS = 200

   #for bisection type methods
   self.search_width = 0.1


  def __call__(self, a, omega=None):

    if omega==None:
      omega=self.omega0
 

    if self.params==None:
      yinit = self.inner_bc(self,self.amin)
    else:
      yinit = self.inner_bc(self,self.amin,self.params)
 
    #aspan = [self.amin,self.amax]

    # integrate to existing position
   
    if np.isscalar(a):

      aspan = [self.amin,a]

      # weird this is the way that suggested to implent this
      sol = solve_ivp(lambda ax, y: self.equation_of_motion(ax,y,omega),aspan,yinit) 

      return sol.y[:,-1]

    else:

      if self.amin!=a[0]:

        aspaninit = [self.amin,a[0]]

      # weird this is the way that suggested to implent this
        solinit = solve_ivp(lambda ax, y: self.equation_of_motion(ax,y,omega),aspaninit,yinit) #,t_eval=a)

        y0 = solinit.y[:,-1]

      else:
        y0=yinit

      aspan = [a[0],a[-1]]

      sol = solve_ivp(lambda ax, y: self.equation_of_motion(ax,y,omega),aspan,y0,t_eval=a)
 
      return sol.y


  def shoot_once(self,omega):

    if self.params==None:
      y0 = self.inner_bc(self,self.amin)
    else:
      y0 = self.inner_bc(self,self.amin,self.params)
 
    aspan = [self.amin,self.amax]
 
    sol = solve_ivp(lambda a, y: self.equation_of_motion(a,y,omega),aspan,y0)
    #sol = solve_ivp(lambda a, y: self.equation_of_motion(a,y,omega),aspan,y0,rtol=1.0e-10,atol=1.0e-10)
    #sol = solve_ivp(lambda a, y: self.equation_of_motion(a,y,omega),aspan,y0,method='LSODA')


    # if inner and outer bc are the same then
    # an unsuccesful integration can yield a success
    # as this will return the inner bc values
    if (sol.success==True) and (self._valid_solution):
      return self.outer_bc(self,self.amax,sol.y[:,-1])
    else:
      return self._fail

  # want to solve for given ea0 on inner bc?
  def solve(self,setSol=True,emax=None,method='newton'):

     # want to solve for the mode with a given emax
    if emax!=None:

      a_s=np.exp(np.linspace(np.log(self.amin),np.log(self.amax),self._N_A_COLPOINTS))

      omega00=self.omega0

      def _maxe_objective(inner_bc_param):

        self.params=inner_bc_param
        _res = opt.newton(self.shoot_once,omega00) # not clear which we want to do

        self.omega0=_res

         # does this assue set sol?
        e, ey = self(a_s)

        return np.max(e)-emax

        # not clear how to generalise - also is 0 singular?
      resEmax = opt.newton(_maxe_objective,0.0)

      # this doesn't make sense generally but trying to get this to work
      #resEmax = opt.bisect(_maxe_objective,0.0,1.0)

      self.params=resEmax

      res = opt.newton(self.shoot_once,omega00)  #,full_output=True)

      self.omega0=omega00 # then switch back if needed

    else:

      #need some form of switch
      
      #try: # bad practice need to do proper error catch
        #could try widening the search width until it crosses the root?    

      #res = opt.bisect(self.shoot_once,self.omega0*(1.0 - self.search_width),self.omega0*(1.0 + self.search_width))
      #except:
      #  res=self.omega0
      #  self.valid_solution=False

       #possibly check if scipy root (or similar) takes this as an argument 
      if method=='newton':
        res = opt.newton(self.shoot_once,self.omega0) #,full_output=True)
      elif method=='bisect':
        res = opt.bisect(self.shoot_once,self.omega0*(1.0 - self.search_width),self.omega0*(1.0 + self.search_width))

      # atempting with a different scheme
      #res_result = opt.root_scalar(self.shoot_once,x0=self.omega0,x1=0.0)
      #res=res_result.root

    if setSol:
      self.omega0=res
        #self.ea0=ea0   

    #print res

    return res



if __name__ == '__main__':

  import matplotlib.pyplot as plt

  cs0 = 0.1

  ea0=0.1

  def circular_inner(self,a): return [0.0,ea0]

  def circular_outer(self,a,y): return y[0] 

  def forcing_freq(a,e): return 0.0

  def csound2(a): return cs0*cs0

  def dlnHcirc_dlna(a): return 1.0

  solver=EccentricDiscIso3DSolver()
  solver.forcing_frequency = forcing_freq
  solver.csound2 = csound2
  solver.dlnHcirc_dlna = dlnHcirc_dlna

  solver.inner_bc = circular_inner
  solver.outer_bc = circular_outer

  solver.omega0=1.0*cs0*cs0

  #res = solver.shoot_once(solver.omega0,0.1)

  res = solver.solve(0.1)

  print(res)

  a_s=np.linspace(1.0,2.0,100)

  print(a_s)

  soly = solver(a_s)

  plt.plot(a_s,soly[0],'k-',linewidth=2)  

  solver.omega0=50.0*cs0*cs0

  res = solver.solve(0.1)
  soly = solver(a_s)

  plt.plot(a_s,soly[0],'k--',linewidth=2) 

  solver.omega0=1.0*cs0*cs0

  ea0=0.4

  solver.omega0=1.0*cs0*cs0

  res = solver.solve(0.1)
  soly = solver(a_s)

  plt.plot(a_s,soly[0],'k:',linewidth=2) 

  plt.show()

  e0 = 0.1

  def zero_gradient_inner(self,a): return [e0,0.0]

  def zero_gradient_outer(self,a,y): return y[1]

  solver.inner_bc = zero_gradient_inner
  solver.outer_bc = zero_gradient_outer

  solver.omega0=0.1*cs0*cs0

  res = solver.solve(0.1)
  soly = solver(a_s)

  #print 1

  plt.plot(a_s,soly[0],'k-',linewidth=2)  

  solver.omega0=5.0*solver.omega0  #0.2*cs0*cs0

  res = solver.solve(0.1)
  soly = solver(a_s)

  #print 2

  plt.plot(a_s,soly[0],'k--',linewidth=2)
  
  #solver.omega0=0.1*cs0*cs0

  e0=0.4
  #solver.omega0=1.0*cs0*cs0  

  res = solver.solve(0.1)
  soly = solver(a_s)

  #print 3

  plt.plot(a_s,soly[0],'k:',linewidth=2)

  plt.show()



