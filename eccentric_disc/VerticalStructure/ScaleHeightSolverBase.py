import numpy as np
from scipy.integrate import simps, solve_ivp
import scipy.optimize as opt
from ..geometry import EccentricGeometry


class ScaleHeightSolverBase( EccentricGeometry ):
 
  def __init__( self ):
    super(ScaleHeightSolverBase,self).__init__()
    #EccentricGeometry.__init__(self)    
 
    self.h0 = 1.0

  def __call__( self, x ):

   # note this requires x[0] = 0
   # will generalise later


    y0 = [self.h0,0.0]

    tspan=[x[0],x[-1]]

    sol = solve_ivp(self._eom,tspan,y0,t_eval=x)

    return sol.y

  def _eom(self,x,y):
    return None

  def solve(self,setSol=True):

    res = opt.newton(self.shoot_once,self.h0)

    # again putting this in to correct the issue with 
    # isothermal discs having a h, -h symmetry,
    res = np.abs(res)

    if setSol:
      self.h0 = res

    return res

  def shoot_once(self,h0):

    # putting this in to correct the issue with 
    # isothermal discs having a h, -h symmetry,

    h0 = np.abs(h0)
 
    y0 = [h0,0.0]

    # actually eccentric anomaly
    tspan = [0.0,np.pi]

    sol = solve_ivp(self._eom,tspan,y0)

    # possibly to way of avoiding negative h

    if np.min(sol.y[0])<0.0:
      return np.exp(-np.min(sol.y[0]))
    else:
      return sol.y[1,-1]



