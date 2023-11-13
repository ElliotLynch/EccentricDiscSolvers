import numpy as np
from scipy.integrate import simps, solve_bvp
import scipy.optimize as opt
from .. import geometry as egeo
from .VerticalStructureBase import VerticalStructureBase
from .ScaleHeightSolverIsothermal import ScaleHeightSolverIsothermal


class VerticalStructureIsothermal( VerticalStructureBase, ScaleHeightSolverIsothermal ):
 
  def __init__(self):
    super(VerticalStructureIsothermal,self).__init__() 

    #VerticalStructureBase.__init__(self)
    #ScaleHeightSolverIsothermal.__init__(self)

  # going to use the Barker and Ogilvie/Ogilvie 01 terminology
  # for this

    # density vertical profile
  def Frho(self,X):

    ztld = self.stretched_z(X)

    return np.exp(-ztld*ztld/2.0)/np.sqrt(2.0*np.pi)

   # pressure vertical profile
  def Fp(self,X):

    ztld = self.stretched_z(X)

    return np.exp(-ztld*ztld/2.0)/np.sqrt(2.0*np.pi)

  def Ftemp(self,X):
    return np.ones_like(X[0])

  def Fsie(self,X):
    return np.ones_like(X[0])





