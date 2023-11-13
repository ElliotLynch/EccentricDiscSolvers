import numpy as np
import scipy.special as spl
from .. import geometry as egeo
from .VerticalStructureBase import VerticalStructureBase
from .ScaleHeightSolverIdeal import ScaleHeightSolverIdeal


class VerticalStructureIdeal( VerticalStructureBase, ScaleHeightSolverIdeal ):
 
  def __inti__(self):
    super(VerticalStructureIdeal,self).__init__()

    #VerticalStructureBase.__init__(self)
    #ScaleHeightSolverIdeal.__init__(self)

  # going to use the Barker and Ogilvie/Ogilvie 01 terminology
  # for this

    # density vertical profile
  def Frho(self,X):

    polytropeindex = 1.0/(1.0 + self.gamma)

    ztld = self.stretched_z(X)

    return self.Cn(polytropeindex)*((1 - ztld*ztld/(2*polytropeindex  + 3.0))**polytropeindex)

   # pressure vertical profile
  def Fp(self,X):

    polytropeindex = 1.0/(1.0 + self.gamma)

    ztld = self.stretched_z(X)

    return ((2.0*polytropeindex+2.0)/(2.0*(polytropeindex + 1.0)))*self.Cn(polytropeindex)*((1 - ztld*ztld/(2*polytropeindex  + 3.0))**(polytropeindex+1))

  # double check what these should be
  def Ftemp(self,X):
    return self.Fp(X)/self.Frho(X)

  def Fsie(self,X):
    return self.Fp(X)/self.Frho(X)

   # normalisation constant
  def Cn(self,n): 
    return (spl.gamma(n+1/5)/spl.gamma(n+1))/np.sqrt((2.0*n+3.0)*np.pi)    







