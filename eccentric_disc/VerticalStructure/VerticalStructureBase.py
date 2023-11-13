import numpy as np
from scipy.integrate import simps, solve_bvp
import scipy.optimize as opt
from .. import geometry as egeo



class VerticalStructureBase( object ):
 
  def __init__(self):
    super(VerticalStructureBase,self).__init__()

    self.coord_system = 'CART'

    # probably best to do it this way
    self.orbits = None
  
     # circular referance scale height
    self.scale_heightcirc = None                                      

    # for now can only do it this way
    self.orbital_regular_grid = True

  def _to_aE(self,X):
    if self.coord_system == 'CART':
      a, eccanom, z = egeo.coord_transforms.cart2aE(self.orbits,X)
    elif self.coord_system == 'CYL':
      a, eccanom, z = egeo.coord_transforms.cyl2aE(self.orbits,X)
    elif self.coord_system == 'AE':
      a, eccanom, z = X
    elif self.coord_system == 'APHI':
      a, eccanom, z = egeo.coord_transforms.aphi2aE(self.orbits,X,increasingEccanom=True)

    return a, eccanom, z

   # dimensionless scale height
  def get_h(self,X):

    # for testing until I think of a better way to do this
    #self.h0 = 1.0
    # we want to make sure all these functions preserve h0
    h00=self.h0

    a, eccanom, z = self._to_aE(X)

    h_horizontal = np.zeros_like(a[:,:,0])
    
    h = np.zeros_like(a)

    ecc = self.orbits.e(a[:,:,0])
    nonlinearity_q = self.orbits.q(a[:,:,0])
    alpha = self.orbits.alpha(a[:,:,0])

    if self.orbital_regular_grid:
      for i in range(a.shape[0]):
        self.e = ecc[i,0]
        self.q = nonlinearity_q[i,0]

        # strictly this has to be zero for
        # the solver to work
        self.alpha = alpha[i,0]

        self.solve()

        #print self.h0


        #print a[i,0,0]
        #print eccanom[i,:,0]

        h_horizontal[i,:] = self(eccanom[i,:,0])[0,:]
        #print h_horizontal[i,:]


    else:
      raise Error('orbit iregular grid not implimented')

     # being lazy
    for i in range(a.shape[2]):
      h[:,:,i] = h_horizontal

    # reset h0
    self.h0=h00

    return h

   # dimensionless scale height
  def get_dh(self,X):

    h00=self.h0

    a, eccanom, z = self._to_aE(X)

    dh_horizontal = np.zeros_like(a[:,:,0])

    dh = np.zeros_like(a)

    ecc = self.orbits.e(a[:,:,0])
    nonlinearity_q = self.orbits.q(a[:,:,0])
    alpha = self.orbits.alpha(a[:,:,0])

    if self.orbital_regular_grid:
      for i in range(a.shape[0]):
        self.e = ecc[i,0]
        self.q = nonlinearity_q[i,0]

        # strictly this has to be zero for
        # the solver to work
        self.alpha = alpha[i,0]

        self.solve()

        dh_horizontal[i,:] = self(eccanom[0,:,0])[1,:]

    else:
      raise Error('orbit iregular grid not implimented')

     # being lazy
    for i in range(a.shape[2]):
      dh[:,:,i] = dh_horizontal

    # reset h0
    self.h0=h00

    return dh

  # dimensional scale height
  # use function above to avoid code duplication
  def get_dH(self,X):

    h00=self.h0

    a, eccanom, z = self._to_aE(X)

    dh_horizontal = np.zeros_like(a[:,:,0])

    dh = np.zeros_like(a)

    ecc = self.orbits.e(a[:,:,0])
    nonlinearity_q = self.orbits.q(a[:,:,0])
    alpha = self.orbits.alpha(a[:,:,0])

    if self.orbital_regular_grid:
      for i in range(a.shape[0]):
        self.e = ecc[i,0]
        self.q = nonlinearity_q[i,0]

        # strictly this has to be zero for
        # the solver to work
        self.alpha = alpha[i,0]

        self.solve()

        dh_horizontal[i,:] = self(eccanom[0,:,0])[1,:]

    else:
      raise Error('orbit iregular grid not implimented')

     # being lazy
    for i in range(a.shape[2]):
      dh[:,:,i] = dh_horizontal

    dH = self.scale_heightcirc(a)*dh

    # reset h0
    self.h0=h00

    return dH


   # dimensional scale height
  def get_H(self,X):
    
    # for testing
    h00=self.h0

    a, eccanom, z = self._to_aE(X)

    h_horizontal = np.zeros_like(a[:,:,0])

    h = np.zeros_like(a)
    
    ecc = self.orbits.e(a[:,:,0])
    nonlinearity_q = self.orbits.q(a[:,:,0])
    alpha = self.orbits.alpha(a[:,:,0])

    if self.orbital_regular_grid:
      for i in range(a.shape[0]):
        self.e = ecc[i,0]
        self.q = nonlinearity_q[i,0]

        # strictly this has to be zero for
        # the solver to work
        self.alpha = alpha[i,0]

        self.solve()


        #print self.h0

        h_horizontal[i,:] = self(eccanom[0,:,0])[0,:]

    else:
      raise Error('orbit iregular grid not implimented')

    for i in range(a.shape[2]):
      h[:,:,i] = h_horizontal

    H = self.scale_heightcirc(a)*h
 
    # reset h0
    self.h0=h00

    return H

   # stretched vertical coordinate
   # not sure if this should be somewhere else?
   # or possibly just the base class
  def stretched_z(self, X):

    h00=self.h0

    a, eccanom, z = self._to_aE(X) 

    h_horizontal = np.zeros_like(a[:,:,0])
 
    h = np.zeros_like(a)

    ecc = self.orbits.e(a[:,:,0])
    nonlinearity_q = self.orbits.q(a[:,:,0])
    alpha = self.orbits.alpha(a[:,:,0])

    if self.orbital_regular_grid: 
      for i in range(a.shape[0]):
        self.e = ecc[i,0]
        self.q = nonlinearity_q[i,0]

        # strictly this has to be zero for
        # the solver to work
        self.alpha = alpha[i,0]

        self.solve()

        h_horizontal[i,:] = self(eccanom[0,:,0])[0,:]
 
    else:
      raise Error('orbit iregular grid not implimented')

    for i in range(a.shape[2]):
      h[:,:,i] = h_horizontal

    H = self.scale_heightcirc(a)*h

    # reset h0
    self.h0=h00

    return z/H


  # going to use the Barker and Ogilvie/Ogilvie 01 terminology
  # for this

    # density vertical profile
  def Frho(self,X):

    return np.ones_like(X[0])

   # pressure vertical profile
  def Fp(self,X):

    return np.ones_like(X[0])







