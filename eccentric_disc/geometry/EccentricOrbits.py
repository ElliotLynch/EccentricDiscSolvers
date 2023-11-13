import numpy as np
from . import EccentricTensorCalc
from . import coord_transforms as crdt
import eccentric_disc.datatools as dtools 

#egeom.EccentricGeometry()

class EccentricOrbits( object ):
  """
  class defining a set of orbits in an orbital coordinate system  
  
  Attributes
  ----------
  e : function
      eccentricity profile
  omega : function
      longitude of ascending node profile
  G : float
      Gravitational constant
  coord_system : str
      output coordinate system 
  dy : float
      y step for computing derivatives
  solver : class
      Mode solver for computing orbits
  
  Methods
  -------
  meanmotion(a)
  ey(a)
  e_wy(a)
  q(a)
  alpha(a)
  radius(a,Ecc)
  angular_vel(a,Ecc)
  jacobian(a,Ecc)
  jacobian_aE(a,Ecc)
  aE_2_lamphi(a,Ecc)
  X(x)
  V(x)
  save(s_grid,save_grid,**kwargs)
  
  """

  # we want to use EccentricGeometry to deal with conversions between the orbital
  # coordinate systems
  # for now just impliment the a E coorinate system implement others when needed

  def __init__( self ):

    # function for eccentricity
    self._e = None

    # function for longitude of pericentre 
    # check what symbol/term I want to use here
    # probably make consistent with REBOUND 
    self._omega = None

    self.coord_system = 'CART'

    # specify an eccentric disc solver to 
    # yeild the eccentric orbits
    # ideally want an easy way of just setting orbital elements
    self._solver = None

    self.auto_solve = True

    ####
    # 'CART' : cartesian coordinates
    # 'CYL' : cylindrical coordinates
    # 'AE' : a E orbital coordinate system - corresponds to the identity map
    # 'LPHI' : NOT IMPLIMENTED YET
    ####

    # assume grid of coordinates is a regular 
    # grid of orbital coordinates
    self.orbital_regular_grid = True

    # use EccentricTensorCalc as an underlying converter
    self._eccentrictensor_inst = EccentricTensorCalc()
 
    self.dy = 1.0e-6

    self.G=1

    # set ey to default function
    self.ey = self._ey_default


  def meanmotion(self,a):
    """Mean motion
  
    Parameters
    ---------- 
    a : float
        semimajor axis 
 
    """
   
    return self.G*(a**(-1.5))

  @property
  def e(self):
    return self._e 

  @e.setter
  def e(self, e):
    self._e = e

  #@property
  #def omega(self,a):
  #  return self._omega(a)

  @property
  def omega(self):
    return self._omega
 
  @omega.setter
  def omega(self, omega):
    self._omega = omega

  @property
  def ey(self):
    return self._ey

  @ey.setter
  def ey(self, ey):
    self._ey = ey

  def _ey_default(self, a):
    """Calculates a de/da, do this if function isn't specified
  
    Parameters
    ---------- 
    a : float
        semimajor axis 
 
    """
    return (self._e(a+a*self.dy)-self._e(a-a*self.dy))/(2.0*self.dy)

  def e_wy(self, a):
    """Calculates a e d (longitude of pericentre)/da
  
    Parameters
    ---------- 
    a : float
        semimajor axis 
 
    """

     # not sure which of these we 
     # need   
    if type(self._omega) == type(None):
      return np.zeros_like(a)
    elif type(self._omega(a+a*self.dy)) == type(None):
      return np.zeros_like(a)
    else:
      return self._e(a)*(self._omega(a+a*self.dy)-self._omega(a-a*self.dy))/(2.0*self.dy)

  # ok let's impliment all these using eccentricgeometry
  # as for the underlying calculations

  def q(self, a):
    """Calculates orbital intersection/nonlinearity parameter
  
    Parameters
    ---------- 
    a : float
        semimajor axis 
 
    """
   
    # need to do it this way to get around
    # issues with initialising arrays of eccentricity
    # and calling _recalculate_derivatives

    self._eccentrictensor_inst._e = self.e(a)
    self._eccentrictensor_inst._ey = self.ey(a)
    self._eccentrictensor_inst._e_wy = self.e_wy(a)
  
    self._eccentrictensor_inst._recalculate_qalpha()

    return self._eccentrictensor_inst.q

  def alpha(self, a):
    """Calculates alpha, the angle determining the relative contribution of twist and eccentricity gradients to q
  
    Parameters
    ---------- 
    a : float
        semimajor axis 
 
    """
   
    self._eccentrictensor_inst._e = self.e(a)
    self._eccentrictensor_inst._ey = self.ey(a)
    self._eccentrictensor_inst._e_wy = self.e_wy(a)

    self._eccentrictensor_inst._recalculate_qalpha()

    return self._eccentrictensor_inst.alpha

  def radius(self,a,Ecc):
   """Calculates the cylindrical radius from the central object
  
    Parameters
    ---------- 
    a : float
        semimajor axis 
    Ecc : float
        Eccentric Anomaly
 
    """
   
   self._eccentrictensor_inst._e = self.e(a)
   self._eccentrictensor_inst._ey = self.ey(a)
   self._eccentrictensor_inst._e_wy = self.e_wy(a)  
 
   self._eccentrictensor_inst._recalculate_qalpha()

   return a*self._eccentrictensor_inst.radius(Ecc)

  def angular_vel(self,a,Ecc):
    """Calculates the angular velocity
  
    Parameters
    ---------- 
    a : float
        semimajor axis
    Ecc : float
        Eccentric Anomaly
 
    """

    self._eccentrictensor_inst._e = self.e(a)
    self._eccentrictensor_inst._ey = self.ey(a)
    self._eccentrictensor_inst._e_wy = self.e_wy(a)

    self._eccentrictensor_inst._recalculate_qalpha()

    # need to scale correctly

    return self.meanmotion(a)*self._eccentrictensor_inst.angular_vel(Ecc)

   # note this is the Jacobian of the a M coordinate system not
   # the a E coordinate system that we've been using through most
   # of this
  def jacobian(self,a,Ecc):
    """Calculates the jacobian of the (Lambda,lambda) cannonical coordinate system
  
    Parameters
    ---------- 
    a : float
        semimajor axis
    Ecc : float
        Eccentric Anomaly
 
    """

    self._eccentrictensor_inst._e = self.e(a)
    self._eccentrictensor_inst._ey = self.ey(a)
    self._eccentrictensor_inst._e_wy = self.e_wy(a)

    self._eccentrictensor_inst._recalculate_qalpha()

    return a*self._eccentrictensor_inst.jacobian(Ecc)


  def jacobian_aE(self,a,Ecc):
    """Calculates the jacobian of the (a,E) cannonical coordinate system
  
    Parameters
    ---------- 
    a : float
        semimajor axis
    Ecc : float
        Eccentric Anomaly
 
    """

    self._eccentrictensor_inst._e = self.e(a)
    self._eccentrictensor_inst._ey = self.ey(a)
    self._eccentrictensor_inst._e_wy = self.e_wy(a)

    self._eccentrictensor_inst._recalculate_qalpha()

    return a*(1.0-self.e(a)*np.cos(Ecc))*self._eccentrictensor_inst.jacobian(self,Ecc)

  def horizontal_divergence(self,a,Ecc):
    """Calculates the horizontal velocity divergence
  
    Parameters
    ---------- 
    a : float
        semimajor axis
    Ecc : float
        Eccentric Anomaly
 
    """
 
    self._eccentrictensor_inst._e = self.e(a)
    self._eccentrictensor_inst._ey = self.ey(a)
    self._eccentrictensor_inst._e_wy = self.e_wy(a)

    self._eccentrictensor_inst._recalculate_qalpha()

    return self.meanmotion(a)*self._eccentrictensor_inst.horizontal_divergence(Ecc)

  def aE_2_lamphi(self,a,Ecc):
    """Conversion of (a,E) to (lam,phi) coordinate system
  
    Parameters
    ---------- 
    a : float
        semimajor axis
    Ecc : float
        Eccentric Anomaly
 
    """

    self._eccentrictensor_inst._e = self.e(a)
    self._eccentrictensor_inst._ey = self.ey(a)
    self._eccentrictensor_inst._e_wy = self.e_wy(a)

    self._eccentrictensor_inst._recalculate_qalpha()

    return self._eccentrictensor_inst.aE_2_lamphi(Ecc)

  @property
  def solver(self):
    return self._solver

  # require user to call solve?

  @solver.setter
  def solver(self,solver):
    self._solver = solver

    self._e = self._e_from_solver
    self._ey = self._ey_from_solver
    self._omega = self._omega_from_solver

    # not sure where this is supposed to inherit from
    #self.G=self._solver.G

    if self.auto_solve:
      self._solver.solve()

  def _e_from_solver(self,x):
    # assume that x is a grid of semimajor axis

    # can use a faster method if the grid is 
    # regularly spaced in orbital coords
    if self.orbital_regular_grid:
       
      if len(x.shape)==3:
        semimajor = x[:,0,0]

        #print semimajor

        # the reason this version doesn't work is it changes the location of 
        # the boundary - but surely scalor should stop that?
        #eccentricity = [self._solver(xval)[0] for xval in semimajor]



        eccentricity = self._solver(semimajor)[0]

        #print 'e = ' , eccentricity
        
        #raise

        eccentricity, dummy1, dummy2 = np.meshgrid(eccentricity,x[0,:,0],x[0,0,:],indexing='ij')

        #print eccentricity[:,0,0]

        #print 'n=3'

        #raise

        return eccentricity

      #2D
      elif len(x.shape)==2:
        semimajor = x[:,0]

        
        #eccentricity = [self._solver(xval)[0] for xval in semimajor]

        # need to go back and check this against the examples
        eccentricity = self._solver(semimajor)[0]

        #print 'n=2'

        # copy along eccentric anomaly direction using numpy outer
        return np.outer(eccentricity,np.ones_like(x[0]))
     
      else:
        semimajor = x[:]
      
        #print 'n=1'

        #try this?
        #eccentricity = [self._solver(xval)[0] for xval in semimajor]

        # this should be faster - I don't think there any
        # disadvantage to doing it this way assuming the correct ordering
        eccentricity = self._solver(semimajor)[0]

        return np.array(eccentricity)


    else:
      # also assume arraylike
      shape = x.shape

      # this will be slow but do it this way for now 
      # and write something more efficient if needed

      # still need to loop over this? - but not this way
      # not sure if this is the correct way of looping over this
      res = [self._solver(xval)[0] for xval in x.flatten()]

      return np.reshape(np.array(res),shape)



  def _ey_from_solver(self,x):
    # assume that x is a grid of semimajor axis

    # can use a faster method if the grid is 
    # regularly spaced in orbital coords
    if self.orbital_regular_grid:

      if len(x.shape)==3:
        semimajor = x[:,0,0]
        #print semimajor

        #eccentricity = [self._solver(xval)[0] for xval in semimajor]
        ey = semimajor*self._solver(semimajor)[1]

        #print 'ey = ' , ey


        ey, dummy1, dummy2 = np.meshgrid(ey,x[0,:,0],x[0,0,:],indexing='ij')

        #print ey[:,0,0]

        #print 'n=3'

        #raise

        return ey

      #2D
      elif len(x.shape)==2:
        semimajor = x[:,0]

        # this is probably wrong
        #ey = [self._solver(xval)[1] for xval in semimajor]

        ey = semimajor*self._solver(semimajor)[1]

        #print 'n=2'

        # copy along eccentric anomaly direction using numpy outer
        return np.outer(ey,np.ones_like(x[0]))

      else:
        semimajor = x[:]

        #print 'n=1'

        #try this?
        #ey = [self._solver(xval)[1] for xval in semimajor]

        ey = semimajor*self._solver(semimajor)[1]

        return np.array(ey)


    else:
      # also assume arraylike
      shape = x.shape

      # this will be slow but do it this way for now 
      # and write something more efficient if needed

      # still need to loop over this
      # again this isn't correct

      res = [xval*self._solver(xval)[1] for xval in x.flatten()]

      return np.reshape(np.array(res),shape)

  def _omega_from_solver(self,x):
    
    # not clear how to obtain this from solver
    # notably for modes this isn't present as it's trivially zero

    pass

  def X(self,x):
    """Calculates position in coordinate system specified by coord_system attribute
  
    Parameters
    ---------- 
    x : (array, array)
        grid of semimajor axis and grid of eccentric anomaly
      
    """

    #not sure if need to play around with _recalculate_qalpha()


    if self.coord_system == 'CART':
      return crdt.aE2cart(self,x)
    elif self.coord_system == 'CYL':    
      return crdt.aE2cyl(self,x)
    elif self.coord_system == 'AE':
      return x

    # probably should return an error somewhere if a coordinate system
    # is chosen that's not implemented (with the list of implimented systems)


    # This yeilds the covarient velocity for different coordinate systems 
    # ideally probably want to centralise these functions and import them
    # 
  def V(self,x):
    """Calculates velocity in coordinate system specified by coord_system attribute
  
    Parameters
    ---------- 
    x : (array, array)
        grid of semimajor axis and grid of eccentric anomaly
      
    """

    # for now leave off implimenting vz
    # return vz=0 if 3d

    self._eccentrictensor_inst._e = self.e(x[0])
    self._eccentrictensor_inst._ey = self.ey(x[0])
    self._eccentrictensor_inst._e_wy = self.e_wy(x[0])

    self._eccentrictensor_inst._recalculate_qalpha()

    #angular_velo = self.meanmotion(x[0])*self._eccentrictensor_inst.angular_vel(x[1]) 

    meann = self.meanmotion(x[0])

    if self.coord_system == 'CART':
      eccentricity = self.e(x[0])

      vr = x[0]*meann*eccentricity*np.sin(x[1])/(1.0 - eccentricity*np.cos(x[1]))
      vphi = np.sqrt(1 - eccentricity*eccentricity)*meann*((1.0 - eccentricity*np.cos(x[1]))**(-2))

      if len(x)==3:
        r, phi, z = crdt.aE2cyl(self,x)

        vx = vr*np.cos(x[1]) - r*np.sin(x[1])*vphi
        vy = vr*np.sin(x[1]) + r*np.cos(x[1])*vphi
      
        return vx, vy, np.zeros_like(vx)

      else:
 
        r, phi = crdt.aE2cyl(self,x)

        vx = vr*np.cos(x[1]) - r*np.sin(x[1])*vphi
        vy = vr*np.sin(x[1]) + r*np.cos(x[1])*vphi
      
        return vx, vy

    elif self.coord_system == 'CYL':
      eccentricity = self.e(x[0])

      vr = x[0]*meann*eccentricity*np.sin(x[1])/(1.0 - eccentricity*np.cos(x[1]))
      vphi = np.sqrt(1 - eccentricity*eccentricity)*meann*((1.0 - eccentricity*np.cos(x[1]))**(-2))

      if len(x)==3:
        return vr, vphi, np.zeros_like(vr)
      else:
        return vr, vphi

    elif self.coord_system == 'AE':
      
      va = 0.0
      vE = meann/(1.0 - self.e(x[0])*np.cos(x[1]))

      if len(x)==3:
        return va, vE, np.zeros_like(va)
      else:
        return va, vE

      # probably should return an error somewhere if a coordinate system
      # is chosen that's not implemented (with the list of implimented systems)


  def transform_covariant(self,V,X,coord_system = 'AE'):

    self._eccentrictensor_inst._e = self.e(X[0])
    self._eccentrictensor_inst._ey = self.ey(X[0])
    self._eccentrictensor_inst._e_wy = self.e_wy(X[0])

    self._eccentrictensor_inst._recalculate_qalpha()

    VNEW = self._eccentrictensor_inst.transform_covariant(V,X,coord_system = coord_system,omega=self._omega(X[0])) # pressumably the treatment of omega works

    return VNEW

  def transform_contravarient(self,V,X,coord_system = 'AE'):

    self._eccentrictensor_inst._e = self.e(X[0])
    self._eccentrictensor_inst._ey = self.ey(X[0])
    self._eccentrictensor_inst._e_wy = self.e_wy(X[0])

    self._eccentrictensor_inst._recalculate_qalpha()

    VNEW = self._eccentrictensor_inst.transform_contravarient(V,X,coord_system = coord_system,omega=self._omega(X[0])) # pressumably the treatment of omega works

    return VNEW

  #def save(self,s_grid,save_grid,*kwargssave_format='PLUTO',*args):
  def save(self,s_grid,save_grid,**kwargs):
    """Save to a specified data format
  
    Parameters
    ---------- 
    s_grid : array
        grid of coordinates points to compute solutions on
    save_grid : array
        grid to coordinates points to save to
 
    """

    if 'save_format' not in kwargs:
      kwargs['save_format'] = 'PLUTO'


    # should we really have exceptions here as save_format PLUTO doesn't 
    # make sense for the class EccentricOrbits
    # we'll probably do it that way once we add more than one save_format!
    return dtools.save(self,s_grid,save_grid,**kwargs) 




