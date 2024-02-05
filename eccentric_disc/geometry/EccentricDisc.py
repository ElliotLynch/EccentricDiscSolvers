import numpy as np
from . import EccentricGeometry
from . import EccentricOrbits
import scipy.integrate as intg

# not sure if I want this in geometry
# it's sort of a data structure

#egeom.EccentricGeometry()


 # class that allows for calculation of various fluid properties
 # in eccentric disc, on top of the disc orbital geometry inherited
 # from eccentric orbits.
class EccentricDisc( EccentricOrbits ):

  # for now we are ignoring vertical structure until 
  # we a ready to impliment 3D stuff

  def __init__( self ):
    EccentricOrbits.__init__(self)
 
    self.eos = None

    # ideally want to connect these to 
    # what's used to setup the solver
    # for now do manually
    self.Ma = None
 
    self.is3D = True
 
    # use entropy as this is a Lagrangian
    # scalor in ideal theory
    self.entropy = None

     # map between property names
     # for __getitem__
     # and functions defining those properties
     # will mostly make use of PLUTOs naming scheme
    self._propery_map = {"RHO" : self.den, "PRS" : self.pre, "DEN" : self.den, "H" : self.scale_height, "dH" : self.dH, "SDEN" : self.sden, "TMP" : self.temp, "TEMP" : self.temp, "V" : self.V, "X" : self.X , "SIE" : self.internal_energy, "B" : self.B}

    self.vertical_structure = None


  def __getitem__(self,key):

    # ok for now forget about implimenting each coordinate/velocity component seperately
 
    if key in self._propery_map:
      return self._propery_map[key]

    raise KeyError(key)
 
  def den(self,X):
    a = X[0]
    eccanom = X[1]

    if self.is3D:
      h = self.vertical_structure.get_H(X)
    else:
      h = 1.0 # for now ignoring the vertical structure

    #print self.Ma(a)
    #print self.vertical_structure.Frho(X)
    #print self.jacobian(a,eccanom)
    #print 'h, ', h

    # the jacobian is the actual jacobian, not the dimensionless one
    # this is the jacobian inherited from EccentricOrbits - which is dimensionfull
 
    if self.is3D:
      return self.Ma(a)*self.vertical_structure.Frho(X)/(2.0*np.pi*self.jacobian(a,eccanom)*h)
    else:
      return self.Ma(a)/(2.0*np.pi*self.jacobian(a,eccanom)*h)

  def scale_height(self,X):

    if self.is3D:
      h = self.vertical_structure.get_H(X)
    else:
      h = 1.0

    return h

  def dH(self,X):
 
    if self.is3D:
      dh = self.vertical_structure.get_dH(X)
    else:
      dh = 0.0

    return dh

  def sden(self,X):
    a = X[0]
    eccanom = X[1]

    return self.Ma(a)/(2.0*np.pi*self.jacobian(a,eccanom))

  def pre(self,X):
    a = X[0]
    eccanom = X[1]

    #print a.shape

    rho=self.den(X)
 

    if self.is3D:
      # still need to be carefuly about correct vertical coordinate
      return self.eos.pre(self.entropy(X),rho)*self.vertical_structure.Fp(X)
    else:
      return self.eos.pre(self.entropy(X),rho)


  def temp(self,X):
    a = X[0]
    eccanom = X[1]

    rho=self.den(X)

    if self.is3D:
      return self.eos.temp(self.entropy(X),rho)*self.vertical_structure.Ftemp(X)
    else:
      return self.eos.temp(self.entropy(X),rho)      

   # internal energy or sie?
  def internal_energy(self,X):
    a = X[0]
    eccanom = X[1]

    rho=self.den(X)

    if self.is3D:
      return self.eos.internal_energy(self.entropy(X),rho)*self.vertical_structure.Fsie(X)
    else:
      return self.eos.internal_energy(self.entropy(X),rho)

    # for now include the magnetic field here - later will probably
    # seperate into it's own module.
    # this currently only works for unstratified 3D
  def B(self,X):

    a=X[0]
    eccanom = X[1]

    #print a
    #print eccanom


    # considuring units where mu0=1

    j = self.jacobian(a,eccanom)/a

    # note this doesn't actually work for 3D discs - the vertical structure
    # and mode shapes need to be fixed to do this
    if self.is3D:
      h0=self.vertical_structure.Hcirc(a)
      h = self.vertical_structure.get_H(X)
    else:
      h0=1.0 #this works with the above treatment of the 2D models 
      h = 1.0

      # this is reliant on these being set
      # at the solver level - really they should be
      # set higher up
    if self._solver.Vt==None:
      # beta0 is explicitly for the vertical field
      Vt=0.0

      if self._solver.Vz==None:
        Vz=np.sqrt(2.0/self.beta0)
      else:
        Vz=self._solver.Vz(a)
    else:
      Vt=self._solver.Vt(a)

      if self._solver.Vz==None:
        Vz=0.0
      else:
        Vz=self._solver.Vz(a)

    rho0=self.Ma(a)/(2.0*np.pi*a*h0)
    p0 = self.eos.pre(self.entropy(X),rho0)*h0

    #Bz0 = np.sqrt(p0*one_over_beta0)

    # for now only dealing with isothermal
    # as we haven't set this property anywhere
    gamma=1.0

    Bz0=np.sqrt(gamma*p0)*Vz
    BE0=np.sqrt(gamma*p0)*Vt/a

    eccentricity = self.e(a)

    #Bx = np.zeros_like(a)
    #By = np.zeros_like(a)

    # dimensionless orbital velocity
    vE=1.0/(1.0 - eccentricity*np.cos(eccanom))

    # no orbit crossing field
    Ba=0.0
    BE=BE0*vE*(h0/h)/j
    Bz = Bz0/j

    # transform the horizontal coordinate system


    if self.coord_system == 'CART':

      # need to check these again
 
      Bx=-a*BE*(np.cos(omega)*np.sin(eccanom) + np.sqrt(1.0-eccentricity**2)*np.cos(eccanom)*np.sin(omega))
      By=a*BE*(np.cos(omega)*np.cos(eccanom)*np.sqrt(1.0-eccentricity**2) + np.sin(eccanom)*np.sin(omega))

    elif self.coord_system == 'CYL':
 
      if type(self.omega(a))==type(None):
        omega=0.0
      else:
        omega=self.omega(a)

      #Br
      Bx=a*eccentricity*np.sin(eccanom)*BE

      #Bphi
      #By=((1 - 2.0*eccentricity*np.cos(eccanom)-eccentricity**2)/(1 - eccentricity*np.cos(eccanom)))*BE/np.sqrt(1-eccentricity**2)

      By=(np.sqrt(1.0 - eccentricity*eccentricity)/(1.0 - eccentricity*np.cos(eccanom)))*BE

    elif self.coord_system == 'AE':
      Bx=Ba
      By=BE

    return Bx, By, Bz

    # covarient components of the vector potential
    # mostly for setting up MHD simulations in cylindrical coords
    # as interpolating the B field leads to non-zero div B
  def Avec(self,X):

    a=X[0]
    eccanom = X[1]

    # considuring units where mu0=1

    j = self.jacobian(a,eccanom)/a

    # note this doesn't actually work for 3D discs - the vertical structure
    # and mode shapes need to be fixed to do this
    if self.is3D:
      h0=self.vertical_structure.Hcirc(a)
      h = self.vertical_structure.get_H(X)
    else:
      h0=1.0 #this works with the above treatment of the 2D models 
      h = 1.0

    def _magnetic_integrands(a, y):

      if self._solver.Vt==None:
        # beta0 is explicitly for the vertical field
        Vt=0.0

        if self._solver.Vz==None:
          Vz=np.sqrt(2.0/self.beta0)
        else:
          Vz=self._solver.Vz(a)
      else:
        Vt=self._solver.Vt(a)

        if self._solver.Vz==None:
          Vz=0.0
        else:
          Vz=self._solver.Vz(a)

      X0=[a,[0.0]]

      rho0=self.Ma(a)/(2.0*np.pi*a*h0)
      p0 = self.eos.pre(self.entropy(X0),rho0)*h0

      # for now only dealing with isothermal
      # as we haven't set this property anywhere
      gamma=1.0

      # check the normalisation of this? - assume this is right for now (I think it is)
      Bz0=np.sqrt(gamma*p0)*Vz
      BE0=np.sqrt(gamma*p0)*Vt/a

      # for generality if we decide to generalise to full 3D
      Jcirc=a
      Hcirc=1.0

      return -Hcirc*Jcirc*BE0, Jcirc*Bz0

    # need to check what happens on a grid etc
    
    # check dimensions, should be 2-3d
    y0=[0.0,0.0]
    aspan = [a[0,0],a[-1,0]]

    print(a[0,0].shape)
    print(a[:,0].shape)
    print(aspan)

    print(y0)

    
 
    sol = intg.solve_ivp(_magnetic_integrands,aspan,y0, t_eval=a[:,0])

    #need to outer this?
    eccentricity = np.outer(self.e(a[:,0]),np.ones_like(a[0]))
    aea = np.outer(self.ey(a[:,0]),np.ones_like(a[0]))

    Az=np.outer(sol.y[0],np.ones_like(a[0]))
    Am=np.outer(sol.y[1],np.ones_like(a[0]))


     # transform to a,E coords
    Aa = -aea*np.sin(eccanom)*Am/a
    AE = Am/(1.0 - eccentricity*np.cos(eccanom))
    Az = Az

     # transform to required coordinate system
    Avec = [Aa,AE]

    A1, A2 = self.transform_covariant(Avec,X,coord_system = self.coord_system)
    A3=Az

    #if self.coord_system == 'CART':

    #  raise NotImplementedError("At some point I'll get around to implementing this in Cartesian coords")

      # need to check these again

      #Bx=-a*BE*(np.cos(omega)*np.sin(eccanom) + np.sqrt(1.0-eccentricity**2)*np.cos(eccanom)*np.sin(omega))
      #By=a*BE*(np.cos(omega)*np.cos(eccanom)*np.sqrt(1.0-eccentricity**2) + np.sin(eccanom)*np.sin(omega))

    #elif self.coord_system == 'CYL':

    #  if type(self.omega(a))==type(None):
    #    omega=0.0
  
        #Ar
    #    A1 = -j*((1.0 - eccentricity*np.cos(eccanom))**2)*(aea*np.sin(eccanom)*np.sin(eccanom)/(np.sqrt(1.0-eccentricity*eccentricity)*(1.0 - eccentricity*np.cos(eccanom))))*AE

    #  else:
    #    omega=self.omega(a)

    #    wy = self.e_wy(a)/eccentricity

    #    A1 = -j*((1.0 - eccentricity*np.cos(eccanom))**2)*(wy + aea*np.sin(eccanom)*np.sin(eccanom)/(np.sqrt(1.0-eccentricity*eccentricity)*(1.0 - eccentricity*np.cos(eccanom))))*AE


      #Aphi
    #  A2 = a*j*((1.0 - eccentricity*np.cos(eccanom))**2)*(1.0 - (eccentricity+aea)*np.cos(eccanom))*AE

      #Az
    #  A3 = Az

    #elif self.coord_system == 'AE':
    #  A1=Aa
    #  A2=AE
    #  A3=Az

    return A1, A2, A3    






