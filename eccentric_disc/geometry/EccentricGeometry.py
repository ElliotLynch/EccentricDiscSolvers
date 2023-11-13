"""
EccentricGeometry.py

File containing classes that deal with geometrical properties of 
eccentric discs and orbital coordinate systems. Includes (limited 
at present) capability to transform between orbital coordinate 
systems. And various functions for tensor calculus in eccentric 
coordinate systems.

TODO:
 * put common geometric functions in parent class.
 * write an orbital coordinates module, in cython for speed ?
 * local orbit v.s. disc coords?

"""
import numpy as np

 # check this doesn't result in a circular referance?
from . import coord_transforms as crdt


# definining these here so I can use mathematica
# cform inputs for functions
def Power(a,p):
  return a**p

 # class containing (mostly) non-tensoral geometrical
 # properties, for use in ideal fluid solvers which
 # typically do not require full tensor calculus
class EccentricGeometry( object ):
  """
  class containing (mostly) non-tensoral geometrical properties
  
  This class contains (mostly) non-tensoral geometrical
  properties, for use in ideal fluid solvers which
  typically do not require full tensor calculus
  
  Attributes
  ----------
    e : float
        orbit eccentricity
    q : float
        orbit intersection/nonlinearity parameter
    alpha : float
        angle controling relative contribution of eccentricity gradient and twist to q
    ey : float
        eccentricity gradient, a de/da
    e_wy : float
        disc twist, a e d omega/da

  Methods
  -------
  
  radius(Ecc)
      Calculates the cylindrical radius from the coordinate origin to the point on the orbit 
      with eccentric anomaly Ecc
  angular_vel(self,Ecc)
      Calculates the angular velocity at the point on the orbit with eccentric anomaly Ecc
  jacobian(Ecc)
      Calculates the Jacobian determinant of the (Lambda,lambda) canonical coordinate system 
      at the point on the orbit with eccentric anomaly Ecc     
  horizontal_divergence(Ecc)   
     Calculates the horizontal velocity divergence of the orbits at the point on the orbit with 
     eccentric anomaly Ecc
  aE_2_lamphi(self,Ecc)
      coverts from the (a,E) coordinate system of Ogilvie and Lynch (2019) 
      to the (lambda,phi) coordinate system of Ogilvie (2001)
     
  """
  
  def __init__(self):

    #default to circular
    self._e = 0.0
    self._q = 0.0
    self._alpha = 0.0
    self._ey = 0.0
    self._e_wy = 0.0

  # make sure things are recalculated
  # correctly. Using getter/setter methods
  # to ensure all geometrical properties are
  # correctly updated when the user sets one
  @property
  def e(self):
    return self._e

  @property
  def q(self):
    return self._q

  @property
  def alpha(self):
    return self._alpha


  # need to change how _recalculate_derivatives
  # works as this has issues when used with arrays

  @e.setter
  def e(self, e):
    self._e = e
    self._recalculate_derivatives()

  @q.setter
  def q(self, q):
    self._q = q
    self._recalculate_derivatives()

  @alpha.setter
  def alpha(self, alpha):
    self.__alpha = alpha
    self._recalculate_derivatives()


  #I don't think these need to be any different in practice
  @property
  def ey(self):
    return self._ey

  @property
  def e_wy(self):
    return self._e_wy

  @ey.setter
  def ey(self, ey):
    self._ey = ey
    self._recalculate_qalpha()

  @e_wy.setter
  def e_wy(self, e_wy):
    self._e_wy = e_wy
    self._recalculate_qalpha()

  def _recalculate_derivatives(self):

    self._ey = (1 - self.e**2)*self.q*np.cos(self.alpha)/(1 + self.e*self.q*np.cos(self.alpha))
    self._e_wy = self.ey*np.tan(self.alpha)/(np.sqrt(1 - self.e**2))

  # want to include checks for
  # bounded/intersecting orbits?
  def _recalculate_qalpha(self):

    qnum2 = self.ey**2 + (1 - self.e*self.e)*self.e_wy*self.e_wy
    qdenum = 1 - self.e*(self.e + self.ey)

    #should there be some tolerance interval for this?
    
    if np.isscalar(qnum2):
      if qnum2==0.0:
        self._q = 0.0
        self._alpha = 0.0
      else:
        self._q = np.sqrt(qnum2)/qdenum
        self._alpha = np.arccos(self.q*self.ey*qdenum/qnum2)
    else:
      self._q = np.where(qnum2==0,0.0,np.sqrt(qnum2)/qdenum)
      self._alpha = np.where(qnum2==0,0.0,np.arccos(self.q*self.ey*qdenum/qnum2))


  def radius(self,Ecc):
    """Calculates the cylindrical radius of the point with eccentric anomaly Ecc.

    Calculates the cylindrical radius from the coordinate origin to the point on the orbit 
    with eccentric anomaly Ecc. Nondimensionalised such that the semimajor axis a=1

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    # setting a = 1
    return 1 - self.e*np.cos(Ecc)

  def true_anom(self,Ecc):
    
    xrot = np.cos(Ecc) - self.e
 
    yrot = np.sqrt(1.0 - self.e*self.e)*np.sin(Ecc)

    return np.arctan2(yrot,xrot)
 


  def angular_vel(self,Ecc):
     """Calculates the angular velocity of the point with eccentric anomaly Ecc.

    Calculates the angular velocity at the point on the orbit with eccentric anomaly Ecc.
    Nondimensionalised such that the mean motion n=1.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """
     
     return np.sqrt(1 - self.e*self.e)/((1 - self.e*np.cos(Ecc))**2) 
 
  

  def jacobian(self,Ecc):
    """Calculates the (dimensionless) Jacobian determinant of the point with eccentric anomaly Ecc.

    Calculates the (dimensionless) Jacobian determinant of the (Lambda,lambda) canonical 
    coordinate system at the point on the orbit with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    return (1 - self.e*(self.e + self.ey))/np.sqrt(1 - self.e**2) - self.ey*np.cos(Ecc)/np.sqrt(1 - self.e**2) - self.e_wy*np.sin(Ecc)

  def horizontal_divergence(self,Ecc):
    """Calculates the horizontal velocity divergence of the point on the orbit with eccentric anomaly Ecc.

    Calculates the horizontal velocity divergence of the orbits at the point on the orbit with 
    eccentric anomaly Ecc. Nondimensionalised such that the mean motion n=1

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """
    
    j = self.jacobian(Ecc)

    
    return (1.0/(1.0 - self.e*np.cos(Ecc)))*(self.ey*np.sin(Ecc)/np.sqrt(1 - self.e**2) - self.e_wy*np.cos(Ecc))/j

   # coverts from the (a,E) coordinate system of Ogilvie and Lynch (2019) 
   # to the (lambda,phi) coordinate system of Ogilvie (2001)
  def aE_2_lamphi(self,Ecc):
     """Converts from the (a,E) coordinate system to the (lambda,phi) coordinate system.

     coverts from the (a,E) coordinate system of Ogilvie and Lynch (2019) 
     to the (lambda,phi) coordinate system of Ogilvie (2001).

     Parameters
     ----------
     Ecc : float
       Eccentric anomaly of the point
        
     """
      
     ey = self.ey
     e_wy = self.e_wy
     e = self.e

     # check regulerisation

     lal = ey*(1 - e*e)/(1.0 - e*e - 2.0*e*ey)
     e_lol = e_wy*(1 - e*e)/(1.0 - e*e - 2.0*e*ey)

     # above appears to give correct jacobian but with an extra (1 - e**2 - 2 e ey)
     # factor possibly from dlam/da ?

     # returning cos/sin of true anomaly
     # instead of true anomaly or the polar angle phi
     # as is more useful

     cosf = (np.cos(Ecc) - e)/(1.0 - e*np.cos(Ecc))
     sinf = np.sqrt(1 - e*e)*np.sin(Ecc)/(1 - e*np.cos(Ecc))

     return e, lal, e_lol, cosf, sinf




 # list of methods is a bit too long to include?
class EccentricTensorCalc( EccentricGeometry ):
  """
  class containing tensoral geometrical properties for an eccentric orbit.
  
  This class includes methods to calculate various
  quantities from tensor calculus in an eccentric
  orbital coordinate system (e.g. metric tensor components). 
  Inherits geometrical properties from EccentricGeometry.
  
  """
  
  def dradius_dlam(self,Ecc):
    """Calculates dR/d(semilatis rectum) of the point with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    e, lel, e_lol, cosf, sinf = self.aE_2_lamphi(Ecc)

    return (1.0 + (e - lel)*cosf - e_lol*sinf)/((1 + e*cosf)**2)

  # setting a = n = 1 throughout

  def dradius_dphi(self,Ecc):
    """Calculates dR/d(azimuthal angle) of the point with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    e, lel, e_lol, cosf, sinf = self.aE_2_lamphi(Ecc)

    return e*(1- e*e)*sinf/((1 + e*cosf)**2)

  #velocity derivatives
  #Omega=angular velocity
  
  def dradius_dphiphi(self,Ecc):
    """Calculates d^2R/d(azimuthal angle)^2 of the point with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    e, lel, e_lol, cosf, sinf = self.aE_2_lamphi(Ecc)

    return (1 - e*e)*e*(cosf + 2*e*sinf*sinf/(1 + e*cosf))*((1 + e*cosf)**-2)

  def dradius_dphilam(self,Ecc):
    """Calculates d^2R/d(azimuthal angle)d(semilatis rectum) of the point with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    e, lel, e_lol, cosf, sinf = self.aE_2_lamphi(Ecc)

    return ((e+lel)*sinf-e_lol*cosf)*((1 + e*cosf)**-2) - 2*e*sinf*(lel*cosf + e_lol*sinf)*((1 + e*cosf)**-3)


  def dradius_dlamlam(self,Ecc):
    """Calculates d^2R/d(semilatis rectum)^2 of the point with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    e, lel, e_lol, cosf, sinf = self.aE_2_lamphi(Ecc)

    #return ((e+lel)*sinf-e_lol*cosf)*((1 + e*cosf)**-2) - 2*e*sinf*(lel*cosf + e_lol*sinf)*((1 + e*cosf)**-3)
    #possibly change to a warning a calculate with finite difference?
    raise Exeption("requires second derivatives")

  # this includes the eccentricity dependance from the (GM/lam^3)^1/2 term
  def dOmega_dsemailatus(self,Ecc):
    """Calculates d(azimuthal velocity)/d(azimuthal angle) of the point with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    e, lel, e_lol, cosf, sinf = self.aE_2_lamphi(Ecc)

    #return ((Power(-1 + Power(e,2),2)*(3 - 3*Power(e,2) - 2*(e*ey + 2*e_wy) - 4*(ey - e*e_wy)*np.cos(Ecc)))/(2.*(-1 + Power(e,2) + 2*e*ey)*Power(-1 + e*np.cos(Ecc),2)))

    return (-1.5*((1 + e*cosf)**2) + 2.0*(1 + e*cosf)*(lel*cosf + e_lol*sinf))*((1 - e*e)**-2.5)

  def dOmega_dTrueAnom(self,Ecc):
    """Calculates d(azimuthal velocity)/d(true anomaly) of the point with eccentric anomaly Ecc.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    e, lel, e_lol, cosf, sinf = self.aE_2_lamphi(Ecc)

    #return ((-2*e*Power(1 - Power(e,2),1.5)*np.sin(Ecc))/Power(-1 + e*np.cos(Ecc),2))

    #modified
    return -2.0*e*sinf*(1.0 + e*cosf)*((1 - e*e)**-1.5)

  # Metric Tensor (note in semilatus rectum coords so somewhat inconsistent)
  def glamlam(self,Ecc):
    """Calculates lam lam metric tensor of the (lam,phi) coordinate system.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    Rlam = self.dradius_dlam(Ecc)
    #return ((Power(-1 + e*np.cos(Ecc),2)*Power(-1 + Power(e,2) + e*ey + ey*np.cos(Ecc) + np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc),2))/(Power(-1 + Power(e,2),2)*Power(-1 + Power(e,2) + 2*e*ey,2)))
    
    return Rlam*Rlam


  def glamphi(self,Ecc):
    """Calculates lam phi component of the metric tensor of the (lam,phi) coordinate system.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    Rlam = self.dradius_dlam(Ecc)
    Rphi = self.dradius_dphi(Ecc)

    #return (-((e*Power(-1 + e*np.cos(Ecc),2)*np.sin(Ecc)*(-1 + Power(e,2) + e*ey + ey*np.cos(Ecc) + np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc)))/(Power(1 - Power(e,2),1.5)*(-1 + Power(e,2) + 2*e*ey)))*(1-e*e))

    return Rlam*Rphi 

  def gphiphi(self,Ecc):
    """Calculates phi phi component of the metric tensor of the (lam,phi) coordinate system.

    Parameters
    ----------
    Ecc : float
      Eccentric anomaly of the point
        
    """

    R = self.radius(Ecc)
    Rphi = self.dradius_dphi(Ecc)    

    #return ((Power(-1 + e*np.cos(Ecc),3)*(1 + e*np.cos(Ecc)))/Power(-1 + Power(e,2),2))
    return R*R + Rphi*Rphi

    #Inverse Metric Tensor

  def invglamlam(self,Ecc):
   """Calculates lam lam component of the inverse metric tensor of the (lam,phi) coordinate system.

   Parameters
   ----------
   Ecc : float
      Eccentric anomaly of the point
        
   """
   
   R = self.radius(Ecc)
   Rlam = self.dradius_dlam(Ecc)
   Rphi = self.dradius_dphi(Ecc)

    #return (((-1 + Power(e,2))*Power(-1 + Power(e,2) + 2*e*ey,2)*(1 + e*np.cos(Ecc)))/((-1 + e*np.cos(Ecc))*Power(-1 + Power(e,2) + e*ey + ey*np.cos(Ecc) + np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc),2)))

   return (R*R + Rphi*Rphi)/(R*R*Rlam*Rlam)

  def invglamphi(self,Ecc):
    """Calculates lam phi component of the inverse metric tensor of the (lam,phi) coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """
 
    R = self.radius(Ecc)
    Rlam = self.dradius_dlam(Ecc)
    Rphi = self.dradius_dphi(Ecc)

    #return (1-e*e)*(-((e*np.sqrt(1 - Power(e,2))*(-1 + Power(e,2) + 2*e*ey)*np.sin(Ecc))/(Power(-1 + e*np.cos(Ecc),2)*(-1 + Power(e,2) + e*ey + ey*np.cos(Ecc) + np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc)))))

    return -Rphi/(R*R*Rlam)

  def invgphiphi(self,Ecc):
    """Calculates phi phi component of the inverse metric tensor of the (lam,phi) coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """
   
    R = self.radius(Ecc)
 
    #return (Power(1-e*e,2)/(Power(-1 + e*np.cos(Ecc),2)))

    return 1.0/(R*R)

  #Christoffel Symbols

  def lamGamma_lamlam(self,Ecc):
    Rlam = self.dradius_dlam(Ecc)
 
    Rlamlam = self.dradius_dlamlam(Ecc)

    return Rlamlam/Rlam

  def lamGamma_lamphi(self,Ecc):
    R = self.radius(Ecc)
    Rlam = self.dradius_dlam(Ecc)
    Rphi = self.dradius_dphi(Ecc)
 
    Rlamphi = self.dradius_dphilam(Ecc)

    return Rlamphi/Rlam - Rphi/R

  def lamGamma_phiphi(self,Ecc):
    R = self.radius(Ecc)
    Rlam = self.dradius_dlam(Ecc)
    Rphi = self.dradius_dphi(Ecc)
 
    Rphiphi = self.dradius_dphiphi(Ecc)

    return -(R*R + 2.0*Rphi*Rphi - R*Rphiphi)/(R*Rlam)

  def phiGamma_lamlam(self,Ecc):
    return np.zeros_like(Ecc)

  def phiGamma_lamphi(self,Ecc):
    R = self.radius(Ecc)
    Rlam = self.dradius_dlam(Ecc)

    return Rlam/R

  def phiGamma_phiphi(self,Ecc):
    R = self.radius(Ecc)
    Rphi = self.dradius_dphi(Ecc)

    return 2.0*Rphi/R

  # note these are the dimensional forms of the shear tensor with n=1

  #contravarient shear
  def slamlam(self,Ecc):
    """Calculates lam lam component of the rate of strain tensor in the (lam,phi) coordinate system.
   
    Calculates lam lam component of the contravarient rate of strain tensor in the (lam,phi) 
    coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """
 
    R = self.radius(Ecc)
    Rlam = self.dradius_dlam(Ecc)
    Rphi = self.dradius_dphi(Ecc)

    Rlamphi = self.dradius_dphilam(Ecc)
    Rphiphi = self.dradius_dphiphi(Ecc)

    Omega = self.angular_vel(Ecc)
    

    #return (-((Power(-1 + Power(e,2),2)*Power(-1 + Power(e,2) + 2*e*ey,2)*(-(e*(-1 + Power(e,2))*e_wy) + e*(-1 + Power(e,2))*e_wy*Power(np.cos(Ecc),2) + Power(e,2)*(-1 + Power(e,2))*e_wy*Power(np.cos(Ecc),3) + np.sqrt(1 - Power(e,2))*(-e + Power(e,3) - ey + 2*Power(e,2)*ey)*np.sin(Ecc)- Power(e,2)*np.sqrt(1 - Power(e,2))*ey*Power(np.sin(Ecc),3) + np.cos(Ecc)*(e_wy - Power(e,2)*e_wy + np.sqrt(1 - Power(e,2))*ey*np.sin(Ecc))))/(Power(-1 + e*np.cos(Ecc),3)*Power(-1 + Power(e,2) + e*ey + ey*np.cos(Ecc) + np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc),3))))

    return (Rlamphi*(R**3) + R*Rlamphi*Rphi*Rphi + Rlam*(Rphi**3) - R*Rlam*Rphi*Rphiphi)*Omega/((R*Rlam)**3)


  def slamphi(self,Ecc):
    """Calculates lam phi component of the rate of strain tensor in the (lam,phi) coordinate system.
  
    Calculates lam phi component of the contravarient rate of strain tensor in the (lam,phi) 
    coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """

    R = self.radius(Ecc)
    Rlam = self.dradius_dlam(Ecc)
    Rphi = self.dradius_dphi(Ecc)

    Rlamphi = self.dradius_dphilam(Ecc)
    Rphiphi = self.dradius_dphiphi(Ecc)

    Omega = self.angular_vel(Ecc)
    Omegalam = self.dOmega_dsemailatus(Ecc)

    #return (1 - e*e)*(((-1 + e)*Power(1 + e,2)*(-1 + Power(e,2) + 2*e*ey)*(3 - 3*e - 3*Power(e,2) + 3*Power(e,3) + 8*Power(e,2)*ey + 4*e*(1 + e)*ey*Power(np.cos(Ecc),2) + 4*(-1 + e)*np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc) + np.cos(Ecc)*(e - Power(e,3) + Power(e,4) - 4*ey - 4*e*ey - Power(e,2)*(1 + 8*ey) + 4*(-1 + e)*e*np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc))))/(4.*Power(-1 + e*np.cos(Ecc),3)*Power(-1 + Power(e,2) + e*ey + ey*np.cos(Ecc) + np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc),2)))

    return ((R*R + Rphi*Rphi)*Omegalam + (Rlam*Rphiphi - Rphi*Rlamphi)*Omega)/(2.0*R*R*Rlam*Rlam)

  def sphiphi(self,Ecc):
    """Calculates phi phi component of the rate of strain tensor in the (lam,phi) coordinate system.

    Calculates phi phi component of the contravarient rate of strain tensor in the (lam,phi) 
    coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """
 
    R = self.radius(Ecc)
    Rlam = self.dradius_dlam(Ecc)
    Rphi = self.dradius_dphi(Ecc)
 
    Omega = self.angular_vel(Ecc)
    Omegalam = self.dOmega_dsemailatus(Ecc)

    #return Power(1 - e*e,2)*((e*Power(1 - Power(e,2),1.5)*np.sin(Ecc)*(-1 + Power(e,2) + 2*ey*np.cos(Ecc) + 2*np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc)))/(2.*Power(-1 + e*np.cos(Ecc),4)*(-1 + Power(e,2) + e*ey + ey*np.cos(Ecc) + np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc))))
  
    return -Rphi*(R*Omegalam + Rlam*Omega)/(R*R*R*Rlam)

  #covarient shear
  def invslamlam(self,Ecc):
    """Calculates lam lam component of the inverse rate of strain tensor in the (lam,phi) coordinate system.

    Calculates lam lam component of the covarient (inverse) rate of strain tensor in the (lam,phi) 
    coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """

    glamlam = self.glamlam(Ecc)
    glamphi = self.glamphi(Ecc)
 
    slamlam = self.slamlam(Ecc)
    slamphi = self.slamphi(Ecc)
    sphiphi = self.sphiphi(Ecc)

    #return ((-3*e*e_wy + 6*Power(e,3)*e_wy - 3*Power(e,5)*e_wy + 2*ey*e_wy - 2*Power(e,4)*ey*e_wy + (-1 + Power(e,2))*(-e + Power(e,3) + 4*ey)*e_wy*Power(np.cos(Ecc),2) - np.sqrt(1 - Power(e,2))*(2*ey - Power(e,2)*ey - Power(e,4)*ey - e*(1 + 2*Power(ey,2)) - 2*e*(-e*e + Power(e_wy,2)) + Power(e,3)*(-e*e + 2*Power(e_wy,2)))*np.sin(Ecc) + np.cos(Ecc)*(2*Power(-1 + Power(e,2),2)*e_wy + np.sqrt(1 - Power(e,2))*(e*(-1 + Power(e,2))*ey - 2*Power(e_wy,2))*np.sin(Ecc)) + np.sqrt(1 - Power(e,2))*Power(ey,2)*np.sin(2*Ecc) + Power(e,2)*np.sqrt(1 - Power(e,2))*Power(e_wy,2)*np.sin(2*Ecc))/(2.*(-1 + Power(e,2))*Power(-1 + Power(e,2) + 2*e*ey,2)))

    return glamlam*glamlam*slamlam + 2.*glamlam*glamphi*slamphi + glamphi*glamphi*sphiphi    



  def invslamphi(self,Ecc):
    """Calculates lam phi component of the inverse rate of strain tensor in the (lam,phi) coordinate system.
 
    Calculates lam phi component of the covarient (inverse) rate of strain tensor in the (lam,phi) 
    coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """

    glamlam = self.glamlam(Ecc)
    glamphi = self.glamphi(Ecc)
    gphiphi = self.gphiphi(Ecc) 

    slamlam = self.slamlam(Ecc)
    slamphi = self.slamphi(Ecc)
    sphiphi = self.sphiphi(Ecc)

   # return (((-1 + e*np.cos(Ecc))*(-3 + 3*Power(e,2) + 4*e*ey + (-e + Power(e,3) + 4*ey)*np.cos(Ecc) + 4*np.sqrt(1 - Power(e,2))*e_wy*np.sin(Ecc)))/(4.*(1 - e*e)*(-1 + Power(e,2) + 2*e*ey)))

    return glamlam*glamphi*slamlam + (glamphi*glamphi + glamlam*gphiphi)*slamphi + glamphi*gphiphi*sphiphi 

  def invsphiphi(self,Ecc):
    """Calculates phi phi component of the rate of strain tensor in the (lam,phi) coordinate system.

    Calculates phi phi component of the covarient (inverse) rate of strain tensor in the (lam,phi) 
    coordinate system.

    Parameters
    ----------
    Ecc : float
       Eccentric anomaly of the point
        
    """
 
    glamphi = self.glamphi(Ecc)
    gphiphi = self.gphiphi(Ecc)

    slamlam = self.slamlam(Ecc)
    slamphi = self.slamphi(Ecc)
    sphiphi = self.sphiphi(Ecc)
    #return (e*np.sqrt(1 - Power(e,2))*(-1 + e*np.cos(Ecc))*np.sin(Ecc)/Power(1 - e*e,2))

    return glamphi*glamphi*slamlam + 2.0*gphiphi*glamphi*slamphi + gphiphi*gphiphi*sphiphi 


    # transforms covarient components of a vector from a,E to desired coordinate system
  def transform_covariant(self,V,X,coord_system = 'AE',omega=None):

    Va = V[0]
    VE = V[1] 

    a = X[0]
    Ecc = X[1]

    ey = self.ey
    #e_wy = self.e_wy
    e = self.e

    J = a*self.jacobian(Ecc)

      # strictly this doesn't work for constant omege ne 0
    if omega==None:

      dxda = (np.cos(Ecc) - e) - ey
      dxdE = -a*np.sin(Ecc)

      dyda = np.sqrt(1 - e*e)*np.sin(Ecc) - e*ey*np.sin(Ecc)/np.sqrt(1 - e*e)
      dydE = a*np.sqrt(1 - e*e)*np.cos(Ecc)

    else:
      wy = self.e_wy/e

      dxda = (np.cos(Ecc) - e)*np.cos(omega) -np.sqrt(1 - e*e)*np.sin(Ecc)*np.sin(omega) - ey*np.cos(omega) + e*ey*np.sin(Ecc)*np.sin(omega)/np.sqrt(1 - e*e) -wy*(np.cos(Ecc) - e)*np.sin(omega) - wy*np.sqrt(1 - e*e)*np.sin(Ecc)*np.cos(omega)
      dxdE = -a*np.sin(Ecc)*np.cos(omega) - a*np.sqrt(1 - e*e)*np.cos(Ecc)*np.sin(omega)

      dyda = np.sqrt(1 - e*e)*np.sin(Ecc)*np.cos(omega) + (np.cos(Ecc) - e)*np.sin(omega) - e*ey*np.sin(Ecc)*np.cos(omega)/np.sqrt(1 - e*e) - ey*np.sin(omega) - wy*np.sqrt(1 - e*e)*np.sin(Ecc)*np.sin(omega) + wy*(np.cos(Ecc) - e)*np.cos(omega)
      dydE = a*np.sqrt(1 - e*e)*np.cos(Ecc)*np.cos(omega) - a*np.sin(Ecc)*np.sin(omega)

    if coord_system=='AE':
      V1=Va
      V2=VE
    elif coord_system=='CART':
      V1 = (dydE*Va-dyda*VE)/J
      V2 = (-dxdE*Va+dxda*VE)/J
    elif coord_system=='CYL':
      #r, phi = crdt.aE2cyl(self,X)

      r=a*self.radius(Ecc)
      f=self.true_anom(Ecc)

      # probably not the sensible way of doing this
      if type(omega)==type(None):
        phi = f
      else:
        phi = f + omega

      V1 = (np.cos(phi)*dydE-np.sin(phi)*dxdE)*Va/J + (-np.cos(phi)*dyda+np.sin(phi)*dxda)*VE/J
      V2 = -r*(np.sin(phi)*dydE+np.cos(phi)*dxdE)*Va/J + r*(np.sin(phi)*dyda+np.cos(phi)*dxda)*VE/J

    return V1, V2

    # transforms contravarient components of a vector from a,E to desired coordinate system
  def transform_contravarient(self,V,X,coord_system = 'AE',omega=None):

    Va = V[0]
    VE = V[1] 

    a = X[0]
    Ecc = X[1]

    ey = self.ey
    #e_wy = self.e_wy
    e = self.e

      # strictly this doesn't work for constant omege ne 0
    if omega==None:
    
      dxda = (np.cos(Ecc) - e) - ey
      dxdE = -a*np.sin(Ecc)

      dyda = np.sqrt(1 - e*e)*np.sin(Ecc) - e*ey*np.sin(Ecc)/np.sqrt(1 - e*e)
      dydE = a*np.sqrt(1 - e*e)*np.cos(Ecc)

    else:
      wy = self.e_wy/e

      dxda = (np.cos(Ecc) - e)*np.cos(omega) -np.sqrt(1 - e*e)*np.sin(Ecc)*np.sin(omega) - ey*np.cos(omega) + e*ey*np.sin(Ecc)*np.sin(omega)/np.sqrt(1 - e*e) -wy*(np.cos(Ecc) - e)*np.sin(omega) - wy*np.sqrt(1 - e*e)*np.sin(Ecc)*np.cos(omega)
      dxdE = -a*np.sin(Ecc)*np.cos(omega) - a*np.sqrt(1 - e*e)*np.cos(Ecc)*np.sin(omega)

      dyda = np.sqrt(1 - e*e)*np.sin(Ecc)*np.cos(omega) + (np.cos(Ecc) - e)*np.sin(omega) - e*ey*np.sin(Ecc)*np.cos(omega)/np.sqrt(1 - e*e) - ey*np.sin(omega) - wy*np.sqrt(1 - e*e)*np.sin(Ecc)*np.sin(omega) + wy*(np.cos(Ecc) - e)*np.cos(omega)
      dydE = a*np.sqrt(1 - e*e)*np.cos(Ecc)*np.cos(omega) - a*np.sin(Ecc)*np.sin(omega)


    if coord_system=='AE':
      V1=Va
      V2=VE
    elif coord_system=='CART':
      V1 = dxda*Va+dxdE*VE
      V2 = dyda*Va+dydE*VE
    elif coord_system=='CYL':
      #r, phi = crdt.aE2cyl(self,X)
       
      r=a*self.radius(Ecc)
      f=self.true_anom(Ecc)

      # probably not the sensible way of doing this
      if type(omega)==type(None):
        phi = f
      else:
        phi = f + omega

      V1 = (np.cos(phi)*dxda+np.sin(phi)*dyda)*Va + (np.cos(phi)*dxdE+np.sin(phi)*dydE)*VE
      V2 = (-np.sin(phi)*dxda+np.cos(phi)*dyda)*Va/r + (-np.sin(phi)*dxdE+np.cos(phi)*dydE)*VE/r

    return V1, V2



if __name__ == '__main__':

  import matplotlib.pyplot as plt

  # best way of turning this into proper tests?

  obj = EccentricTensorCalc()
  
  obj.e = 0.9
  obj.q = 0.9

  eccanom = np.linspace(0.0,2.0*np.pi,200) 

  # divergence scaling test

  div = obj.horizontal_divergence(eccanom)
  div2 = obj.glamlam(eccanom)*obj.slamlam(eccanom) + 2.0*obj.glamphi(eccanom)*obj.slamphi(eccanom) + obj.gphiphi(eccanom)*obj.sphiphi(eccanom)
  div3 = obj.invglamlam(eccanom)*obj.invslamlam(eccanom) + 2.0*obj.invglamphi(eccanom)*obj.invslamphi(eccanom) + obj.invgphiphi(eccanom)*obj.invsphiphi(eccanom)
  

  plt.plot(eccanom,div,'k-',linewidth=2,label=r'$\Delta$')
  plt.plot(eccanom,div2,'b--',linewidth=2,label=r'$g_{i j} S^{i j}$')
  plt.plot(eccanom,div3,'g:',linewidth=2,label=r'$g^{i j} S_{i j}$')

  plt.legend(loc='best',fontsize=18)
 
  plt.show()


  # should be 2
  traceg = obj.invglamlam(eccanom)*obj.glamlam(eccanom) + 2.0*obj.invglamphi(eccanom)*obj.glamphi(eccanom) + obj.invgphiphi(eccanom)*obj.gphiphi(eccanom)

  plt.plot(eccanom,traceg,'k-',linewidth=2)

  plt.show()





