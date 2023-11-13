import numpy as np
#from . import EccentricDiscGamma2_2D
#from . import  EccentricDiscIso3D 
from . import dgeoFIsodf, dgeoFIsode, ddgeoFIsodf2, ddgeoFIsodfde
from . import dgeoFGamma2_2Ddf, dgeoFGamma2_2Dde, ddgeoFGamma2_2Ddf2, ddgeoFGamma2_2Ddfde
from . import dgeoFGamma2_2Ddf_AntiAligned, dgeoFGamma2_2Dde_AntiAligned, ddgeoFGamma2_2Ddf2_AntiAligned, ddgeoFGamma2_2Ddfde_AntiAligned

# class containing the eccentric mode equation in standard form for a shooting method
# self.forcing_frequency, self.csound2, self.dlnHcirc_dlna should be functions
# which specify the force precession frequency, circular counds speed squared and background
# pressure gradient (dlnHcirc_dlna = d ln (a P)/d ln a)
class EccentricDiscIsoMHD_2D( object ):

  def __init__(self):

    self.forcing_frequency = None # function pointer, function of (a,e)
    self.csound2 = None

    # function descibing the slope
    # of the internal energy in the disc
    # i.e. d ln(H_a^circ)/d ln a
    self.dlnHcirc_dlna = None
 
    self.GM = 1.0

    # referance plasma beta - treat as constant for now
    self.beta0 = 1.0

    # alternatively can set as a ratio of alfven velocity to sound speed
    # Vz = v_a^z/c_s
    # Vt = a v_a^{E}/c_s 
    self.Vz = None 
    self.Vt = None
    
    # derivative of the above
    self.dVz_da = None
    self.dVt_da = None


  # equation of motion, accepts precession frequency as an
  # argument
  def equation_of_motion(self,a,X,omega):

    e = X[0]
    ea = X[1]
    f = e+a*ea

    # note only considuring isothermal so gamma=1

    if self.Vt==None:
      # beta0 is explicitly for the vertical field
      one_over_beta0t=0.0
      
      if self.Vz==None:
        one_over_beta0z=1.0/self.beta0
      else:
        one_over_beta0z=0.5*self.Vz(a)*self.Vz(a)
    else:
      one_over_beta0t=0.5*self.Vt(a)*self.Vt(a)
      
      if self.Vz==None:
        one_over_beta0z=0.0
      else:
        one_over_beta0z=0.5*self.Vz(a)*self.Vz(a)
       
    mag_weightt=one_over_beta0t/(1.0 + one_over_beta0t + one_over_beta0z)
    mag_weightz=one_over_beta0z/(1.0 + one_over_beta0t + one_over_beta0z)

    # don't want 3D as not clear it's consistent
    dFde = (1.0-mag_weightt-mag_weightz)*dgeoFIsode(e,f,is3D=False) + mag_weightt*dgeoFGamma2_2Dde_AntiAligned(e,f) + mag_weightz*dgeoFGamma2_2Dde(e,f)

    # these are the same in 3d and 2d theory
    dFdf = (1.0-mag_weightt-mag_weightz)*dgeoFIsodf(e,f) + mag_weightt*dgeoFGamma2_2Ddf_AntiAligned(e,f) + mag_weightz*dgeoFGamma2_2Ddf(e,f)
    ddFdf2 = (1.0-mag_weightt-mag_weightz)*ddgeoFIsodf2(e,f) + mag_weightt*ddgeoFGamma2_2Ddf2_AntiAligned(e,f) + mag_weightz*ddgeoFGamma2_2Ddf2(e,f)
    ddFdfde = (1.0-mag_weightt-mag_weightz)*ddgeoFIsodfde(e,f) + mag_weightt*ddgeoFGamma2_2Ddfde_AntiAligned(e,f) + mag_weightz*ddgeoFGamma2_2Ddfde(e,f)

    n = np.sqrt(self.GM/(a**3))
    
    omegaf = self.forcing_frequency(a,e)

    fastspeed2 = self.csound2(a)*(1.0 + one_over_beta0t + one_over_beta0z)

    eaa = -(2.0/a)*ea + ((omega - omegaf)*n*a*a*e/(fastspeed2*np.sqrt(1 - e*e)) + dFde - a*ea*ddFdfde  - self.dlnHcirc_dlna(a)*dFdf)/(a*a*ddFdf2)

    if self.Vt!=None:

      dFmagdf_AntiAligned=dgeoFGamma2_2Ddf_AntiAligned(e,f)

      vt_gradient_coeff=a*self.Vt(a)*self.dVt_da(a)/(1.0 + one_over_beta0t + one_over_beta0z)

      eaa+= - vt_gradient_coeff*dFmagdf_AntiAligned/(a*a*ddFdf2)


    if self.Vz!=None:

      dFmagdf=dgeoFGamma2_2Ddf(e,f)

      vz_gradient_coeff=a*self.Vz(a)*self.dVz_da(a)/(1.0 + one_over_beta0t + one_over_beta0z)

      eaa+= - vz_gradient_coeff*dFmagdf/(a*a*ddFdf2)

    return [ea,eaa]




