import numpy as np

 # write minimum working 
 # example for this right now

class EOSIsothermal( object ):

  # really we don't want to set this up this way

  def pre(self,s,rho):
    return np.exp(s)*rho

  def internal_energy(self,s,rho):
    return np.exp(s)*np.log(rho)

  def temp(self,s,rho):
    return np.exp(s)*np.ones_like(rho)

