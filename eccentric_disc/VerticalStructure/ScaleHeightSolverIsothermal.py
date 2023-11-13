import numpy as np
from scipy.integrate import simps, solve_ivp
import scipy.optimize as opt
from .ScaleHeightSolverBase import ScaleHeightSolverBase

class ScaleHeightSolverIsothermal( ScaleHeightSolverBase ):
 
  def __init__( self ):
    #ScaleHeightSolverBase.__init__(self)
    super(ScaleHeightSolverIsothermal,self).__init__()

    self.Irradiated = False

    self.visc=False
    self.alphab=0.0

  def _eom(self,x,y):
    h = y[0]
    dh = y[1]

    temp=1.0
    if self.Irradiated:
      temp=1.0/np.sqrt(1.0 - self.e*np.cos(x))

     #be able to seperately calculate viscous terms
    if self.visc:
      div=self.horizontal_divergence(x)
      temp=temp*(1.0 - self.alphab*(div + dh/h))
    
    # x is eccentric anomaly
    ddh = (temp/h) - ((1.0 - self.e*np.cos(x))**(-3))*h

    return np.array([(1.0 - self.e*np.cos(x))*dh,(1.0 - self.e*np.cos(x))*ddh]) 

# check this 

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    h00=1.0

    solver = ScaleHeightSolverIsothermal()

    solver.q = 0.1

    ecc = np.linspace(0.0,2.0*np.pi,200)

    for e in [0.5,0.6,0.61]:

      try:

        solver.e = e

        solver.h0 = h00

        res = solver.solve()

        soly = solver(ecc)

        plt.plot(ecc,soly[0],'-',linewidth=2,label='e = '+str(e))

        h00 = soly[0,0]
      except:
        continue

    plt.xlabel('E',fontsize=18)
    plt.ylabel('h',fontsize=18)

    plt.legend(loc='best')

    plt.show()




