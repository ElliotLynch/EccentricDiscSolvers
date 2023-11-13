import numpy as np
from scipy.integrate import simps, solve_ivp
import scipy.optimize as opt
from .ScaleHeightSolverBase import ScaleHeightSolverBase


N_INTEGRAL_POINTS = 100


class ScaleHeightSolverIdeal( ScaleHeightSolverBase ):
 
  def __init__( self ):
    super(ScaleHeightSolverIdeal,self).__init__()
    #ScaleHeightSolverBase.__init__(self)    

    self.gamma = 2.0

  def _eom(self,x,y):
    h = y[0]
    dh = y[1]

    j = self.jacobian(x)

    # x is eccentric anomaly
    ddh = (h**(-self.gamma))*(j**(1.0-self.gamma)) - ((1.0 - self.e*np.cos(x))**(-3))*h

    return np.array([(1.0 - self.e*np.cos(x))*dh,(1.0 - self.e*np.cos(x))*ddh]) 



# check this 

 # has to be run as "python -m eccentric_disc.VerticalStructure.ScaleHeightSolverIdeal"

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    h00=1.0

    solver = ScaleHeightSolverIdeal()

    solver.q = 0.1

    ecc = np.linspace(0.0,2.0*np.pi,200)

    for e in [0.1,0.5,0.7,0.8,0.9]:

      try:

        solver.e = e

        solver.h0 = h00

        res = solver.solve()

        soly = solver(ecc)

        plt.plot(ecc,soly[0],'-',linewidth=2,label='e = '+str(e)+', q = '+str(solver.q))

        h00 = soly[0,0]
      except:
        continue

    solver.q = 0.5

    for e in [0.1,0.5,0.7,0.8,0.9]:

      try:

        solver.e = e

        solver.h0 = h00

        res = solver.solve()

        soly = solver(ecc)

        plt.plot(ecc,soly[0],'--',linewidth=2,label='e = '+str(e)+', q = '+str(solver.q))

        h00 = soly[0,0]
      except:
        continue

    plt.xlabel('E',fontsize=18)
    plt.ylabel('h',fontsize=18)

    plt.legend(loc='best')

    plt.show()




