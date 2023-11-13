import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from .. import geometry as egeo 


class BinaryOrbitPlot( object ):

  # want this to be common
  NAZIMUTH_PLOTPOINTS = 500

  def __init__(self,ax=None):
  
    self.e=0.0
    self.a=1.0
    self.omega=0.0 

    #self.linestyle =

    self.qbin = 1.0

    self.ax = None

    self.binary_scale=0.25

  def plot(self,*args,**kwargs):

    eccanom = np.linspace(0.0,2.0*np.pi,self.NAZIMUTH_PLOTPOINTS)

    r = self.a*(1.0 - self.e*np.cos(eccanom))

    xrot = self.a*(np.cos(eccanom) - self.e)

    yrot = self.a*np.sqrt(1.0 - self.e*self.e)*np.sin(eccanom)

    phi = np.arctan2(yrot,xrot) + self.omega

    x = r*np.cos(phi)
    y = r*np.sin(phi)

    r_peri = self.a*(1.0 - self.e*np.cos(0.0))

    xrot_peri = self.a*(np.cos(0.0) - self.e)

    yrot_peri = self.a*np.sqrt(1.0 - self.e*self.e)*np.sin(0.0)

    phi_peri = np.arctan2(yrot_peri,xrot_peri) + self.omega

    x_peri = r_peri*np.cos(phi_peri)
    y_peri = r_peri*np.sin(phi_peri)

    r_apo = self.a*(1.0 - self.e*np.cos(np.pi))

    xrot_apo = self.a*(np.cos(np.pi) - self.e)

    yrot_apo = self.a*np.sqrt(1.0 - self.e*self.e)*np.sin(np.pi)

    phi_apo = np.arctan2(yrot_apo,xrot_apo) + self.omega
    
    x_apo = r_apo*np.cos(phi_apo)
    y_apo = r_apo*np.sin(phi_apo)

    print (x_peri, y_peri)
    print (x_apo, y_apo)

    # not sure what size to do

    mprimary = 1.0/(1.0+self.qbin)
    msecondary = self.qbin/(1.0+self.qbin) 

    # primary
    primary = plt.Circle((x_peri, y_peri), mprimary*self.binary_scale,color='k')

    # secondary
    secondary = plt.Circle((x_apo, y_apo), msecondary*self.binary_scale,color='k')

    print(primary, secondary)

    print(self.ax)

    self.ax.add_patch(primary)
    self.ax.add_patch(secondary)

    # orbit
    return primary, secondary, self.ax.plot(x,y,*args,**kwargs)

    






