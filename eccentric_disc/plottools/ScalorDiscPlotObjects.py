import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.projections import PolarAxes
from matplotlib.figure import Figure

class ScalorDiscAxes( Axes ):

  draw_patch = plt.Circle((0.5, 0.5), 0.5)

  def __init__(self,*args,**kwargs):
    Axes.__init__(self,*args,**kwargs)

    self.crdgrid = None
    self.scalor = None
    self.to_cart_grid = self._to_cart_grid_default

    #self.ax = None
    #self.fig = None

    self.subsample = 1

    self.orbits = None
    self.scalorcontour = None
    self.scalorcontourf = None
    self.scalorpcolormesh = None
 
    #self..set_visible(False)

    args[0].add_axes(self)

  def _gen_axes_patch(self):
    return plt.Circle((0.5, 0.5), 0.5)

  def _to_cart_grid_default(self,grid):
    return grid

  def plot_orbits( self, vmin = None, vmax = None ):

    x, y = self.to_cart_grid(self.crdgrid)

    self.orbits = self.plot(x[::self.subsample].T, y[::self.subsample].T,'k-',linewidth=2)

  # give option to specify scalor here

  def contour( self, Z=None, vmin = None, vmax = None ):

    if type(Z) != type(None):
      self.scalor=Z

    x, y = self.to_cart_grid(self.crdgrid)

    if vmin == None:
      vmin = np.min(self.scalor)
    if vmax == None:
      vmax = np.max(self.scalor)

    self.scalorcontour = Axes.contour(self,x, y, self.scalor) #, vmin=vmin, vmax=vmax)

    self.active_scalor = self.scalorcontour

  def contourf( self, Z=None, vmin = None, vmax = None ):

    if type(Z) != type(None):
      self.scalor=Z

    x, y = self.to_cart_grid(self.crdgrid)

    if vmin == None:
      vmin = np.min(self.scalor)
    if vmax == None:
      vmax = np.max(self.scalor)

    self.scalorcontourf = Axes.contourf(self,x, y, self.scalor) #, vmin=vmin, vmax=vmax)
 
    self.active_scalor = self.scalorcontourf

  def pcolormesh( self, Z=None, vmin = None, vmax = None ):

    if type(Z) != type(None):
      self.scalor=Z

    x, y = self.to_cart_grid(self.crdgrid)

    if vmin == None:
      vmin = np.min(self.scalor)
    if vmax == None:
      vmax = np.max(self.scalor)

    self.scalorpcolormesh = Axes.pcolormesh(self,x, y, self.scalor) #, vmin=vmin, vmax=vmax) 
    self.active_scalor = self.scalorpcolormesh

  def set_colorbar(self,**kwargs):

    # not clear if I'm just reimplimenting an existing functionality of matplotlib

    return self.get_figure().colorbar(self.active_scalor, ax = self, **kwargs)


class ScalorDiscFigure(Figure):

  def __init__(self, *args, **kwargs):
    Figure.__init__(self, *args, **kwargs)

    self.crdgrid = None
    self.to_cart_grid = self._to_cart_grid_default
    self.subsample = 1

  def _to_cart_grid_default(self,grid):
    return grid

  def add_subplot(self,*args,**kwargs):
    ax = Figure.add_subplot(self,*args,frame_on=False,visible=False,**kwargs)
    ax = ScalorDiscAxes(self,ax.get_position())

    ax.crdgrid = self.crdgrid
    ax.to_cart_grid = self.to_cart_grid
    ax.subsample = self.subsample

    return ax

def scalordisc_subplot(nrows=1, ncols=1,sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None,grid=None,to_cart_grid=None,subsample=1,**fg_kwargs):

  fig = plt.figure(FigureClass=ScalorDiscFigure,**fg_kwargs)
  fig.crdgrid = grid

  if type(to_cart_grid)!=type(None):
    fig.to_cart_grid = to_cart_grid
  
  fig.subsample = subsample


  #print to_polar_grid
  #raise

  ax = []

  for plot_number in range(nrows*ncols):
    if subplot_kw != None:
      ax += [fig.add_subplot(nrows,ncols,1 + plot_number,**subplot_kw)]
    else:
      ax += [fig.add_subplot(nrows,ncols,1 + plot_number)]

  return fig, ax


 # has to be run as "python -m eccentric_disc.plottools.ScalorDiscPlotObjects"
 # as __main__ doesn't act within the package otherwise so relative imports
 # do not work.

if __name__ == '__main__':

  from ..geometry import EccentricOrbits  #.aE2cyl(orbits,X)
  from ..geometry.coord_transforms import aE2cart

  amin = 1.0
  amax = 2.0
  e0 = 0.1

  def eorbital(a): 
    return e0*(amax/a) 

  def oorbital(a): 
    return np.zeros_like(a)

  orbits = EccentricOrbits()
  orbits.e = eorbital
  orbits.omega = oorbital

  #return polar coords for plotting

  # ideally write EccentricOrbits in such 
  # a way to interface directly
  def to_cart_grid(grid):
    semimajor = grid[0]
    eccanom = grid[1]
   
    # use this as an oportunity to write the
    # 2d version. also nonasci bs going on for some reason
    return aE2cart(orbits,(semimajor,eccanom))


  grid = np.meshgrid(np.linspace(amin,amax,200),np.linspace(0,2.0*np.pi,200), indexing = 'ij')  

  fig, ax = scalordisc_subplot(2, 3, grid=grid,to_cart_grid=to_cart_grid,subsample=20)

  sigma = 1.0/(1.0 - e0*np.cos(grid[1]))  #disc.surface_density(grid)
  #sie = np.ones_like(grid[0]) #disc.specific_internalenergy(grid)

  #orbits only
  ax[0].scalor = np.log10(sigma)

  ax[0].plot_orbits()
  #ax[0].plot_scalor()
  #ax[0].set_colorbar(orientation='horizontal')

  ax[0].set_title('orbits only')

  #contour
  ax[1].scalor = np.log10(sigma)

  ax[1].plot_orbits()
  ax[1].contour()
  #ax[1].set_colorbar(orientation='horizontal')

  ax[1].set_title('contour')

  #contourf
  ax[2].scalor = np.log10(sigma)

  ax[2].plot_orbits()
  ax[2].contourf()
  ax[2].set_colorbar(orientation='horizontal')

  ax[2].set_title('contourf')

  #pcolormesh
  ax[3].scalor = np.log10(sigma)

  ax[3].plot_orbits()
  ax[3].pcolormesh()
  ax[3].set_colorbar(orientation='horizontal')

  #ax[3].get_xaxis().set_visible(False)

  ax[3].set_title('pcolormesh')

  #orbits only no border
  ax[4].scalor = np.log10(sigma)

  ax[4].plot_orbits()

  ax[4].get_xaxis().set_visible(False)
  ax[4].get_yaxis().set_visible(False)
  ax[4].set_frame_on(False)

  ax[4].set_title('orbits only, no border')

  #pcolormesh no border
  ax[5].scalor = np.log10(sigma)

  ax[5].plot_orbits()
  ax[5].pcolormesh()
  ax[5].set_colorbar(orientation='horizontal')

  ax[5].get_xaxis().set_visible(False)
  ax[5].get_yaxis().set_visible(False)
  ax[5].set_frame_on(False)

  ax[5].set_title('pcolormesh, no border')

  plt.show()
  







