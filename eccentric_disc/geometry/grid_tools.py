import numpy as np
import coord_transforms as crdt


####
# Implement functions for generating grids in different coordinate systems
#
####

 ####
 # Takes 2-3 1D arrays of coordinates in a specified coordinate system
 # and yeilds the a E orbital coordinates of a regularly spaced grid in
 # the inputed coordinate system.
 # 
 # e.g. inputing an array of x and y coordinates produces a regularly spaced
 # grid in cartesian coordinates of semimajor axis and eccentric anomaly
 # 
 ####
def aE_meshgrid(coord='CART',orbits,X):

  if coord=='CART':
    Xgrids = np.meshgrid(X,indexing='ij') # I'm not bothering giving this as an option 
                                          # the default meshgrid behaviour is stupid
    return cart2aE(Xgrids,orbits)


  elif coord == 'CYL':
    Xgrids = np.meshgrid(X,indexing='ij')

    return cyl2aE(Xgrids,orbits)

  elif coord == 'AE':
    return np.meshgrid(X,indexing='ij')

  # implemented some not implemented error



 ###
 # Same as aE_meshgrid  but for cartesian coordinates
 # 
 ####
def cart_meshgrid(coord='CART',orbits,X):

  if coord=='CART':
    return np.meshgrid(X,indexing='ij')
  elif coord == 'CYL':
    Xgrids = np.meshgrid(X,indexing='ij')

    return cyl2cart(Xgrids)
  elif coord == 'AE':
    Xgrids = np.meshgrid(X,indexing='ij')

    return aE2cart(Xgrids,orbits)
    
  # implemented some not implemented error



 ###
 # Same as aE_meshgrid  but for cylindrical coordinates
 # 
 ####
def cyl_meshgrid(coord='CART',X):

  if coord=='CART':
    Xgrids = np.meshgrid(X,indexing='ij')
    
    return cart2cyl(Xgrids)
  elif coord == 'CYL':
    return np.meshgrid(X,indexing='ij')
  elif coord == 'AE':
    Xgrids = np.meshgrid(X,indexing='ij')
    
    return aE2cyl(Xgrids,orbits)
    
  # implemented some not implemented error


 ##
 # Add functions to deal with grids which are not orbit regular
 # 
 ##
 
def ordered_flatten(X):

  shape = X.shape

  # there must be a numpy helper function which does this
  origin_indices = np.meshgrid(*(np.arange(N) for N in X.shape),indexing='ij')

  flat_X = X.flatten()
  flat_origin_indices = (origin_index.flatten() for origin_index in origin_indices)  

  sort_indices = np.argsort(flat_X)

  sorted_X = flat_X[sort_indices]

  retar = []
  indexar = []

  count=0
  for i in xrange(len(sorted_X)):
 
    indexar += [tuple([count]+[flat_origin_index[i] for flat_origin_index in flat_origin_indices])]

    if retar[-1] != sorted_X[i]:
      retar+=[sorted_X[i]]
      count+=1
    
  return retar, indexar
    

def flat_data_to_grid(X,index_array):

  ndim = len(index_array[0])-1
 
  shape = [0]*ndim

  # really bad way of doing this
  for i in xrange(ndim):
    for j in xrange(len(index_array)):
      if index_array[j][i+1] > shape[i]:
        shape[i]=index_array[j][i+1]

  retar = np.zeros(shape)

  for i in xrange(len(X)):

    retar[*index_array[i][1:]] = X[index_array[i][0]] 

  return retar






 # should probably put as proper unit tests
 # at some point
if __name__ == '__main__':

  pass

