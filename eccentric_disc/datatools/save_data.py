import numpy as np
import scipy.interpolate as interp

#def save(s,s_grid,save_grid,save_format='PLUTO',*args):

def save(s,s_grid,save_grid,**kwargs):

  if kwargs['save_format']=='PLUTO':
    # idealy want to put parse args to put in the appropriate PLUTO defaults
    return savePLUTO(s,s_grid,save_grid,**kwargs)
  if kwargs['save_format']=='RADMC':
    return saveRADMC(s,s_grid,save_grid,**kwargs)
  else:  
    raise NotImplementedError("This save format hasn't been implimented yet.")

  # save eccentric mode to a dust density input file
  # to be read into radmc
  # based off script provided by Josh Lovell, originally writen by Seba Marino
def saveRADMC(s,s_grid,save_grid,**kwargs):

  if 'Nspec' in kwargs:
    Nspec=kwargs['Nspec']
  else:
    Nspec=1

  if 'Nth' in kwargs:
    Nth=kwargs['Nth']
  else:
    Nth=2

  if 'Ms' in kwargs:
    Ms=kwargs['Ms']
  else:
    Ms=[1]

  if 'filename' in kwargs:
    filename=kwargs['filename']
  else:
    filename="dust_density.inp"

  #save_dens_nonaxisym(Nspec, Redge, R, Thedge, Th, Phiedge, Phi, Ms, h, sigmaf, *args):

  # args has the arguments that sigmaf needs in the right order

  agrid = s_grid[0]
  rgrid = save_grid[0]
  phigrid = save_grid[1]

  Nr=rgrid.shape[0]+1
  Nphi=phigrid.shape[1]

  
  # how do we setup phiedge?

  dphi=Phiedge[1]-Phiedge[0]
  rho_d=np.zeros((Nspec,(Nth-1)*2,Nphi,Nr-1)) # density field


  #fieldar = np.zeros_like(agrid)
  den = np.zeros((Nphi,Nr))

  es = s.e(agrid)

  eccgrid = np.zeros_like(agrid)

   # need to do it this way to deal with the branch cuts
  for i in range(agrid.shape[0]):
    
    if type(s.omega(agrid[i]))==type(None):
      trueanom=phigrid[i]
    else:
      trueanom=phigrid[i]+s.omega(agrid[i])

    xrot=(np.cos(trueanom)+es[i])/(1.0 + es[i]*np.cos(trueanom))
    yrot=np.sqrt(1.0 - es[i]*es[i])*np.sin(trueanom)/(1.0 + es[i]*np.cos(trueanom))

    eccgrid[i] = np.arctan2(yrot,xrot)

    #eccgrid[i] = np.arccos((np.cos(phigrid[i]) + es[i])/(1.0 + es[i]*np.cos(phigrid[i])))

  grid = [agrid,eccgrid]

  denfield = s["DEN"](grid)

  phis=phigrid[0]
  rs=rgrid[:,0]

  for i in range(phigrid.shape[1]):

    den[i,:] = interp.interp1d(agrid[:,i],denfield[:,i])(rs)

   # need to think about how to incoperate the vertical structure here
  for ia in range(Nspec):
    M_dust_temp= 0.0
    for k in range(Nz):

      rho_d[ia,k,:,:]=den[:,:]
      rho_d[ia,2*(Nth-1)-1-k,:,:]=den[:,:]

      # ok need to think how to preform this integration
      M_dust_temp+=2.0*rho_d[ia,k,j,i]*(dphi)*rho*(Redge[i+1]-Redge[i])*(Thedge[Nth-2-k+1]-Thedge[Nth-2-k])*R[i]*au**3.0

    rho_d[ia,:,:,:]=rho_d[ia,:,:,:]*Ms[ia]/M_dust_temp

  #for ia in xrange(Nspec):
  #  M_dust_temp= 0.0 #np.zeros(Nspec) 

  #  for k in xrange(Nth-1):
  #    for i in xrange(Nr-1):  
  #      for j in xrange(Nphi):
  #        phi=Phi[j]
  #        rho_d[ia,k,j,i]=rho_3d_dens(rho,phi, z,h, sigmaf, *args )
  #        rho_d[ia,2*(Nth-1)-1-k,j,i]=rho_d[ia,k,j,i]

          #M_dust_temp+=2.0*rho_d[ia,k,j,i]*(dphi)*rho*(Redge[i+1]-Redge[i])*(Thedge[Nth-2-k+1]-Thedge[Nth-2-k])*R[i]*au**3.0

    #rho_d[ia,:,:,:]=rho_d[ia,:,:,:]*Ms[ia]/M_dust_temp

  # Save
  dust_d=open(filename,'w')

  dust_d.write('1 \n') # iformat  
  dust_d.write(str((Nr-1)*2*(Nth-1)*(Nphi))+' \n') # iformat n cells 
  dust_d.write(str(Nspec)+' \n') # n species

  for ai in xrange(Nspec):
    for j in range(Nphi):
      for k in range(2*(Nth-1)):
        for i in range(Nr-1):
          dust_d.write(str(rho_d[ai,k,j,i])+' \n')

  dust_d.close()





  # save eccentric mode to pluto file
  # for now just saving 2D, isothermal data
  # later want to add full 3D as an option
def savePLUTO(s,s_grid,save_grid,**kwargs):

  """
   Author:        Loren E. Held, Elliot Lynch
   Date:          WedJun092021

   Description:   This function allows the user to save and eccentric mode as 
                  PLUTO initial conditions (3D arrays of hydro fluid variables 
                  in PLUTO grid coordinates) in a binary file format (dbl) that 
                  can be imported into PLUTO v4.2. The data must be in code units.
  
  """

  # not really the best way of doing this but o well

  if 'Nz' in kwargs:
    Nz=kwargs['Nz']
  else:
    Nz=16

  if 'filename' in kwargs:
    filename=kwargs['filename']
  else:
    filename="data.dbl"

   # may not be the best way of doing this
  gridded_e=False
  if 'gridded_e' in kwargs:
    gridded_e=kwargs['gridded_e']

  fill_copy=False

  # dealing with the edge of the disc, for free bc's (if they work)
  # setting to zero is the best option, but that's not good for
  # circular boundaries.
  if 'fill_value' in kwargs:
    if kwargs['fill_value']=='copy':
      fill_copy=True    
      fill_value=-1.0
    else:
      fill_value=kwargs['fill_value']
  else:
    fill_value='extrapolate'

  #print kwargs['isMHD']

  isMHD=False
  if 'isMHD' in kwargs:
    isMHD=kwargs['isMHD'] 

  vectorPotential=True
  if 'vectorPotential' in kwargs:
    vectorPotential=kwargs['vectorPotential']

  agrid = s_grid[0]
  rgrid = save_grid[0]
  phigrid = save_grid[1]

  Nx=rgrid.shape[0]
  Ny=phigrid.shape[1]

  #fieldar = np.zeros_like(agrid)
  den = np.zeros((Nx,Ny))
  vx1 = np.zeros((Nx,Ny))
  vx2 = np.zeros((Nx,Ny))
  vx3 = np.zeros((Nx,Ny))

  if isMHD:
    bx1 = np.zeros((Nx,Ny))
    bx2 = np.zeros((Nx,Ny))
    bx3 = np.zeros((Nx,Ny))

    Ax1 = np.zeros((Nx,Ny))
    Ax2 = np.zeros((Nx,Ny))
    Ax3 = np.zeros((Nx,Ny))

  if gridded_e:
    es = s.e
  else:
    es = s.e(agrid)

  eccgrid = np.zeros_like(agrid)

   # need to do it this way to deal with the branch cuts
  for i in range(agrid.shape[0]):
    
    if gridded_e:
     if type(s.omega)==type(None):
       trueanom=phigrid[i]
     else:
       trueanom=phigrid[i]+s.omega
    else:
      if type(s.omega(agrid[i]))==type(None):
        trueanom=phigrid[i]
      else:
        trueanom=phigrid[i]+s.omega(agrid[i])

    xrot=(np.cos(trueanom)+es[i])/(1.0 + es[i]*np.cos(trueanom))
    yrot=np.sqrt(1.0 - es[i]*es[i])*np.sin(trueanom)/(1.0 + es[i]*np.cos(trueanom))

    eccgrid[i] = np.arctan2(yrot,xrot)

  grid = [agrid,eccgrid]

  # radius in the orbit regular grid
  radius_grid=agrid*(1.0 - es*np.cos(eccgrid))

  denfield = s["DEN"](grid)
  
  # cylndrical coords
  vr, vphi = s.V(grid)
  
  # changing to rectangular velocities
  vx1field = vr
  vx2field = radius_grid*vphi


  vx3field = np.zeros_like(denfield)

  # if we introduce horizontal field we need to
  # check the components
  if isMHD:
    #bxfields = s.B(grid)
    br, bphi, bz = s.B(grid)

    Ar, Aphi, Az  = s.Avec(grid)

    # presumably rectangular variables
    Afields = [Ar, Aphi/radius_grid, Az]

    #Afields = [Ar, Aphi, Az]

    bxfields=[br, radius_grid*bphi, bz]
 
  #raise

  phis=phigrid[0]
  rs=rgrid[:,0]

  # try this to see how much it effects
  #ey=s.ey(agrid)
  #Le = ey*(1 - es*es)/(1.0 - es*es - 2.0*es*ey)
  #denfield = (1. - es**2)**(3./2.)*(1. + es*np.cos(phigrid))/(1. + (es - Le)*np.cos(phigrid))

  for i in range(phigrid.shape[1]):

    # extrapolate out of bounds?
    den[:,i] = interp.interp1d(radius_grid[:,i],denfield[:,i], kind='linear',bounds_error=False, fill_value=fill_value)(rs)
    vx1[:,i] = interp.interp1d(radius_grid[:,i],vx1field[:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)
    vx2[:,i] = interp.interp1d(radius_grid[:,i],vx2field[:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)
    vx3[:,i] = interp.interp1d(radius_grid[:,i],vx3field[:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)

    #den[:,i] = interp.CubicSpline(radius_grid[:,i],denfield[:,i])(rs)
    #vx1[:,i] = interp.CubicSpline(radius_grid[:,i],vx1field[:,i])(rs)
    #vx2[:,i] = interp.CubicSpline(radius_grid[:,i],vx2field[:,i])(rs)
    #vx3[:,i] = interp.CubicSpline(radius_grid[:,i],vx3field[:,i])(rs) 

    if fill_copy:
      indices=np.array(range(den[:,i].size))
      
      # assume the edge is the same for all variables
      lowindex=np.min(np.where(den[:,i]>0,indices,indices.size))
      highindex=np.max(np.where(den[:,i]>0,indices,0))
 
      # check this selects the correct index
      den[:lowindex,i] = den[lowindex,i]
      den[highindex:,i] = den[highindex,i]
     
      vx1[:lowindex,i] = vx1[lowindex,i]
      vx1[highindex:,i] = vx1[highindex,i] 

      vx2[:lowindex,i] = vx2[lowindex,i]
      vx2[highindex:,i] = vx2[highindex,i]

      vx3[:lowindex,i] = vx3[lowindex,i]
      vx3[highindex:,i] = vx3[highindex,i]

    if isMHD:

      # interpolating on the Bfield leads to spurious non-zero div B so using the A field instead
  
      bx1[:,i] = interp.interp1d(radius_grid[:,i],bxfields[0][:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)    
      bx2[:,i] = interp.interp1d(radius_grid[:,i],bxfields[1][:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)
      bx3[:,i] = interp.interp1d(radius_grid[:,i],bxfields[2][:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)

      Ax1[:,i] = interp.interp1d(radius_grid[:,i],Afields[0][:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)    
      Ax2[:,i] = interp.interp1d(radius_grid[:,i],Afields[1][:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)
      Ax3[:,i] = interp.interp1d(radius_grid[:,i],Afields[2][:,i],kind='linear',bounds_error=False, fill_value=fill_value)(rs)

      if fill_copy:
        bx1[:lowindex,i] = bx1[lowindex,i]
        bx1[highindex:,i] = bx1[highindex,i]

        bx2[:lowindex,i] = bx2[lowindex,i]
        bx2[highindex:,i] = bx2[highindex,i]

        bx3[:lowindex,i] = bx3[lowindex,i]
        bx3[highindex:,i] = bx3[highindex,i]

        Ax1[:lowindex,i] = Ax1[lowindex,i]
        Ax1[highindex:,i] = Ax1[highindex,i]

        Ax2[:lowindex,i] = Ax2[lowindex,i]
        Ax2[highindex:,i] = Ax2[highindex,i]

        Ax3[:lowindex,i] = Ax3[lowindex,i]
        Ax3[highindex:,i] = Ax3[highindex,i]



  rho3D = np.zeros((Nx,Ny,Nz))
  vx3D = np.zeros((Nx,Ny,Nz))
  vy3D = np.zeros((Nx,Ny,Nz))
  vz3D = np.zeros((Nx,Ny,Nz))

  if isMHD:
    bx3D = np.zeros((Nx,Ny,Nz))
    by3D = np.zeros((Nx,Ny,Nz))
    bz3D = np.zeros((Nx,Ny,Nz))

    Ax3D = np.zeros((Nx,Ny,Nz))
    Ay3D = np.zeros((Nx,Ny,Nz))
    Az3D = np.zeros((Nx,Ny,Nz))

    # note - just doing unstratified/2.5D for now  
  for i in range(Nz):

    rho3D[:,:,i]=den
    vx3D[:,:,i]=vx1
    vy3D[:,:,i]=vx2
    vz3D[:,:,i]=vx3

    if isMHD:
      bx3D[:,:,i]=bx1
      by3D[:,:,i]=bx2
      bz3D[:,:,i]=bx3

      Ax3D[:,:,i]=Ax1
      Ay3D[:,:,i]=Ax2
      Az3D[:,:,i]=Ax3


  #rho = data[0]
  #vx = data[1]
  #vy = data[2]
  #vz = data[3]

  rho_rs =  np.transpose(rho3D, (2,1,0))
  vx_rs =  np.transpose(vx3D, (2,1,0))
  vy_rs =  np.transpose(vy3D, (2,1,0))
  vz_rs =  np.transpose(vz3D, (2,1,0))

  if isMHD:
    bx_rs = np.transpose(bx3D, (2,1,0))
    by_rs = np.transpose(by3D, (2,1,0))
    bz_rs = np.transpose(bz3D, (2,1,0))

    Ax_rs = np.transpose(Ax3D, (2,1,0))
    Ay_rs = np.transpose(Ay3D, (2,1,0))
    Az_rs = np.transpose(Az3D, (2,1,0))


   # this will also change for
   # non-isothermal discs
  if isMHD:
    Nfields=7
  else:
    Nfields=4  

  PLUTOCustomIC = np.zeros([Nfields, int(Nz), int(Ny), int(Nx)])

   # possibly we want to switch to labels?
  for n in range(0, Nfields):
    for i in range(0, int(Nz)):
      for j in range(0, int(Ny)):
        for k in range(0, int(Nx)):
          if (n==0):
            PLUTOCustomIC[n][i][j][k] = rho_rs[i][j][k]
          elif (n==1):
            PLUTOCustomIC[n][i][j][k] = vx_rs[i][j][k]
          elif (n==2):
            PLUTOCustomIC[n][i][j][k] = vy_rs[i][j][k]
          elif (n==3):
            PLUTOCustomIC[n][i][j][k] = vz_rs[i][j][k]
                #elif (n==4):
                #    PLUTOCustomIC[n][i][j][k] = prs_reshaped[i][j][k]
          elif (n==4):
            PLUTOCustomIC[n][i][j][k] = bx_rs[i][j][k]
            #PLUTOCustomIC[n][i][j][k] = Ax_rs[i][j][k]
          elif (n==5):
            #PLUTOCustomIC[n][i][j][k] = Ay_rs[i][j][k]
            PLUTOCustomIC[n][i][j][k] = by_rs[i][j][k]
          elif (n==6):
            #PLUTOCustomIC[n][i][j][k] = Az_rs[i][j][k]
            PLUTOCustomIC[n][i][j][k] = bz_rs[i][j][k]

  output_file = open(filename, 'wb')
  PLUTOCustomIC.tofile(output_file)
  output_file.close()

  if isMHD:

    if vectorPotential:
      return rho_rs, vx_rs, vy_rs, vz_rs, Ax_rs, Ay_rs, Az_rs
    else:
      print('Warning: Interpolation between grids does not preserve the solinoidal condition. This should only be used for a purely vertical field.') 
      return rho_rs, vx_rs, vy_rs, vz_rs, bx_rs, by_rs, bz_rs
  else:
    return rho_rs, vx_rs, vy_rs, vz_rs


