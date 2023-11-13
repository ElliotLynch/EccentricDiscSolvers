import numpy as np


#will do it dynamically later
NMAXITER=1000
TOL=2.2e-16

 # SpecialIp = (1/(2 Pi)) int (1 - q cos x)^(-p) d x
def SpecialIp(q,p):

   # can define terms recursively
  Ip=0.0
  Inp=1.0
  for n in range(1,NMAXITER):
    Ip+=Inp
    Inp *= 0.25*(p + 2.0*n - 1)*(p + 2.0*n - 2.0)*q*q/(n*n)
    
    if Inp/Ip < TOL:
      break

  return Ip

 # SpecialIp1 = (1/(2 Pi)) int cos x (1 - q cos x)^(-p) d x
def SpecialIp1(q,p):

   # can define terms recursively
  Ip1=0.0
  Inp1=0.5*p
  for n in range(1,NMAXITER):
    Ip1+=Inp1
    Inp1 *= 0.25*(p + 2*n)*(p + 2*n - 1)*q*q/((n+1)*n)

    if Inp1/Ip1 < TOL:
      break

  return q*Ip1  

 # d SpecialIp/d q
def DSpecialIp(q,p):

   # can define terms recursively
  dIp=0.5*(p+1)*p
  dInp=0.25*(p+1)*p
  for n in range(2,NMAXITER):
    #Ip+=2*n*Inp
    dInp *= 0.25*(p + 2.0*n - 1)*(p + 2.0*n - 2.0)*q*q/(n*n)
    dIp+=2*n*dInp

    if dInp/dIp < TOL:
      break

  return dIp*q


 # d SpecialIp1/d q
 # double check this
def DSpecialIp1(q,p):

   # can define terms recursively
  dIp1=0.0
  dInp1=0.5*p
  for n in range(1,NMAXITER):
    dIp1+=(2*n-1)*dInp1
    dInp1 *= 0.25*(p + 2*n)*(p + 2*n - 1)*q*q/((n+1)*n)

    if dInp1/dIp1 < TOL:
      break

  return dIp1



 # d^2 SpecialIp/d q^2
def DDSpecialIp(q,p):

   # can define terms recursively
  ddIp=0.5*(p+1)*p
  ddInp=0.25*(p+1)*p
  for n in range(2,NMAXITER):
    #Ip+=2*n*Inp
    ddInp *= 0.25*(p + 2.0*n - 1)*(p + 2.0*n - 2.0)*q*q/(n*n)
    ddIp+=2*n*(2*n-1)*ddInp

    if ddInp/ddIp < TOL:
      break

  return ddIp


 # d SpecialIp/d q
def DDSpecialIp1(q,p):

   # can define terms recursively
  dInp1=0.5*p
  dIp1=0.0
  dInp1 *= 0.25*(p + 2)*(p + 1)/2
  #dIp1=2.0*dInp1
  for n in range(2,NMAXITER):
    dIp1+=(2*n-1)*(2*n-2)*dInp1
    dInp1 *= 0.25*(p + 2*n)*(p + 2*n - 1)*q*q/((n+1)*n)
    

    if dInp1/dIp1 < TOL:
      break

  return dIp1*q



if __name__ == '__main__':
  import matplotlib.pyplot as plt

  qs=np.linspace(0.0,1.0,100)

  p=5.0/3.0

  Ips=np.array([SpecialIp(q,p) for q in qs])

  plt.semilogy(qs,Ips,'k-')
  plt.show()

  Im1p=np.array([SpecialIp(q,-1) for q in qs])
  I0p=np.array([SpecialIp(q,0) for q in qs])
  I1p=np.array([SpecialIp(q,1) for q in qs]) 
  I2p=np.array([SpecialIp(q,2) for q in qs])


  plt.semilogy(qs,Im1p,'k-')
  plt.semilogy(qs,I0p,'r-')
  plt.semilogy(qs,I1p,'b-')
  plt.semilogy(qs,I2p,'g-')

  plt.semilogy(qs,np.ones_like(qs),'kx--')
  plt.semilogy(qs,np.ones_like(qs),'rx--')
  plt.semilogy(qs,1.0/np.sqrt(1.0 - qs*qs),'bx--')
  plt.semilogy(qs,(1.0 - qs*qs)**(-1.5),'gx--')

  plt.show()

  DIm1p=np.array([DSpecialIp(q,-1) for q in qs])
  DI0p=np.array([DSpecialIp(q,0) for q in qs])
  DI1p=np.array([DSpecialIp(q,1) for q in qs]) 
  DI2p=np.array([DSpecialIp(q,2) for q in qs])

  print(DIm1p,DI0p)

  plt.semilogy(qs,DIm1p,'k-')
  plt.semilogy(qs,DI0p,'r-')
  plt.semilogy(qs,DI1p,'b-')
  plt.semilogy(qs,DI2p,'g-')

  plt.semilogy(qs,np.zeros_like(qs),'kx--')
  plt.semilogy(qs,np.zeros_like(qs),'rx--')
  plt.semilogy(qs,qs*((1.0 - qs*qs)**(-1.5)),'bx--')
  plt.semilogy(qs,3.0*qs*((1.0 - qs*qs)**(-2.5)),'gx--')
  
  plt.show()
 
  DDIm1p=np.array([DDSpecialIp(q,-1) for q in qs])
  DDI0p=np.array([DDSpecialIp(q,0) for q in qs])
  DDI1p=np.array([DDSpecialIp(q,1) for q in qs])
  DDI2p=np.array([DDSpecialIp(q,2) for q in qs])

  print(DDIm1p,DDI0p)

  plt.semilogy(qs,DDIm1p,'k-')
  plt.semilogy(qs,DDI0p,'r-')
  plt.semilogy(qs,DDI1p,'b-')
  plt.semilogy(qs,DDI2p,'g-')
 
  plt.semilogy(qs,np.zeros_like(qs),'kx--')
  plt.semilogy(qs,np.zeros_like(qs),'rx--')
  plt.semilogy(qs,((1.0 - qs*qs)**(-1.5)) + 3.0*qs*qs*((1.0 - qs*qs)**(-2.5)),'bx--')
  plt.semilogy(qs,3.0*(1.0 - qs*qs)**(-2.5) + 5.0*3.0*qs*qs*(1.0 - qs*qs)**(-3.5)  ,'gx--')
  
  plt.show()

  Im1p1=np.array([SpecialIp1(q,-1) for q in qs])
  I0p1=np.array([SpecialIp1(q,0) for q in qs])
  I1p1=np.array([SpecialIp1(q,1) for q in qs])
  I2p1=np.array([SpecialIp1(q,2) for q in qs])

  plt.plot(qs,Im1p1,'k-')
  plt.plot(qs,I0p1,'r-') 
  plt.plot(qs,I1p1,'b-')
  plt.plot(qs,I2p1,'g-')
   
  plt.plot(qs,-0.5*qs,'kx--')
  plt.plot(qs,(I0p - Im1p)/qs,'rx--')
  plt.plot(qs,(I1p - I0p)/qs,'bx--')
  plt.plot(qs,(I2p - I1p)/qs,'gx--')

  plt.show()


  DIm1p1=np.array([DSpecialIp1(q,-1) for q in qs])
  DI0p1=np.array([DSpecialIp1(q,0) for q in qs])
  DI1p1=np.array([DSpecialIp1(q,1) for q in qs])
  DI2p1=np.array([DSpecialIp1(q,2) for q in qs])

  plt.plot(qs,DIm1p1,'k-')
  plt.plot(qs,DI0p1,'r-')
  plt.plot(qs,DI1p1,'b-')
  plt.plot(qs,DI2p1,'g-')
   
  plt.plot(qs,-0.5*np.ones_like(qs),'kx--')
  plt.plot(qs,(DI0p - DIm1p)/qs - (I0p - Im1p)/(qs*qs),'rx--')
  plt.plot(qs,(DI1p - DI0p)/qs  - (I1p - I0p)/(qs*qs),'bx--')
  plt.plot(qs,(DI2p - DI1p)/qs  - (I2p - I1p)/(qs*qs),'gx--')

  plt.show()

  DDIm1p1=np.array([DDSpecialIp1(q,-1) for q in qs])
  DDI0p1=np.array([DDSpecialIp1(q,0) for q in qs])
  DDI1p1=np.array([DDSpecialIp1(q,1) for q in qs])
  DDI2p1=np.array([DDSpecialIp1(q,2) for q in qs])
      
  plt.plot(qs,DDIm1p1,'k-')
  plt.plot(qs,DDI0p1,'r-')
  plt.plot(qs,DDI1p1,'b-')
  plt.plot(qs,DDI2p1,'g-')
   
  plt.plot(qs,np.zeros_like(qs),'kx--')
  plt.plot(qs,(DDI0p - DDIm1p)/qs - 2*(DI0p - DIm1p)/(qs*qs) + 2*(I0p - Im1p)/(qs*qs*qs),'rx--')
  plt.plot(qs,(DDI1p - DDI0p)/qs  - 2*(DI1p - DI0p)/(qs*qs)  + 2*(I1p - I0p)/(qs*qs*qs),'bx--')
  plt.plot(qs,(DDI2p - DDI1p)/qs  - 2*(DI2p - DI1p)/(qs*qs)  + 2*(I2p - I1p)/(qs*qs*qs),'gx--')
  
  plt.show()





