import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import special
#calculates helicity expression with total flux normalised to 1
#Uses Finn and Antonsen ("FA") formula involving psi and I
#Includes cases alpha2 > and < 0 also alpha2 = 0 and alpha1 = 0 
R1=0.5
R2=1
a1=np.linspace(0.01,4,600) # Bessel functions only take positive values
a2=np.linspace(-4.01,4,600) #Assures that the arrays have the same number of elements
a1,a2=np.meshgrid(a1,a2)

def helicity(alpha1,alpha2):
   scsp=scipy.special
   tol=0.00001
   R1=0.5
   R2=1
   a1=alpha1*R1
   a2=abs(alpha2*R1)
   ar=abs(alpha2*R2)

   #used for finding B1 for normalisation
   if alpha1==0 and alpha2==0:
     return('helicity is',0)

   ctb1rat=scsp.jv(1,a1)*(1/alpha1-1/alpha2)
   dis=2/(np.pi*a2)

   if alpha2 < 0:
      b2prod=-scsp.jv(1,a1)*scsp.yv(0,a2)-scsp.jv(0,a1)*scsp.yv(1,a2)
      C2B2RAT=(scsp.jv(0,a1)*scsp.jv(1,a2)+scsp.jv(1,a1)*scsp.jv(0,a2))
   else:
      b2prod=scsp.jv(1,a1)*scsp.yv(0,a2)-scsp.jv(0,a1)*scsp.yv(1,a2)
      C2B2RAT=(scsp.jv(0,a1)*scsp.jv(1,a2)-scsp.jv(1,a1)*scsp.jv(0,a2))
   #Ratio from if statments which cancels value for B1

   C2B2RAT=C2B2RAT/b2prod 

   def F0(x):
      f=scsp.jv(0,x)+C2B2RAT*scsp.yv(0,x)
      return f
   #Function F representing normalised field in
   #outer layer of loop - F0    function f=F1(x)
   def F1(x):
      #global C2B2RAT
      f=scsp.jv(1,x)+C2B2RAT*scsp.yv(1,x)
      return f


   #normalise field on axis to give flux = 1 from if statments abs val of a2 used
   B1=1/(2*np.pi*(R2*b2prod*F1(ar)/(abs(alpha2)*dis)+R1*ctb1rat))

   #match field components at R1 to give constants for field in outer layer
   B2=B1*b2prod/dis
   C2=B2*C2B2RAT

   #continuity constants for vector potential to simplify calculations, K2 term
   CTHETA = B1*R1*ctb1rat

   #calculates first term K1  in helicity (helicity of inner layer) function k1 
   def k1(alpha1,b1,r1):
      j0=scsp.jv(0,alpha1*r1)
      j1=scsp.jv(1,alpha1*r1)
      k1=(2*np.pi*b1**2)*((r1**2)*(j0**2+j1**2)-2*r1*j0*j1/alpha1)/alpha1
      return k1
#Calculates term k2 in helicity
#Works for postive and negative alpha2 - "sgn" is sign of alpha2 function k2 
   def k2(alpha2,b2,r1,r2):
      #global  CTHETA
      xc=abs(alpha2*r1)
      x=abs(alpha2*r2)
      f1c=F1(xc)
      f1=F1(x)
      f0c=F0(xc)
      f0=F0(x)
      sgn=alpha2/abs(alpha2)

      k2=(r2**2)*(f0**2+f1**2)-2*r2*f0*f1/abs(alpha2)-(r1**2)*(f0c**2+f1c**2)
      k2=(k2+2*r1*f0c*f1c/abs(alpha2))*(2*np.pi*b2**2)/abs(alpha2)
      k2=k2+4*np.pi*b2*CTHETA*(f0c-f0)/abs(alpha2)
      k2=sgn*k2 # apart from K1 the other part of Ks sign depends on the sign of a2
      return k2
   helicity=k1(alpha1,B1,R1)+k2(alpha2,B2,R1,R2) #main case needed
   return helicity
   #Calculates term k2 in helicity when alpha2 = 0 function k2_vac
   '''def k2_vac(alpha1,b1,r1,r):
      a1=alpha1*r1
      j1=scsp.jv(1,a1)
      j0=scsp.jv(0,a1)
      k2_vac=(2*r1**2*j1.**2/alpha1-r1**3*j0*j1)*log(r/r1)+0.5*r1*j1*j0*(r**2-r1**2)
      k2_vac=2*np.pi*b1.**2*k2_vac
      return k2_vac'''
   #Function F representing normalised field in
   #outer layer of loop - F0,  function f=F0(x)

   #calculate terms for K
 

print(helicity(a1,a2))



