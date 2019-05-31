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
def helicity(alpha1,alpha2):
global B1 B2 C2 C2B2RAT CTHETA RC R
#Calculate constants required for fields and vector potential
   scsp=scipy.special
   tol=0.00001
   a1=alpha1*RC
   a2=abs(alpha2*RC)
   ar=abs(alpha2*R)
   # a1~0 and a2>0
   if abs(alpha1) <  tol and abs(alpha2) > tol 
      C2B2RAT=-scsp.jv(1,a2)./scsp.yv(1,a2)
      #normalise field on axis to give flux = 1
      B1=1./(pi*RC*(RC-pi*scsp.yv(1,a2).*R.*F1(ar)))

      #match field components at RC to give constants for field in outer layer
      B2=-pi*0.5*a2.*scsp.yv(1,a2).*B1
      C2=pi*0.5*a2.*scsp.jv(1,a2).*B1
      #continuity constants for vector potential to simplify calculations, K2 term
      CTHETA = 0.5*RC^2*B1
   
   elif   abs(alpha2) < tol: #a2~0
      #normalise field on axis to give flux = 1
      j1=scsp.jv(1,a1)
      j0=scsp.jv(0,a1)
      B1=1./(pi*(2*RC*j1./alpha1+j0*(R*R-RC*RC)))
      #Only use ctb1rat if a1 or a2 not zero
   else:
      ctb1rat=scsp.jv(1,a1).*(1./alpha1-1./alpha2)
      dis=2./(pi*a2)
   if alpha2 < 0: #a1 and a2 cont. Rc
      b2prod=-scsp.jv(1,a1).*scsp.yv(0,a2)-scsp.jv(0,a1).*scsp.yv(1,a2)
      C2B2RAT=(scsp.jv(0,a1).*scsp.jv(1,a2)+scsp.jv(1,a1).*scsp.jv(0,a2))
   else:
      b2prod=scsp.jv(1,a1).*scsp.yv(0,a2)-scsp.jv(0,a1).*scsp.yv(1,a2)
      C2B2RAT=(scsp.jv(0,a1).*scsp.jv(1,a2)-scsp.jv(1,a1).*scsp.jv(0,a2))
   
   C2B2RAT=C2B2RAT./b2prod

   #normalise field on axis to give flux = 1
   B1=1./(2*pi*(R*b2prod.*F1(ar)./(abs(alpha2).*dis)+RC*ctb1rat))

   #match field components at RC to give constants for field in outer layer
   B2=B1.*b2prod./dis
   C2=B2.*C2B2RAT

   #continuity constants for vector potential
   CTHETA = B1*RC.*ctb1rat
   

   #calculate terms for K
   if  abs(alpha1) < tol and abs(alpha2) < tol:
      helicity = 0
   else:
      helicity = k2(alpha2,B2,RC,R)
      
   elif abs(alpha2) < tol:
      helicity = k1(alpha1,B1,RC)+k2_vac(alpha1,B1,RC,R)
   else:
      helicity=k1(alpha1,B1,RC)+k2(alpha2,B2,RC,R)
   

   #calculates first term K1  in helicity (helicity of inner layer) function k1 
   def k1(alpha1,b1,rc):
      j0=scsp.jv(0,alpha1*rc)
      j1=scsp.jv(1,alpha1*rc)
      k1=2*pi*b1^2*(rc^2*(j0.^2+j1.^2)-2*rc*j0.*j1./alpha1)./alpha1
      return k1
#Calculates term k2 in helicity
#Works for postive and negative alpha2 - "sgn" is sign of alpha2 function k2 
   def k2(alpha2,b2,rc,r):
      #global  CTHETA
      xc=abs(alpha2*rc)
      x=abs(alpha2*r)
      f1c=F1(xc)
      f1=F1(x)
      f0c=F0(xc)
      f0=F0(x)
      sgn=alpha2/abs(alpha2)

      k2=r^2*(f0.^2+f1.^2)-2*r*f0.*f1./abs(alpha2)-rc^2*(f0c.^2+f1c.^2)
      k2=(k2+2*rc*f0c.*f1c./abs(alpha2))*2*pi*b2.^2./abs(alpha2)
      k2=k2+4*pi*b2.*CTHETA.*(f0c-f0)./abs(alpha2)
      k2=sgn*k2
      return k2

   #Calculates term k2 in helicity when alpha2 = 0 function k2_vac
   def k2_vac(alpha1,b1,rc,r):
      a1=alpha1*rc
      j1=scsp.jv(1,a1)
      j0=scsp.jv(0,a1)
      k2_vac=(2*rc^2*j1.^2./alpha1-rc^3*j0.*j1)*log(r/rc)+0.5*rc*j1.*j0*(r^2-rc^2)
      k2_vac=2*pi*b1.^2.*k2_vac
      return k2_vac
   #Function F representing normalised field in
   #outer layer of loop - F0,  function f=F0(x)
   def F0(x):
      #global C2B2RAT
      f=scsp.jv(0,x)+C2B2RAT.*scsp.yv(0,x)
      return f
   #Function F representing normalised field in
   #outer layer of loop - F0    function f=F1(x)
   def f=F1(x):
      #global C2B2RAT
      f=scsp.jv(1,x)+C2B2RAT.*scsp.yv(1,x)
      return f
   def f=f1x(x):
      f1 = x.*F1(x)
      return f




