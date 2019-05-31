'working program for k with a2=0 (ksingle) 2layer'
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import special
R=0.5
a=np.linspace(0.1,3,800) # Bessel functions only take positive values
scsp=scipy.special
def helicity(a):
    'function to return k for constant alpha in all layers'
    ar=a*R
    def j(v,a):
    	'Function to return bessel of 1st kind wrt order v,and R,a'
    	BF11=scsp.jv(v,a)
    	return (BF11)
    def y(v,a):
    	'Function to return bessel of 2st kind wrt order v,and R,a'
    	BF22=scsp.yv(v,a)
    	return (BF22)

    psi0=(2*np.pi*R/a)*j(1,ar)+np.pi*(1-R**2)*j(0,ar)
    r=1
    k=(2*np.pi/a)*((R**2)*(j(0,ar)**2)+(R**2)*(j(1,ar)**2)-2*(R/a)*j(0,ar)*j(1,ar))
    k=k+(2*(R**2)*((j(1,ar)**2)/a)-(R**3)*(j(0,ar)**2)*(j(1,ar)**2))*np.log(r/R)+0.5*R*j(0,ar)*j(1,ar)*(r**2-R**2)
    k=2*np.pi*k
    return k#/(psi0**2)
'for k_single check'
plt.plot(a,helicity(a))
plt.xlabel('a')
plt.ylabel('K')
#plt.savefa2r2('Ktest.png') 
plt.show()