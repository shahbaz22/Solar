import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import special

R1=np.linspace(0,0.5,1000)
R2=np.linspace(0.5,1,1000) # like this both R1 and R2 are the same size
a1=4
a2=-2

scsp=scipy.special
def J(v,a,R):
	'Function to return bessel of 1st kind wrt order v,and R,a'
	BF11=scsp.jv(v,a*R)
	return (BF11)
def Y(v,a,R):
	'Function to return bessel of 2st kind wrt order v,and R,a'
	BF22=scsp.yv(v,a*R)
	return (BF22)

Delta=2/(np.pi*abs(a2)*0.5) 
#Only Rc used set to 0.5
C2=(J(0,a1,0.5)*J(1,abs(a2),0.5)+J(1,a1,0.5)*J(0,abs(a2),0.5))/Delta 
B2=-(J(1,a1,0.5)*Y(0,abs(a2),0.5)+J(0,a1,0.5)*Y(1,abs(a2),0.5))/Delta
BzR1=J(0,a1,R1)
BtR1=J(1,a1,R1)
#xrange on Bt just allows the bessel function to vary

BzR2=B2*J(0,abs(a2),R2)+C2*Y(0,abs(a2),R2)
BtR2=-(B2*J(1,abs(a2),R2)+C2*Y(1,abs(a2),R2))

plt.plot(R1,BzR1,R2,BzR2,color='blue',label='Bz')
plt.plot(R1,BtR1,R2,BtR2,color='green',label='Bt',linestyle='--')
plt.axvline(x=0.5,color='red')
plt.xlabel('r')
plt.ylabel('B')
plt.title('a1=4  ,a2=-2')
plt.legend()
plt.show()
