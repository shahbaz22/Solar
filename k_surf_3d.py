import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import special
# This is a check for a git upload
R1=0.5
R2=1
a1=np.linspace(0.01,3,600) # Bessel functions only take positive values
a2=np.linspace(-4.01,4,600) #Assures that the arrays have the same number of elements
a1,a2=np.meshgrid(a1,a2)

scsp=scipy.special
def J(v,a,R):
	'Function to return bessel of 1st kind wrt order v,and R,a'
	BF11=scsp.jv(v,a*R)
	return (BF11)
def Y(v,a,R):
	'Function to return bessel of 2st kind wrt order v,and R,a'
	BF22=scsp.yv(v,a*R)
	return (BF22)

def C(a11,a22):
	C2=(J(0,a11,R1)*J(1,abs(a22),R1)+J(1,a11,R1)*J(0,abs(a22),R1))*(np.pi*abs(a22)*R1)*(1/2)
	return(C2)

def B(a11,a22):
	B2=(-(J(1,a11,R1)*Y(0,abs(a22),R1)+J(0,a11,R1)*Y(1,abs(a22),R1))*(np.pi*abs(a22)*R1)*(1/2))
	return(B2)

#Delta=2/(np.pi*abs(a2)*R1) 

def F(v,a,R):
	F0a2R1=J(v,a,R)+(C(a1,a2)/B(a1,a2))*Y(v,a,R)
	return (F0a2R1)

K1=(2*np.pi/a1)*((R1**2)*(J(0,a1,R1)**2 + J(1,a1,R1)**2)
	-2*R1*(1/a1)*J(0,a1,R1)*J(1,a1,R1))
K2=((2*np.pi*B(a1,a2)**2)/abs(a2))*(((F(0,abs(a2),R2)**2) +(F(1,abs(a2),R2)**2) -
	(2*R2/abs(a2))*F(0,abs(a2),R2)*F(1,abs(a2),R2)))
K3=((2*np.pi*B(a1,a2)**2)/abs(a2))*((R1**2)*((F(0,abs(a2),R1)**2) +(F(1,abs(a2),R1)**2)- 
	(2*R1/abs(a2))*F(0,abs(a2),R1)*F(1,abs(a2),R1)))
K4=(4*np.pi*B(a1,a2)/abs(a2))*(F(0,abs(a2),R1)-F(0,abs(a2),R2))*R1*J(1,a1,R1)*(1/a1 - 1/a2)


norm=((2*np.pi*B(a1,a2))/abs(a2))*R2*F(1,abs(a2),R2)+2*np.pi*R1*J(1,a1,R1)*(1/a1 - 1/a2)

K=+K1+K2-K3+K4
kn=K/(norm**2)

fig=plt.figure()
ax=fig.gca(projection='3d')
surf=ax.plot_surface(a1,a2,kn, cmap=cm.coolwarm,
	linewidth=0,antialiased=False)
ax.set_xlabel('a1')
ax.set_ylabel('a2')
ax.set_zlabel('K')
#Add color bar
fig.colorbar(surf,shrink=0.5,aspect=5) 
ax.set_zlim(0, 2.5)
#plt.ylim(-4,4)
#plt.xlim(3,0)
plt.show()