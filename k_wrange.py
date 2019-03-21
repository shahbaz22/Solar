import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import special

R1=0.5
R2=1
a1=np.linspace(0.01,4,600) # Bessel functions only take positive values
a2=np.linspace(0.01,4,600) #Assures that the arrays have the same number of elements
#a1,a2=np.meshgrid(a1,a2)
#a1 = np.asarray(np.meshgrid(a1,a2))
#a2 = a1
#print('type(a2): ', type(a2))
#print('type(a1): ', type(a1))
#print(a1.shape)
#print(a2.shape)
#print('a1[4,4]: ', a1[4,4])
#print('a2[4,4]: ', a2[4,4])
#plt.scatter(a1,a2)
#plt.show()

scsp=scipy.special
def J(v,a,R):
	'Function to return bessel of 1st kind wrt order v,and R,a'
	BF11=scsp.jv(v,a*R)
	return (BF11)
def Y(v,a,R):
	'Function to return bessel of 2st kind wrt order v,and R,a'
	BF22=scsp.yv(v,a*R)
	return (BF22)

Delta=2/(np.pi*abs(a2)*R1) 

#np.multiply element wise multiplication
#np.dot
C2=J(0,a1,R1)*J(1,abs(a2),R1)+J(1,a1,R1)*J(0,abs(a2),R1)/Delta 
B2=-(J(1,a1,R1)*Y(0,abs(a2),R1)+J(0,a1,R1)*Y(1,abs(a2),R1)/Delta)

F0a2R1=J(0,abs(a2),R1)+(C2/B2)*Y(0,abs(a2),R1)
F1a2R1=J(1,abs(a2),R1)+(C2/B2)*Y(1,abs(a2),R1)
F0a2R2=J(0,abs(a2),R2)+(C2/B2)*Y(0,abs(a2),R2)
F1a2R2=J(1,abs(a2),R2)+(C2/B2)*Y(1,abs(a2),R2)


K1=(2*np.pi/a1)*((R1**2)*(J(0,a1,R1)**2 + J(1,a1,R1)**2)
	-2*R1*(1/a1)*J(0,a1,R1)*J(1,a1,R1))
K2=((2*np.pi*B2**2)/a2)*(((F0a2R2**2) +(F1a2R2**2) -
	(2*R2/a2)*F0a2R2*F1a2R2))
K3=((2*np.pi*B2**2)/a2)*((R1**2)*((F0a2R1**2) +(F1a2R1**2)- 
	(2*R1/a2)*F0a2R1*F1a2R1))
K4=(4*np.pi*B2/a2)*(F0a2R1-F0a2R2)*R1*J(1,a1,R1)*(1/a1 - 1/a2)


norm=((2*np.pi*B2)/a2)*R2*F1a2R2+2*np.pi*R1*J(1,a1,R1)*(1/a1 - 1/a2)

#W1=((J(0,a1,R1)**2 + J(1,a1,R1)**2)
#	-(1/a1)*J(0,a1,R1)*J(1,abs(a2),R1))
#W2= ((B2**2)*R2**2)*((F0a2R2**2 +F1a2R2**2 -
#	(R2/a2)*F0a2R2*F1a2R2))
#W3=((F0a2R1**2 +F1a2R1**2)+F0a2R1*F1a2R1*(1/a2))

K=+K1+K2-K3+K4
Knorm=K/(norm**2)
#W=W1+W2+W3

#fig=plt.figure()
#ax=fig.gca(projection='3d')
#surf=ax.plot_surface(a1,a2,K, cmap=cm.coolwarm,
#	linewidth=0,antialiased=False)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#surf=ax.plot_surface(a1,a2,W,cmap=cm.coolwarm,
	#linewidth=0,antialiased=False)
#ax.set_xlabel('a1')
#ax.set_ylabel('a2')
#ax.set_zlabel('K')
#Add color bar
#fig.colorbar(surf,shrink=0.5,aspect=5) 
#ax.set_zlim(0, 2.5)
#plt.ylim(-4,4)
#plt.xlim(3,0)
plt.plot(a1,K)
plt.show()
