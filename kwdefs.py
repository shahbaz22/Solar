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
a1=1.5
a2=-2

#a1=np.linspace(0.01,4,600) # Bessel functions only take positive values
#a2=np.linspace(-4.01,4,600) #Assures that the arrays have the same number of elements

#a1,a2=np.meshgrid(a1,a2)

scsp=scipy.special
def J(v,a,R):
	'Function to return bessel of 1st kind wrt order v,and R,a'
	BF11=scsp.jv(v,a*R)
	return (BF11)
def Y(v,a,R):
	'Function to return bessel of 2st kind wrt order v,and R,a'
	BF22=scsp.yv(v,a*R)
	return (BF22)
	


#Delta=2/(np.pi*abs(a2)*R1) 

def kn(a1,a2):
#a1>0,a2<0;a1>0,a2>0
	#norm0=((2*np.pi*B2)/abs(a2))*R2*F(1,abs(a2),R2)+2*np.pi*R1*J(1,a1,R1)*(1/a1 - 1/a2)
#a1>0,a2=0
	#norm1=(2*np.pi*R1/a1)*J(1,a1,R1)+np.pi*(R2**2-R1**2)*J(0,a1,R1)
#a1=0,a2>0
	#norm2=np.pi*(R1**2-np.pi*R1*R2*Y(1,a2,R1)*F(1,a2,R1))
	C2=(J(0,a1,R1)*J(1,abs(a2),R1)+J(1,a1,R1)*J(0,abs(a2),R1))*(np.pi*abs(a2)*R1)*(1/2)
	B2=(-(J(1,a1,R1)*Y(0,abs(a2),R1)+J(0,a1,R1)*Y(1,abs(a2),R1))*(np.pi*abs(a2)*R1)*(1/2))
	
	def F(v,a,R):
		F0a2R1=J(v,a,R)+(C2/B2)*Y(v,a,R)
		return (F0a2R1)

	K1=(2*np.pi/a1)*((R1**2)*(J(0,a1,R1)**2 + J(1,a1,R1)**2)
	-2*R1*(1/a1)*J(0,a1,R1)*J(1,a1,R1))
	K2=((2*np.pi*B2**2)/abs(a2))*(((F(0,abs(a2),R2)**2) +(F(1,abs(a2),R2)**2) -
	(2*R2/abs(a2))*F(0,abs(a2),R2)*F(1,abs(a2),R2)))
	K3=((2*np.pi*B2**2)/abs(a2))*((R1**2)*((F(0,abs(a2),R1)**2) +(F(1,abs(a2),R1)**2)- 
	(2*R1/abs(a2))*F(0,abs(a2),R1)*F(1,abs(a2),R1)))
	K4=(4*np.pi*B2/abs(a2))*(F(0,abs(a2),R1)-F(0,abs(a2),R2))*R1*J(1,a1,R1)*(1/a1 - 1/a2)

	def cond1(K1,K2,K3,K4):
		K=+K1-K2+K3-K4
		norm0=((2*np.pi*B2)/abs(a2))*R2*F(1,abs(a2),R2)+2*np.pi*R1*J(1,a1,R1)*(1/a1 - 1/a2)
		kn=K/(norm0**2)
		return kn

	def cond2(K1):
		K=K1+(2*np.pi/a1)*(((2*R1**2)/a1)*(J(1,a1,R1)**2)*np.log(R2/R1) 
		-(R1**3)*J(0,a1,R1)*J(1,a1,R1)*np.log(R2/R1)
		+0.5*R1*J(0,a1,R1)*J(1,a1,R1)*(R1**2-R2**2))
		norm1=(2*np.pi*R1/a1)*J(1,a1,R1)+np.pi*(R2**2-R1**2)*J(0,a1,R1)
		kn=K/(norm1**2)
		return kn

	def cond3(K2,K3):
		K=K2-K3+((2*np.pi*B2*R1**2)/a2)*(F(0,a2,R1)-F(0,a2,R2))
		norm2=np.pi*(R1**2-np.pi*R1*R2*Y(1,a2,R1)*F(1,a2,R1))
		kn=K/(norm2**2)
		return kn
	
	kn = np.zeros(np.shape(a1))
	
	kn = np.where((a1>0) & (a2<0), cond1(K1,K2,K3,K4),kn)
	kn = np.where((a1>0) & (a2>0), cond1(K1,-1*K2,-1*K3,-1*K4), kn)
	kn = np.where((a1>0) & (a2==0), cond2(K1), kn)
	kn = np.where((a1==0) & (a2>0), cond3(K2,K3), kn)
	return(kn)

	#ind1 = a1>0 and a2<0
	#kn[ind1] = cond1(K1[ind1],K2[ind1],K3[ind1],K4[ind1],norm0)
	#ind2 = a1>0 and a2>0
	#kn[ind2] = cond1(K1,-1*K2,-1*K3,-1*K4,norm0)

	"""
	if np.any(a1>0) and np.any(a2<0):
		K=+K1-K2+K3-K4
		kn=K/(norm0**2)
	elif np.any(a1>0) and np.any(a2>0):
		K=K1+K2-K3+K4
		kn=K/(norm0**2)
	elif np.any(a1>0) and np.any(a2==0):
		K=K1+(2*np.pi/a1)*(((2*R1**2)/a1)*(J(1,a1,R1)**2)*np.log(R2/R1) 
		-(R1**3)*J(0,a1,R1)*J(1,a1,R1)*np.log(R2/R1)
		+0.5*R1*J(0,a1,R1)*J(1,a1,R1)(R1**2-R2**2))
		kn=K/(norm1**2)
	elif np.any(a1==0) and np.any(a2>0):
		K=K2-K3+((2*np.pi*B2*R1**2)/a2)*(F(0,a2,R1)-F(0,a2,))
		kn=K/(norm2**2)
	"""
	return(kn)


'''fig=plt.figure()
ax=fig.gca(projection='3d')
surf=ax.plot_surface(a1,a2,kn(a1,a2), cmap=cm.coolwarm,
	linewidth=0,antialiased=False)
ax.set_xlabel('a1')
ax.set_ylabel('a2')
ax.set_zlabel('K')
#Add color bar
fig.colorbar(surf,shrink=0.5,aspect=5) 
ax.set_zlim(0, 2.5)
plt.show()
print(np.min(kn))'''
print(kn(a1,a2))