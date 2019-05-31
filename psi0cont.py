'working porgram for a1,a2>0, for psi0 with surface plot'
import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.cm as cm
import matplotlib.colors as cc
from scipy import special
R1=0.5
R2=1
'a1 already used in kn function'
al1=np.linspace(0.1,3,800) # Bessel functions only take positive values
al2=np.linspace(0.1,4,800) #Assures that the arrays have the same number of elements
alpha1,alpha2=np.meshgrid(al1,al2)
scsp=scipy.special
'Function to take a1 and a2 vales can return 2D K'
'Used for making surface plot of a1,a2,K'
'check if def kn can take array'
def psi0(alpha1,alpha2):
    
    'Variables for alpha and R values'
    'Simplified meshgid vals (a1R1==a1,a2==a2R1,a2R2==ar in PB code)'     
    a1=alpha1*R1
    a2=abs(alpha2)*R1
    ar=abs(alpha2*R2)
    def j(v,a):
    	'Function to return bessel of 1st kind wrt order v,and R,a'
    	BF11=scsp.jv(v,a)
    	return (BF11)
    def y(v,a):
    	'Function to return bessel of 2st kind wrt order v,and R,a'
    	BF22=scsp.yv(v,a)
    	return (BF22)
    'Now need to make array for each constant over certain range'
    'Ranges needed a2>0,a2<0,a1 always >0'

    'for a2>0'
    b2prod=j(1,a1)*y(0,a2)-j(0,a1)*y(1,a2)
    c2prod=j(0,a1)*j(1,a2)-j(1,a1)*j(0,a2)
    ctb1rat=j(1,a1)*(1/alpha1-1/alpha2)
    dis=2/(np.pi*a2)

    'ratio of terms needd to calculate F, B1 cancels'
    c2b2rat=c2prod/b2prod 
    def F0(x):
    	f=j(0,x)+c2b2rat*y(0,x)
    	return f
    def F1(x):
    	f=j(1,x)+c2b2rat*y(1,x)
    	return f
    'norm=1 to find B1, then using B1 find C2, B2'
    psi0=(2*np.pi*(R2*b2prod*(F1(ar)/(abs(alpha2)*dis))+R1*ctb1rat))   
    return psi0

cp = plt.contourf(alpha1, alpha2, psi0(alpha1,alpha2))
plt.colorbar(cp)
plt.title('Psi0 contour plot')
plt.xlabel('a1')
plt.ylabel('a2')
plt.savefig('Psi0cont_a1_a2.png')
plt.show()
