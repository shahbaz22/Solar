'not working program to plot psi0 vs a1 for a2=0'
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
'alpha1 values'
al1=np.linspace(0.1,3,800) # Bessel functions only take positive values
scsp=scipy.special
def j(v,a):
	'Function to return bessel of 1st kind wrt order v,and R,a'
	BF11=scsp.jv(v,a)
	return (BF11)
def y(v,a):
	'Function to return bessel of 2st kind wrt order v,and R,a'
	BF22=scsp.yv(v,a)
	return (BF22)
'ps0 function for a2=0'
psi0=(2*np.pi*R1/al1)*j(1,al1*R1)
psi0=psi0+np.pi*(R2**2-R1**2)*j(0,al1*R1)
plt.plot(al1,psi0)
plt.title('psi0 vs a1')
plt.xlabel('a1')
plt.ylabel('psi0')
plt.savefig('psi0 for a2=0')
plt.show()