'working porgram for a1,a2>0, for K with surface plot'
'with simplified terms, for trouble shooting'
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
'a1r1 already used in kn function'
al1=np.linspace(0.1,3,800) # Bessel functions only take positive values
al2=np.linspace(0.1,4,800) #Assures that the arrays have the same number of elements
alpha1,alpha2=np.meshgrid(al1,al2)
scsp=scipy.special
'Function to take a1r1 and a1r1 vales can return 2D K'
'Used for making surface pla2r2 of a1r1,a1r1,K'
'check if def kn can take ara2r2y'
def kn(alpha1,alpha2):
    'Variables for alpha and R values'
    'Simplified meshgid vals (a1R1==a1r1,a1r1==a1r1R1,a1r1R2==ar in PB code)'     
    a1r1=alpha1*R1
    a2r1=abs(alpha2)*R1
    a2r2=abs(alpha2*R2)
    def j(v,a):
    	'Function to return bessel of 1st kind wrt order v,and R,a'
    	BF11=scsp.jv(v,a)
    	return (BF11)
    def y(v,a):
    	'Function to return bessel of 2st kind wrt order v,and R,a'
    	BF22=scsp.yv(v,a)
    	return (BF22)
    'Now need to make array for each constant over certain range'
    'for a2>0'
    c2prod=j(0,a1r1)*j(1,a2r1)-j(1,a1r1)*j(0,a2r1)
    b2prod=j(1,a1r1)*y(0,a2r1)-j(0,a1r1)*y(1,a2r1)
    ctb1rat=j(1,a1r1)*(1/alpha1-1/alpha2)
    dis=2/(np.pi*a2r1)

    'ratio of terms needd to calculate F, B1 cancels'
    c2b2rat=c2prod/b2prod 
    def F0(x):
    	f=j(0,x)+c2b2rat*y(0,x)
    	return f
    def F1(x):
    	f=j(1,x)+c2b2rat*y(1,x)
    	return f
    'norm=1 to find B1, then using B1 find C2, B2'
    
    normt=(2*np.pi*(R2*b2prod*(F1(a2r2)/(abs(alpha2)*dis))+R1*ctb1rat))   
    B1=1/normt
    B2=B1*b2prod/dis
    C2=B1*c2prod/dis

 
    'calculates first term K1  in helicity (helicity of inner layer) function k1'
    'Fist term only uses R1,(Rc)'
    def k1(alpha1,r1):
        j0=j(0,alpha1*r1)
        j1=j(1,alpha1*r1)
        k1=(2*np.pi*B1**2)*((r1**2)*(j0**2+j1**2)-2*r1*j0*j1/alpha1)/alpha1
        return k1
#Calculates term k2 in helicity
#Works for postive and negative alpha1r1 - "sgn" is sign of alpha1r1 function k2 
    def k2(alpha2,r1,r2):
        'k2 where all terms are functions of a2'
        f1r1=F1(abs(alpha2*r1))
        f1r2=F1(abs(alpha2*r2))
        f0r1=F0(abs(alpha2*r1))
        f0r2=F0(abs(alpha2*r2))
        sgn=alpha2/abs(alpha2)

        k2=((r2**2)*(f0r2**2+f1r2**2)-2*r2*f0r2*f1r2/abs(alpha2))-(r1**2)*(f0r1**2+f1r1**2)
        k2=(k2+2*r1*f0r1*f1r1/abs(alpha2))*(2*np.pi*B2**2)/abs(alpha2) #here whole term multipied by constant
        k2=k2+(4*np.pi*B2*B1*r1*ctb1rat*(f0r1-f0r2))/abs(alpha2)
        k2=sgn*k2 # apart from K1 the other part of Ks siga2r2depends on the sign of a1r1
        return k2

    helicity=k1(alpha1,R1)+k2(alpha2,R1,R2) #main case needed
    return helicity

fig=plt.figure()
ax=fig.gca(projection='3d')
fn = kn(alpha1,alpha2)
zlim = 2.5
norm = cc.LogNorm(vmin=0.1, vmax=zlim)
surf=ax.plot_surface(alpha1,alpha2, fn, cmap=cm.hot,
                     linewidth=0,antialiased=False,
                     norm=norm)
ax.set_xlabel('a1')
ax.set_ylabel('a2')
ax.set_zlabel('K')
#Add color bar
fig.colorbar(surf,shrink=0.5,aspect=5) 
ax.set_zlim(-0.5, zlim)
ax.set_xlim(0, 3.0)
ax.set_ylim(0, 4.0)
#plt.savefig('pos_a_3dF7.png')
plt.show()
#for k_single check
'''plt.plot(al1,kn(al1,al1))
plt.xlabel('a1r1')
plt.ylabel('a2r2a1r1,a1r1)')
plt.savefa2r2('Ktest.png') 
plt.show()'''
