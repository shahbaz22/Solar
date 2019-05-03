'not working porgram for a1,a2>0 for W with surface plot'
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
al1=np.linspace(0.5,3,200) # Bessel functions only take positive values
al2=np.linspace(0.5,3,200) #Assures that the arrays have the same number of elements
alpha1,alpha2=np.meshgrid(al1,al2)
scsp=scipy.special
'Function to take a1 and a2 vales can return 2D K'
'Used for making surface plot of a1,a2,K'
'check if def kn can take array'
def Wn(alpha1,alpha2):
    
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
    b2prod=-j(1,a1)*y(0,a2)-j(0,a1)*y(1,a2)
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
    B1=1/(2*np.pi*(R2*b2prod*(F1(ar)/(abs(alpha2)*dis))+R1*ctb1rat))   
    B2=B1*b2prod/dis
    C2=B2*c2b2rat

    CTHETA = B1*R1*ctb1rat
    #mu0=4*np.pi*10**-7
    'calculates first term K1  in helicity (helicity of inner layer) function k1'
    'Fist term only uses R1,(Rc) and L set to 1'
    def w1(alpha1,b1,r1):
        
        j0=j(0,alpha1*r1)
        j1=j(1,alpha1*r1)
        w1=(b1**2)*((r1**2)*(j0**2+j1**2)-r1*j0*j1/alpha1)
        return w1
#Calculates term k2 in helicity
#Works for postive and negative alpha2 - "sgn" is sign of alpha2 function k2 
    def w2(alpha2,b2,r1,r2):
        'all functions of a2'
        f1c=F1(abs(alpha2*r1))
        f1=F1(abs(alpha2*r2))
        f0c=F0(abs(alpha2*r1))
        f0=F0(abs(alpha2*r2))
        sgn=alpha2/abs(alpha2)

        w2=(r2**2)*(f0**2+f1**2)-r2*f0*f1/abs(alpha2)-(r1**2)*(f0c**2+f1c**2)
        w2=(w2+r1*f0c*f1c/abs(alpha2))*(b2**2)
        w2=sgn*w2 # apart from K1 the other part of Ks sign depends on the sign of a2
        return w2
    w=(w1(alpha1,B1,R1)+w2(alpha2,B2,R1,R2)) # multi. my constant
    return w

fig=plt.figure()
ax=fig.gca(projection='3d')
fn = Wn(alpha1,alpha2)
print(np.max(fn))
#fn_mod = np.ma.masked_where(fn>2.5, fn)
zlim = np.max(fn)
#fn_mod = np.where(fn>zlim, zlim, fn)
norm = cc.LogNorm(vmin=0.1, vmax=zlim)
surf=ax.plot_surface(alpha1,alpha2, fn, cmap=cm.hot,
   					 linewidth=0,antialiased=False,
   					 norm=norm)
ax.set_xlabel('a1')
ax.set_ylabel('a2')
ax.set_zlabel('W')
#Add color bar
fig.colorbar(surf,shrink=0.5,aspect=5) 
#ax.set_zlim(0, zlim)
#plt.ylim(2,1)
#plt.xlim(1,2)
plt.savefig('pos_a_3d_W.png')
plt.show()
"""
plt.plot(al2,kn(al2,al2))
plt.show()
"""