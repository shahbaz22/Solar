'working program for K surf plot for a2, a1 values, older code using np.where'
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
R2=1.0
'a1r1 already used in kn function'
# al1=np.linspace(0.1,3,800) # Bessel functions only take positive values
# al2=np.linspace(-4.01,4,800) #Assures that the arrays have the same number of elements
#al2test=np.linspace(0.01,-3,800)
#al1test=np.linspace(0.1,1,800)
# alpha1,alpha2=np.meshgrid(al1,al2)
#alpha1t, alpha2t=np.meshgrid(al1test,al2test)
scsp=scipy.special
'Function to take a1r1 and a1r1 vales can return 2D K'
'Used for making surface pla2r2 of a1r1,a1r1,K'
'check if def kn can take ara2r2y'
def kn(alpha1,alpha2):
    'Variables for alpha and R values'
    'Simplified meshgid vals'     
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
    def c2p(a1,a2):
        'c2 constant for positive a2'
        c21=j(0,a1)*j(1,a2)-j(1,a1)*j(0,a2)
        return c21
    def b2p(a1,a2):
        'b2 constant for positive a2'
        b21=j(1,a1)*y(0,a2)-j(0,a1)*y(1,a2)
        return b21
    def c2n(a1,a2):
        'c2 constant for negative a2'
        c22=j(0,a1)*j(1,a2)+j(1,a1)*j(0,a2)
        return c22
    def b2n(a1,a2):
        'b2 constant for negative a2'
        b22=-j(1,a1)*y(0,a2)-j(0,a1)*y(1,a2)
        return b22
    'producing constants for different a2 values, used for B1'
    c2prodp=np.where((alpha1>0) & (alpha2>0),c2p(a1r1,a2r1),0)
    b2prodp=np.where((alpha1>0) & (alpha2>0),b2p(a1r1,a2r1),0)
    c2prodn=np.where((alpha1>0) & (alpha2<0),c2n(a1r1,a2r1),0)
    b2prodn=np.where((alpha1>0) & (alpha2<0),b2n(a1r1,a2r1),0)

    #print(b2prodp)
    ctb1rat=j(1,a1r1)*(1/alpha1-1/alpha2)
    'if there is an error double check sign of dis and ctb1rat'
    dis=2/(np.pi*a2r1)

    'ratio of constants for F, depend on sign of a2'
    c2b2ratp=np.where((alpha1>0) & (alpha2>0),c2p(a1r1,a2r1)/b2p(a1r1,a2r1),0)
    c2b2ratn=np.where((alpha1>0) & (alpha2<0),c2n(a1r1,a2r1)/b2n(a1r1,a2r1),0)

    def F0(x,c2b2rat):
    	f=j(0,x)+c2b2rat*y(0,x)
    	return f
    def F1(x,c2b2rat):
    	f=j(1,x)+c2b2rat*y(1,x)
    	return f
    'norm for positive a2'
    def norm(b2prod,c2b2rat):
        norm0=(2*np.pi*(R2*b2prod*(F1(a2r2,c2b2rat)/(abs(alpha2)*dis))+R1*ctb1rat))   
        return norm0
    'B1 is the same if a2>0 or a2<0'
    'functon within function used to avoid zero terms in b2prod as 1/0 gives error'
    B1p=np.where((alpha1>0) & (alpha2>0),1/norm(b2p(a1r1,a2r1),c2b2ratp),0)
    B2p=B1p*b2prodp/dis
    C2p=B1p*c2prodp/dis
    
    B1n=np.where((alpha1>0) & (alpha2<0),1/norm(b2n(a1r1,a2r1),c2b2ratn),0)
    B2n=B1n*b2prodn/dis
    C2n=B1n*c2prodn/dis   

    def k(B1,B2,c2b2rat,r1,r2):
        'function calculates K for pos and neg a2 depending on constants'
        'for pos a2 B1p,B2p and c2b2ratp'
        'for neg a2 B1n B2n and c2b2ratn'
        'k1 only uses a1 and r1'
        j0=j(0,alpha1*r1)
        j1=j(1,alpha1*r1)
        k1=(2*np.pi*B1**2)*((r1**2)*(j0**2+j1**2)-2*r1*j0*j1/alpha1)/alpha1
        print(k1,'k1')        
        'k2 only uses a2'
        f1r1=F1(abs(alpha2*r1),c2b2rat)
        f1r2=F1(abs(alpha2*r2),c2b2rat)
        f0r1=F0(abs(alpha2*r1),c2b2rat)
        f0r2=F0(abs(alpha2*r2),c2b2rat)
        sgn=alpha2/abs(alpha2)

        k2=((r2**2)*(f0r2**2+f1r2**2)-2*r2*f0r2*f1r2/abs(alpha2))-(r1**2)*(f0r1**2+f1r1**2)
        k2=(k2+2*r1*f0r1*f1r1/abs(alpha2))*(2*np.pi*B2**2)/abs(alpha2) #here whole term multipied by constant
        k2=k2+(4*np.pi*B2*B1*r1*ctb1rat*(f0r1-f0r2))/abs(alpha2)
        k2=k2*sgn
        print(k2,'k2')
        return k1+k2
    kp1=np.where((alpha1>0) &(alpha2>0),k(B1p,B2p,c2b2ratp,R1,R2),0)
    kn1=np.where((alpha1>0) &(alpha2<0),k(B1n,B2n,c2b2ratn,R1,R2),0)
    helicity=kp1+kn1
    return helicity
print(kn(2,-4))
'3d plot code'
# fig=plt.figure()
# ax=fig.gca(projection='3d')
# fn = kn(alpha1,alpha2)
# zmax = 2.5
# zmin= -0.5 #can log norm function take negative values?
# #norm = cc.LogNorm(vmin=zmin, vmax=zmax)
# norm=cc.PowerNorm(gamma=1)
# surf=ax.plot_surface(alpha1,alpha2, fn, cmap=cm.hot,
#                      linewidth=0,antialiased=False,
#                      norm=norm)
# ax.set_xlabel('a1')
# ax.set_ylabel('a2')
# ax.set_zlabel('K')
# #Add color bar
# fig.colorbar(surf,shrink=0.5,aspect=5) 
# ax.set_zlim(zmin, zmax)
# #ax.set_xlim(0, 3.0)
# #ax.set_ylim(0, 4.0)
# print(max(al1),min(al1))
# #plt.savefig('pos_neg_a2_3d.png')
# plt.show()


#for k_single check
# plt.plot(al1,kn(al1,al1))
# plt.xlabel('a')
# plt.ylabel('K')
# plt.savefig('ksingle.png') 
# plt.show()
