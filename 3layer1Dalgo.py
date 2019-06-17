'working program using for loops on a 3 layer model over arange of r'
'with B4=0 and including a3 which can be negative'
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
R2=0.75
R3=1
'a1r1 already used in kn function'
a1=np.linspace(0.1,3,800) # Bessel functions only take positive values
a2=np.linspace(-4.01,4,800) #Assures that the arrays have the same number of elements
a3=np.linspace(0.1,3,800)
#a1,a2=np.meshgrid(a1,a2)

scsp=scipy.special
pi=np.pi
'Function to take a1r1 and a1r1 vales can return 2D K'
'Used for making surface pla2r2 of a1r1,a1r1,K'
'check if def kn can take ara2r2y'
def helicity(alpha1,alpha2,alpha3):
    'Variables for alpha and R values'
    'Simplified meshgid vals' 
    'The order of computation being important'    
    def j(v,a):
    	'Function to return bessel of 1st kind wrt order v,and R,a'
    	BF11=scsp.jv(v,a)
    	return (BF11)
    def y(v,a):
    	'Function to return bessel of 2st kind wrt order v,and R,a'
    	BF22=scsp.yv(v,a)
    	return (BF22)
    'Simplified Variables'
    a1r1=alpha1*R1
    a2r1=abs(alpha2*R1)
    a2r2=abs(alpha2*R2)
    a3r3=abs(alpha3*R3)
    a3r2=abs(alpha3*R2)

    'sig reverses sign for constant c2 and b2 depending on sign if alpha2'
    'sig1 always positve so disregard, conseq. sig12==sig2'
    sig2=np.sign(alpha2) 
    sig3=np.sign(alpha3)
    sig23=np.sign(alpha2)*np.sign(alpha3)

    def F0(x):
        f=j(0,x)+cb2rat*y(0,x)
        return f
    def F1(x):
        f=j(1,x)+cb2rat*y(1,x)
        return f
    def G0(x):
        f=j(0,x)+cb3rat*y(0,x)
        return f
    def G1(x):
        f=j(1,x)+cb3rat*y(1,x)
        return f

    'constants c2prod and b2prod computed first as no dep. on F'
    c2prod=j(0,a1r1)*j(1,a2r1)-sig2*j(1,a1r1)*j(0,a2r1)
    b2prod=sig2*j(1,a1r1)*y(0,a2r1)-j(0,a1r1)*y(1,a2r1)
    'F requires c2b2rat'
    print(c2prod,'here is c2prod')
    print(b2prod,'here is b2prod')
    cb2rat=np.where(b2prod==0,0,c2prod/b2prod)
    print(cb2rat,'cb2rat')
    'c3prod and b3prod lower down as they depend on F function'
    'assume y function always included in higher orders'
    c3prod=F0(a2r2)*j(1,a3r2)-sig23*F1(a2r2)*j(0,a3r2)
    b3prod=sig23*F1(a2r2)*y(0,a3r2)-F0(a2r2)*y(1,a3r2)
    print('b3prod',b3prod)
    print('c3prod',c3prod)
    cb3rat=np.where(b3prod==0,0,c3prod/b3prod)
    print('cb3rat',cb3rat)
    'compared to previous code 1/dis used to simp. calculations'
    dis1=(pi*a2r1)/2
    dis2=(pi*a3r2)/2
    'want to find B1 and set ps0 eventually to to 1'
    'Calculating normalisation constants'
    'Only factos of B1 in psi0'
    'if a1,a2,a3=a1,a1,a1 b3prod=0 s.t B3=0'
    'When B3=0 ps1 (3rd lay cont.) is zero'
    B2=dis1*b2prod
    B3=dis2*B2*b3prod
    ctbdiff=j(1,a1r1)*(1/abs(alpha1)-sig2/abs(alpha2))
    'ctbdiff known as ctb1rat in 2lay code'
    psi0=2*pi*(R2*B2*(F1(a2r2)/(abs(alpha2)))+R1*ctbdiff)
    psi1=(2*pi*B3/abs(alpha3))*(R3*G1(a3r3)-R2*G1(a3r2))      
    print(psi0,'psi0')
    B1=1/(psi0+psi1) 
    B2=dis1*b2prod*B1
    B3=dis2*B2*b3prod
    print(B1,'B1')
    'now calculaing terms for k'
    k1=(2*pi*B1**2)*((R1**2)*(j(0,a1r1)**2+j(1,a1r1)**2)-2*R1*j(0,a1r1)*j(1,a1r1)/alpha1)/alpha1
    print(k1,'k1')
    'k2 only uses a2'
    ctb1rat=B1*R1*j(1,a1r1)*(1/alpha1-sig2/abs(alpha2))    
    k2=((R2**2)*(F0(a2r2)**2+F1(a2r2)**2)-2*R2*F0(a2r2)*F1(a2r2)/abs(alpha2))-(R1**2)*(F0(a2r1)**2+F1(a2r1)**2)
    k2=(k2+2*R1*F0(a2r1)*F1(a2r1)/abs(alpha2))*(2*pi*B2**2)/abs(alpha2) #here whole term multipied by constant
    k2=k2+(4*pi*B2*ctb1rat*(F0(a2r1)-F0(a2r2)))/abs(alpha2)
    k2=k2*sig2
    print(k2,'k2')
    'for k3: a2-->a3,R1-->R2,R2-->R3,sig2-->sig3,B1,B2-->B2,B3 ect. compared to k2'    
    ctb2rat=B2*R2*F1(a2r2)*(1/abs(alpha2)-sig23/abs(alpha3))+ctb1rat 
    k3=((R3**2)*(G0(a3r3)**2+G1(a3r3)**2)-2*R3*G0(a3r3)*G1(a3r3)/abs(alpha3))-(R2**2)*(G0(a3r2)**2+G1(a3r2)**2)
    k3=(k3+2*R2*G0(a3r2)*G1(a3r2)/abs(alpha3))*(2*pi*B3**2)/abs(alpha3) #here whole term multipied by constant
    k3=k3+(4*pi*B3*ctb2rat*(G0(a3r2)-G0(a3r3)))/abs(alpha3)
    k3=k3*sig3
    print(k3,'k3')

    helicity=k1+k2+k3
    return helicity
print(helicity(3,-1,-1))
'''fig=plt.figure()
ax=fig.gca(projection='3d')
fn = helicity(a1,a2,a2)
zmax = 2.5
zmin= -0.5 #can log norm function take negative values?
#norm = cc.LogNorm(vmin=zmin, vmax=zmax)
norm=cc.PowerNorm(gamma=1)
surf=ax.plot_surface(a1,a2, fn, cmap=cm.hot,
                     linewidth=0,antialiased=False,
                     norm=norm)
ax.set_xlabel('a1')
ax.set_ylabel('a2')
ax.set_zlabel('K')
#Add color bar
fig.colorbar(surf,shrink=0.5,aspect=5) 
#ax.set_zlim(zmin, zmax)
#ax.set_xlim(0, 3.0)
#ax.set_ylim(0, 4.0)
#print(max(al1),min(al1))
#plt.savefig('pos_neg_a2_3d.png')
plt.show()
'''
#for k_single check
'''plt.plot(a1,helicity(a1,a1,a1))
plt.xlabel('a1r1')
plt.ylabel('K')
#plt.savefa2r2('Ktest.png') 
plt.show()'''