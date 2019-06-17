'working program using for loops on a 2 layer model'
'Fixing this code was more difficult than expected'
'I missed an a2r2 value in c2prod, PRINTING IS YOUR FRIEND'
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import special
R1=0.5
R2=1.0
'a1r1 already used in kn function'
a1=np.linspace(0.1,3,800) # Bessel functions only take positive values
a2=np.linspace(-4.01,4,800) #Assures that the arrays have the same number of elements
scsp=scipy.special
'Function to take a1r1 and a1r1 vales can return 2D K'
'Used for making surface pla2r2 of a1r1,a1r1,K'
'check if def kn can take ara2r2y'
def helicity(alpha1,alpha2):
    'Variables for alpha and R values'
    'Simplified meshgid vals'     
    def j(v,a):
    	'Function to return bessel of 1st kind wrt order v,and R,a'
    	BF11=scsp.jv(v,a)
    	return (BF11)
    def y(v,a):
    	'Function to return bessel of 2st kind wrt order v,and R,a'
    	BF22=scsp.yv(v,a)
    	return (BF22)
    a1r1=alpha1*R1
    a2r1=abs(alpha2)*R1
    a2r2=abs(alpha2*R2)
    #sig=alpha2/abs(alpha2)  
    sig=np.sign(alpha2)  
    'sig reverses sign for constant c2 and b2 depending on sign if alpha2'


    c2prod=j(0,a1r1)*j(1,a2r1)-sig*j(1,a1r1)*j(0,a2r1)
    b2prod=sig*j(1,a1r1)*y(0,a2r1)-j(0,a1r1)*y(1,a2r1)
    print(c2prod,'c2prod')
    print(b2prod,'b2prod')
    #print(sig)

    #print(np.shape(b2prod))


    ctb1rat=j(1,a1r1)*(1/alpha1-1/alpha2)
    'if there is an error double check sign of dis and ctb1rat'
    dis=2/(np.pi*a2r1)

    c2b2rat=c2prod/b2prod
    print(c2b2rat,'c2b2rat')

    def F0(x):
    	f=j(0,x)+c2b2rat*y(0,x)
    	return f
    def F1(x):
    	f=j(1,x)+c2b2rat*y(1,x)
    	return f
    'norm for positive a2 and negative a2'
    B1=1/(2*np.pi*(R2*b2prod*(F1(a2r2)/(abs(alpha2)*dis))+R1*ctb1rat))   
    print(1/B1,'psi0')
    B2=B1*b2prod/dis
    C2=B2*c2b2rat
    print(B1)       
    'B1 is the same if a2>0 or a2<0'
    'functon within function used to avoid zero terms in b2prod as 1/0 gives error'

    def k(alpha1,alpha2,r1,r2):
        'function calculates K for pos and neg a2 depending on constants'
        'k1 only uses a1 and r1'
        j0=j(0,alpha1*r1)
        j1=j(1,alpha1*r1)
        k1=(2*np.pi*B1**2)*((r1**2)*(j0**2+j1**2)-2*r1*j0*j1/alpha1)/alpha1

        'k2 only uses a2'
        f1r1=F1(a2r1)
        f1r2=F1(a2r2)
        f0r1=F0(a2r1)
        f0r2=F0(a2r2)
        sgn=alpha2/abs(alpha2)

        k2=((r2**2)*(f0r2**2+f1r2**2)-2*r2*f0r2*f1r2/abs(alpha2))-(r1**2)*(f0r1**2+f1r1**2)
        k2=(k2+2*r1*f0r1*f1r1/abs(alpha2))*(2*np.pi*B2**2)/abs(alpha2) #here whole term multipied by constant
        k2=k2+(4*np.pi*B2*B1*r1*ctb1rat*(f0r1-f0r2))/abs(alpha2)
        k2=k2*sgn
        print(k1,'k1')
        print(k2,'k2')

        return k1+k2
    
    helicity=k(alpha1,alpha2,R1,R2)
    return helicity
print(helicity(3,-1))
#'for k_single check'
'''plt.plot(a1,helicity(a1,a1))
plt.xlabel('a1r1')
plt.ylabel('K')
#plt.savefa2r2('Ktest.png') 
plt.show()'''