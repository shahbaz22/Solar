'program to find smallest root of a high order polynomial by using a number of initial starting points'
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import special
from scipy.optimize import fsolve
a2=np.linspace(-4.01,4,800) #Assures that the arrays have the same number of elements
scsp=scipy.special
xv=np.linspace(-5,5,1000)
def quartf(a,b,c,d,x):
	'constants infront of coeficent terms'
	y=a*x**3+b*x**2+c*x+d
	return y
# plt.plot(xv,quartf(1,2,-1,-1,xv))
# plt.axhline(color='r',linestyle=':')
# plt.ylim([-10,10])
# plt.ylabel('y')
# plt.xlabel('x')
# plt.show()
def j(v,a):
    'Function to return bessel of 1st kind wrt order v,and R,a'
    BF11=scsp.jv(v,a)
    return (BF11)
def B(a):
    B0=(a/(2*np.pi*j(1,a)))
    return B0
print(B(0.1),'B')