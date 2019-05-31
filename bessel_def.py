'non working program for K surf plot for'
'a2 values >0,=0,<0, a1 always >0'
import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.cm as cm
import matplotlib.colors as cc
from scipy import special
a2=np.linspace(1,500,100000)*0.5
scsp=scipy.special

def j(v,a):
    'Function to return bessel of 1st kind wrt order v,and R,a'
    BF11=scsp.jv(v,a)
    return (BF11)
def y(v,a):
    'Function to return bessel of 2st kind wrt order v,and R,a'
    BF22=scsp.yv(v,a)
    return (BF22)
delta=y(0,a2)*j(1,a2)-y(1,a2)*j(0,a2)
plt.plot(a2,delta)
plt.show()