import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import special
#in this program I want to seperate the values of a1 and a2 depending on sign and perform some function
# This is a check for a git upload
R1=0.5
R2=1
a1=np.linspace(0,4,600) # Bessel functions only take positive values
a2=np.linspace(-4,4,600) #Assures that the arrays have the same number of elements

a1,a2=np.meshgrid(a1,a2)

#np.where preserves shape of array
mna2=np.zeros(np.shape(a2))
mask_nega2=np.where(a2<0,a2,mna2)
mask_posa2=np.where(a2>0,a2,mna2)
red_mask_posa2=mask_posa2[300:600]
red_mask_nega2=mask_nega2[0:300]
#k=np.concatenate((red_mask_posa2,red_mask_nega2), axis=0)
#print(red_mask_posa2,red_mask_posa2.ndim)
#print('reda2vals',red_mask_nega2,np.shape(red_mask_nega2))
print(np.shape(k))
print(k)
print(a2)