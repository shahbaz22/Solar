'non working program using for N layer model over arange of r'
'use num'
import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.cm as cm
import matplotlib.colors as cc
from scipy import special
scsp=scipy.special
pi=np.pi

#a1,a2,a3,a4,a5=1,2,3,4,5
a1,a2,a3=3,3,2
'Array of helicity "an" used as inputs'

an=[a1,a2,a3]
print('here is an',an)
def get_helicityn(a):
    'function to get n values of helicty for n layers'
    'rn is array for radii such that a[1]->rn[1]...a[n]->rn[n]'
    '''rn=[]
    'for loop to generate array of r values coresp. to n alpha values'
    for i in np.linspace(len(an),1,len(an)):
    r=1/i
    rn.append(r)'''
    rn=[0.25,0.5,1]
    print('here is rn',rn)
    #print(rn[0],rn[2])
    #print(len(rn))
    #print(an[0])
    #print(an[2])


    'if cbrat=0 then j obtained'
    'cbrat determins order of function, eg. G0,G1 when c3/b3 '
    def F0(x,cbrat):
        f=j(0,abs(x))+cbrat*y(0,abs(x))
        return f
    def F1(x,cbrat):
        f=j(1,abs(x))+cbrat*y(1,abs(x))
        return f
    def y(v,a):
        BF22=scsp.yv(v,abs(a))
        return (BF22)
    def j(v,a):
        BF11=scsp.jv(v,abs(a))
        return (BF11)
    def sig(n1,n2):
        siga=np.sign(an[n1])*np.sign(an[n2])
        return siga

    'Fist step is to define constants B2...BN,C2....CN and functions F'
    'Next function value depends on previous value of constant'
    'Next value of constant depends of previous function'
    'indexing of array starts from zero s.t a1=an[0]'
    'j1(a1r1)->F1(a2r2),j1(a2r1)->j1(a3r2)'
    'Pattern; next value of args. +1 '
    'if both args. same order (num.) fuction order increases by 1'
    'otherwise function remains the same'
    'if same args. always +1 then always have diff order'
    'hence some functions always remain same'
    'increasing order F0->G0... ,F1->G1 wrt 3 lay code'
    c=np.zeros(len(an))
    'creates an array of zeros the size of an to store c values'
    b=np.zeros(len(an))
    'creates an array of zeros the size of an to store c values'
    cbrat=np.zeros(len(an))
    'creates an array of zeros the size of an to store c/b values'
    c2prod=j(0,an[0]*rn[0])*j(1,an[1]*rn[0])-sig(0,1)*j(1,an[0]*rn[0])*j(0,an[1]*rn[0])
    b2prod=sig(0,1)*j(1,an[0]*rn[0])*y(0,an[1]*rn[0])-j(0,an[0]*rn[0])*y(1,an[1]*rn[0])
    'first element c[0] and b[0] calculated outside of loop below to minimize # of if statments within loop'
    c[0]=c2prod
    b[0]=b2prod
    'cbrat is the ratio of c and b terms, will have len(an)-2 terms'
    'rememebr python index is one below actual index as it starts from 0'
    'start from 1 in the for loop and end at index len(an)-1'
    'here i is the index value'
    for i in np.arange(1,len(an),1): 
      if i<=(len(an)-2):
        cbrat[i]=np.where(b[i-1]==0,0,c[i-1]/b[i-1])
        'in cprod the previous value of cbrat is required' 
        'eg. for B3 (b[1]), F1 is used and reqs. c2/b2==c[0]/b[0]==cbrat[1]'
        cprod=F0(an[i]*rn[i],cbrat[i])*j(1,an[i+1]*rn[i])
        c[i]=cprod-sig(i,i+1)*F1(an[i]*rn[i],cbrat[i])*j(0,an[i+1]*rn[i])
        bprod=sig(i,i+1)*F1(an[i]*rn[i],cbrat[i])*y(0,an[i+1]*rn[i])
        b[i]=bprod-F0(an[i]*rn[i],cbrat[i])*y(1,an[i+1]*rn[i])
        'i+1 must not exceed len(an)'
    print('here is c',c)
    print('here is b',b)
    'next step after determining constants is as follows'
    'generalise dis terms eg, dis[i]=pi*rn[i]*an[i+1]/2'
    'bcn[i]=dis[i-1]*bn[i]*bcn[i-1], B2 needs to be completed outside for loop'
    'next find pattern in psi and recalc bcn eg B1,B2,B3 normalized'
    'Psi0 will keep increasing with number of layers changing B1'
    'rememebr index zero is first element in python'
    'B2=b[0], B3=b[1], ect.'

    dis=np.zeros(len(an))
    bn=np.zeros(len(an))
    'psi array will be useful for checks'
    psi=np.zeros(len(an))
    'dis[0]==dis2'
    dis[0]=pi*rn[0]*an[1]/2
    bn[0]=dis[0]*b[0]
    psi0=rn[0]*b[0]*(F1(an[0]*rn[0],cbrat[0])/(abs(an[1])))
    psi[0]=2*pi*(psi0+j(1,an[0]*rn[0])*(1/abs(an[0])-1/abs(an[1])))

    'bn term computed outside of loop to allow for'
    'bn[i-1] recursion term'
    'for loop starts at 1'
    for i in np.arange(1,len(an),1):
      if i<=(len(an)-2):
        'for dis[1]==dis3, use an[2]==a3'
        dis[i]=pi*rn[i]*an[i+1]/2
        bn[i]=dis[i]*b[i]*bn[i-1]
        psi[i]=(2*pi*bn[i]/abs(an[i+1]))*(rn[i+1]*(F1(an[i+1]*rn[i+1],cbrat[i+1])-rn[i]*(F1(an[i+1]*rn[i],cbrat[i])))) 
    
    psit=np.sum(psi)
    b1=1/psit
    print('here is b1',b1)
    'element-wise multiplication using np.dot'
    bn=np.dot(b1,bn)
    print('here is bn',bn)
    'next step to compare values for 3layer for bn and b1~psi'

print(get_helicityn(an))


