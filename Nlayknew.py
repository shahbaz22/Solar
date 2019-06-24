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
a1,a2,a3=3,1,-2
'Array of helicity "a" used as inputs'

a=[a1,a2,a3]
r=[0.25,0.5,1]
print('here is a',a)
def get_helicityn(a):
    'function to get n values of helicty for n layers'
    'r is array for radii such that a[1]->r[1]...a[n]->r[n]'
    '''r=[]
    'for loop to generate array of r values coresp. to n alpha values'
    for i in np.linspace(len(a),1,len(a)):
    r=1/i
    r.append(r)'''
    print('here is r',r)
    #print(r[0],r[2])
    #print(len(r))
    #print(a[0])
    #print(a[2])


    'if cbrat=0 then j obtained'
    'cbr determins order of function, eg. G0,G1 when c3/b3 '
    def F0(x,cbrat):
        f=j(0,abs(x))+cbrat*y(0,abs(x))
        return f
    def F1(x,cbrat):
        f=j(1,abs(x))+cbrat*y(1,abs(x))
        return f
    def y(v,x):
        BF22=scsp.yv(v,abs(x))
        return (BF22)
    def j(v,x):
        BF11=scsp.jv(v,abs(x))
        return (BF11)
    def sig(n1,n2):
        siga=np.sign(a[n1])*np.sign(a[n2])
        return siga
    def sig1(n1):
        siga=np.sign(a[n1])
        return siga

    'Fist step is to define constats B2...BN,C2....CN and functions F'
    'Next function value depends on previous value of constat'
    'Next value of constat depends of previous function'
    'indexing of array starts from zero s.t a1=a[0]'
    'j1(a1r1)->F1(a2r2),j1(a2r1)->j1(a3r2)'
    'Patter; next value of args. +1 '
    'if both args. same order (num.) fuction order increases by 1'
    'otherwise function remains the same'
    'if same args. always +1 then always have diff order'
    'hence some functions always remain same'
    'increasing order F0->G0... ,F1->G1 wrt 3 lay code'
    c=np.zeros(len(a))
    'creates an array of zeros the size of a to store c values'
    b=np.zeros(len(a))
    'creates an array of zeros the size of a to store c values'
    cbr=np.zeros(len(a))
    'creates an array of zeros the size of a to store c/b values (cbrat in old code)'
    c2prod=j(0,a[0]*r[0])*j(1,a[1]*r[0])-sig(0,1)*j(1,a[0]*r[0])*j(0,a[1]*r[0])
    b2prod=sig(0,1)*j(1,a[0]*r[0])*y(0,a[1]*r[0])-j(0,a[0]*r[0])*y(1,a[1]*r[0])
    'first element c[0] and b[0] calculated outside of loop below to minimize # of if statments within loop'
    c[0]=c2prod
    b[0]=b2prod
    #cbr[0]=np.where(b[0]==0,0,c[0]/b[0])
    'cbr is the ratio of c and b terms, will have len(a)-2 terms'
    'rememebr python index is one below actual index as it starts from 0'
    'start from 1 in the for loop and end at index len(a)-1'
    'here i is the index value'
    print('np.arange(1,len(a),1)',np.arange(1,len(a),1))
    print('len(a)-2',len(a)-2)
    for i in np.arange(1,len(a),1): 
      if i<=(len(a)-2):
        cbr[i]=np.where(b[i-1]==0,0,c[i-1]/b[i-1])
        print('cbr in loop',cbr[i])
        'in cprod the previous value of cbr is required' 
        'eg. for B3 (b[1]), F1 is used and reqs. c2/b2==c[0]/b[0]==cbr[1]'
        cprod=F0(a[i]*r[i],cbr[i])*j(1,a[i+1]*r[i])
        c[i]=cprod-sig(i,i+1)*F1(a[i]*r[i],cbr[i])*j(0,a[i+1]*r[i])
        bprod=sig(i,i+1)*F1(a[i]*r[i],cbr[i])*y(0,a[i+1]*r[i])
        b[i]=bprod-F0(a[i]*r[i],cbr[i])*y(1,a[i+1]*r[i])
        'i+1 must not exceed len(a)'
    print('here is c',c)
    print('here is b',b)
    'next step after determining constats is as follows'
    'generalise dis terms eg, dis[i]=pi*r[i]*a[i+1]/2'
    'bcn[i]=dis[i-1]*bn[i]*bcn[i-1], B2 needs to be completed outside for loop'
    'next find patter in psi and recalc bcn eg B1,B2,B3 normalized'
    'Psi0 will keep increasing with number of layers chaging B1'
    'rememebr index zero is first element in python'
    'B2=b[0], B3=b[1], ect.'
    'for loop below used to recalculate cbr values as final cbr value outside of'
    'range of above for loop, due to final value of c,b depending on cbr[i-1]'
    for i in np.arange(0,len(a),1):
    	if i<=(len(a)):
    		cbr[i]=np.where(b[i]==0,0,c[i]/b[i])
    print('here is cbr',cbr)

    dis=np.zeros(len(a))
    bn=np.zeros(len(a))
    'psi array will be useful for checks'
    psi=np.zeros(len(a))
    'dis[0]==dis2'
    dis[0]=pi*r[0]*abs(a[1])/2
    bn[0]=dis[0]*b[0]
    psi0=(2*pi*b[0]*dis[0]/abs(a[1]))*(r[1]*F1(a[1]*r[1],cbr[0])-r[0]*F1(a[1]*r[0],cbr[0]))
    psi[0]=psi0+(2*pi*r[0]/abs(a[0]))*j(1,a[0]*r[0])
    
    #psi0=r[0]*b[0]*(F1(a[0]*r[0],cbr[0])/(abs(a[1])))
    #psi[0]=2*pi*(psi0+j(1,a[0]*r[0])*(1/abs(a[0])-1/abs(a[1])))
    print('psi[0]',psi[0])
    'bn term computed outside of loop to allow for'
    'bn[i-1] is the recursion term, cbr[0]==c2/b2'
    'for loop starts at 1'
    for i in np.arange(1,len(a),1):
      if i<=(len(a)-2):
        'for dis[1]==dis3, use a[2]==a3'
        dis[i]=pi*r[i]*abs(a[i+1])/2
        bn[i]=dis[i]*b[i]*bn[i-1]
        psi[i]=(2*pi*bn[i]/abs(a[i+1]))*(r[i+1]*(r[1+i]*F1(a[i+1]*r[i+1],cbr[i])-r[i]*(F1(a[i+1]*r[i],cbr[i])))) 
    psit=np.sum(psi)
    b1=1/psit
    print('here is b1',b1)
    'element-wise multiplication using np.dot'
    bn=np.dot(b1,bn)
    print('here is bn',bn)
    'next step to compare values for 3layer for bn and b1~psi'
    'next section is the helicity first values of helicity ca be calculated outside for loop'
    k=np.zeros(len(a))
    'k1'
    k[0]=(r[0]**2)*(j(0,a[0]*r[0])**2+j(1,a[0]*r[0])**2)
    k[0]=((2*pi*b1**2)/abs(a[0]))*(k[0]-2*(r[0]/abs(a[0]))*j(0,a[0]*r[0])*j(1,a[0]*r[0]))
    'if there is a print error it is due to a syntax error in the called equation'
    print('k[0]',k[0])
    'reason for k2==k[1] perfomed outside the for loop is because the k2 term is not general'
    'k2 term requires the first value of ctbrat and in general the new value and the previous value ctbrat required'
    'cannot start loop from a lower value than 1'
    ctbr=np.zeros(len(a))
    ctbr[0]=0

    #ctb1rat=B1*R1*j(1,a1r1)*(1/alpha1-sig2/abs(alpha2))    
    #k2=((R2**2)*(F0(a2r2)**2+F1(a2r2)**2)-2*R2*F0(a2r2)*F1(a2r2)/abs(alpha2))-(R1**2)*(F0(a2r1)**2+F1(a2r1)**2)
    #k2=(k2+2*R1*F0(a2r1)*F1(a2r1)/abs(alpha2))*(2*pi*B2**2)/abs(alpha2) #here whole term multipied by constant
    #k2=k2+(4*pi*B2*ctb1rat*(F0(a2r1)-F0(a2r2)))/abs(alpha2)
    #k2=k2*sig1(1)

    ctbr[1]=b1*r[0]*j(1,a[0]*r[0])*(1/abs(a[0])-sig(0,1)/abs(a[1]))    
    k2=((r[1]**2)*(F0(a[1]*r[1],cbr[0])**2+(F1(a[1]*r[1],cbr[0])**2)-2*r[1]*F0(a[1]*r[1],cbr[0])*F1(a[1]*r[1],cbr[0])/abs(a[1]))-(r[0]**2)*(F0(a[1]*r[0],cbr[0])**2+F1(a[1]*r[0],cbr[0])**2))
    k2=(k2+2*r[0]*F0(a[1]*r[0],cbr[0])*F1(a[1]*r[0],cbr[0])/abs(a[1]))*(2*pi*b[0]**2)/abs(a[1]) #here whole term multipied by constant
    k2=k2+(4*pi*b[0]*ctbr[1]*(F0(a[1]*r[0],cbr[0])-F0(a[1]*r[1],cbr[0])))/abs(a[1])
    k2=k2*sig1(1)  
    print('k2',k2)
#    for i in np.arange(1,len(a),1):
#      if i<=(len(a)-2):
  


print(get_helicityn(a))


