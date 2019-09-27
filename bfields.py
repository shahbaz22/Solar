'working program using  mutli-layer model to calculate force-free fields for Bz, Bt, alpha vs r'
import numpy as np
import scipy
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.cm as cm
import matplotlib.colors as cc
from scipy import special
scsp=scipy.special
pi=np.pi

'Array of helicity "a" used as inputs'
a1 =[ 1,2,3,-1,4]
r1=[0.25,0.5,0.6,0.75,1]
# r1=np.linspace(0.01,1,len(a1))
# r1=r1[1:]

# print('here is a1',a1)
#print('here is r1',r1)
#r=[0.5,0.75,1]
# r=[]
# 'for loop to generate array of r values coresp. to n alpha values'
# for i in np.linspace(len(a),1,len(a)):
#   ra=1/i
#   r.append(ra)

# print('here is a',a)
# print('here is r',r)
def getb(a,r):
    'function to get n values of helicty for n layers'
    'r is array for radii such that a[1]->r[1]...a[n]->r[n]'
    'if cbrat=0 then j obtained'
    'cbr determins order of function, eg. G0,G1 when c3/b3 '
    
    if len(a)!=len(r):
        raise ValueError('dimension of a=/=r')    
    
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
    def js(v,x):
        BF11=scsp.jv(v,x)
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
    'cbr is the ratio of c and b terms, will have len(a)-2 terms'
    'rememebr python index is one below actual index as it starts from 0'
    'start from 1 in the for loop and end at index len(a)-1'
    'here i is the index value'
    # print('np.arange(1,len(a),1)',np.arange(1,len(a),1))
    # print('len(a)-2',len(a)-2)
    for i in np.arange(1,len(a),1): 
      if i<=(len(a)-2):
        cbr[i]=np.where(b[i-1]==0,0,c[i-1]/b[i-1])
        # print('cbr in loop',cbr[i])
        'in cprod the previous value of cbr is required' 
        'eg. for B3 (b[1]), F1 is used and reqs. c2/b2==c[0]/b[0]==cbr[1]'
        cprod=F0(a[i]*r[i],cbr[i])*j(1,a[i+1]*r[i])
        c[i]=cprod-sig(i,i+1)*F1(a[i]*r[i],cbr[i])*j(0,a[i+1]*r[i])
        bprod=sig(i,i+1)*F1(a[i]*r[i],cbr[i])*y(0,a[i+1]*r[i])
        b[i]=bprod-F0(a[i]*r[i],cbr[i])*y(1,a[i+1]*r[i])
        'i+1 must not exceed len(a)'
    #print('here is c',c)
    #print('here is b',b)
    'next step after determining constats is as follows'
    'generalise dis terms eg, dis[i]=pi*r[i]*a[i+1]/2'
    'bcn[i]=dis[i-1]*bn[i]*bcn[i-1], B2 needs to be completed outside for loop'
    'next find patter in psi and recalc bcn eg B1,B2,B3 normalized'
    'Psi0 will keep increasing with number of layers chaging B1'
    'rememebr index zero is first element in python'
    'B2=b[0], B3=b[1], ect.'
    'for loop below used to recalculate cbr values as final cbr value outside of'
    'range of above for loop, due to final value of c,b depending on cbr[i-1]'
    #ida = np.where(b!=0)
    #idb = np.where(b==0)
    #cbr[ida] = c[ida]/b[ida]
    #cbr[idb] = 0

    ida = b!=0
    cbr[ida]=c[ida]/b[ida]
    cbr[~ida] = 0
    print('here is cbr',cbr)


    'want to specify the intial field Bz ,Bt'
    'compute all other fields inside for loop'
    bzr=np.zeros((len(a),100))
    btr=np.zeros((len(a),100))
    rma=np.zeros((len(a),100))
    #rma[0]=np.linspace(0,0.5,100)
    'need to make mult dim r'
    'will give error if a[0]=0'
    'generating multi dim array for sub r ranges'
    i = -1
    rma[i+1]=np.linspace(0,r[i+1],100)
    for i in np.arange(len(a)-1):
        rma[i+1]=np.linspace(r[i],r[i+1],100)
    #print('rma=',rma)


    dis=np.zeros(len(a))
    'dis[0]==dis2'
    dis[0]=pi*r[0]*abs(a[1])/2


    bzr[0] =  j(0,a[0]*rma[0])
    btr[0] = sig1(0)*j(1,a[0]*rma[0])
    dis = np.zeros(len(a))
    bn = np.zeros(len(a))
    cn = np.zeros(len(a))
    dis[0] = pi*r[0]*abs(a[1])/2
    bn[0] = dis[0]*b[0] #B2
    cn[0] = dis[0]*c[0] 

    'B2 done bellow porblem'
    for i in np.arange(1,len(a),1):
      if i<=(len(a)-2):
        'for dis[1]==dis3, use a[2]==a3'
        dis[i]=pi*r[i]*abs(a[i+1])/2
        bn[i]=dis[i]*b[i]*bn[i-1]
        cn[i]=dis[i]*c[i]*bn[i-1]
    cn = np.roll(cn,1)
    bn = np.roll(bn,1)

    print('b',b)
    print('c',c)
    for i in np.arange(1,len(a),1):
        'for dis[1]==dis3, use a[2]==a3'
        bzr[i]=bn[i]*j(0,a[i]*rma[i])+cn[i]*y(0,a[i]*rma[i])
        btr[i]=sig1(i)*(bn[i]*j(1,a[i]*rma[i])+cn[i]*y(1,a[i]*rma[i]))
    print('dis',dis)
    #print('bzr[2]',bzr[2])
    #print('btr[2]',btr[2])
    plt.plot(rma[0],bzr[0],color='blue',label='Bz')
    plt.plot(rma[0],btr[0],color='green',label='Bt',linestyle='--')
    for i in np.arange(1,len(a),1):
        plt.plot(rma[i],bzr[i],color='blue')
        plt.plot(rma[i],btr[i],color='green',linestyle='--')
    plt.vlines(r,1,-0.2,color='red',linestyle='-')
    plt.xlabel('r')
    plt.ylabel('B')
    #plt.title('a={}'.format(a))
    plt.legend()
    plt.savefig('bfieldscont')
    plt.show()



    return ida

def aprofile(r,l):
    bz=(1-(l**2)/7+((l**2)/7)*(1-r**2)**7 - l**2*r**2*(1-r**2)**6)**0.5
    alp=2*l*(1-r**2)**2*(1-4*r**2)/bz
    bt=l*r*((1-r**2))**3
    return alp
# r1=np.linspace(0.01,1.01,1)
# a1=aprofile(r1,2)
print('a1',a1)
# a1[-2]=0.01
# a1[-1]=0.5
# a1[-3]=0.1
print('r1',r1)
getb(a1,r1)
