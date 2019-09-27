'working program to calculate the transfer of helicity using multi-layer code'
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

'Array of "a", "r" used as inputs'


def get_dk_n(a,r):
    'function to get n values of helicty for n layers'
    'r is array for radii such that a[1]->r[1]...a[n]->r[n]'
    'if cbrat=0 then j obtained'
    'cbr determins order of function, eg. G0,G1 when c3/b3 '
    'check for len(a)==len(r)'
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
        #print('cbr in loop',cbr[i])
        'in cprod the previous value of cbr is required' 
        'eg. for B3 (b[1]), F1 is used and reqs. c2/b2==c[0]/b[0]==cbr[1]'
        cprod=F0(a[i]*r[i],cbr[i])*j(1,a[i+1]*r[i])
        c[i]=cprod-sig(i,i+1)*F1(a[i]*r[i],cbr[i])*j(0,a[i+1]*r[i])
        bprod=sig(i,i+1)*F1(a[i]*r[i],cbr[i])*y(0,a[i+1]*r[i])
        b[i]=bprod-F0(a[i]*r[i],cbr[i])*y(1,a[i+1]*r[i])
        'i+1 must not exceed len(a)'
    # print('here is c',c)
    # print('here is b',b)
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
    # print('here is cbr',cbr)

    dis=np.zeros(len(a))
    bn=np.zeros(len(a))
    'psi array will be useful for checks'
    psi=np.zeros(len(a))
    'dis[0]==dis2'
    dis[0]=pi*r[0]*abs(a[1])/2
    bn[0]=dis[0]*b[0]

    ctbdiff=j(1,a[0]*r[0])*(1/abs(a[0])-sig(0,1)/abs(a[1]))
    # print(ctbdiff,'ctbdiff')
    'ctbdiff known as ctb1rat in 2lay code'
    psi[0]=2*pi*(r[1]*b[0]*F1(a[1]*r[1],cbr[0])*dis[0]/((abs(a[1])))+r[0]*ctbdiff)
    
    # print('psi[0]',psi[0])
    'bn term computed outside of loop to allow for'
    'bn[i-1] is the recursion term, cbr[0]==c2/b2'
    'for loop starts at 1'
    for i in np.arange(1,len(a),1):
      if i<=(len(a)-2):
        'for dis[1]==dis3, use a[2]==a3'
        dis[i]=pi*r[i]*abs(a[i+1])/2
        bn[i]=dis[i]*b[i]*bn[i-1]
        psi[i]=(2*pi*bn[i]/abs(a[i+1]))*(r[1+i]*F1(a[i+1]*r[i+1],cbr[i])-sig(i,i+1)*r[i]*F1(a[i+1]*r[i],cbr[i]))
    #print('array of psi',psi)
    'len(a)-1 values for psi and psi0 containts 1st and 2nd layer'    
    psit=np.sum(psi)
    #print('total psi',psit)
    b1=1/psit
    # print('here is b1',b1)
    'element-wise multiplication using np.dot'
    bn=np.dot(b1,bn)
    # print('here is bn',bn)
    'next step to compare values for 3layer for bn and b1~psi'
    'next section is the helicity first values of helicity ca be calculated outside for loop'
    k=np.zeros(len(a))
    'k1'
    k[0]=(r[0]**2)*(j(0,a[0]*r[0])**2+j(1,a[0]*r[0])**2)
    k[0]=sig1(0)*((2*pi*b1**2)/abs(a[0]))*(k[0]-2*(r[0]/abs(a[0]))*j(0,a[0]*r[0])*j(1,a[0]*r[0]))
    'if there is a print error it is due to a syntax error in the called equation'
    # print('k[0]',k[0])
    'reason for k2==k[1] perfomed outside the for loop is because the k2 term is not general'
    'k2 term requires the first value of ctbrat and in general the new value and the previous value ctbrat required'
    'cannot start loop from a lower value than 1'
    ctbr=np.zeros(len(a))

    'np.roll used to shift inex of bn one to the right and add b1'
    bn=np.roll(bn,1)
    # print('here is shifted bn',bn)
    bn[0]=b1
    # print('here is complete bn',bn)
    'shifting indices of cbr to start the loop from 1!!!!'
    cbr=np.roll(cbr,1)
    # print('cbr',cbr)

    for i in np.arange(1,len(a),1):
      if i<=(len(a)-1):
          F0a2r2=F0(a[i]*r[i],cbr[i])
          F1a2r2=F1(a[i]*r[i],cbr[i])
          F0a2r1=F0(a[i]*r[i-1],cbr[i])
          F1a2r1=F1(a[i]*r[i-1],cbr[i])

          ctbr[i]=bn[i-1]*r[i-1]*F1(a[i-1]*r[i-1],cbr[i-1])*(1/abs(a[i-1])-sig(i-1,i)/abs(a[i]))+ctbr[i-1]    
          ks=(r[i]**2)*(F0a2r2**2+F1a2r2**2)-2*r[i]*F0a2r2*F1a2r2/abs(a[i])-(r[i-1]**2)*(F0a2r1**2+F1a2r1**2)
          ks=(ks+2*r[i-1]*F0a2r1*F1a2r1/abs(a[i]))*(2*pi*bn[i]**2)/abs(a[i]) #here whole term multipied by constant
          ks=ks+(4*pi*bn[i]*ctbr[i]*(F0a2r1-F0a2r2))/abs(a[i])
          k[i]=ks*sig1(i)
    kt=np.sum(k)

    def bsq(a):
       B0=np.where(a==0,1/(np.pi),(a/(2*pi*j(1,a)))**2)
       return B0

    def ksr(al):
      'k single layer for root finding, with r as total radius of loop, B,L=1'
      'using normalisation'   
      ks1=(2*pi/al)*bsq(al)*((js(0,al)**2+js(1,al)**2)-(2/al)*js(0,al)*js(1,al))
      return ks1-kt
    rootal=scipy.optimize.brentq(ksr, -3.83, 3.83)
    
    # rshift=0.5*(r[1]-r[0])
    # rshift=-rshift+r

    # plt.scatter(rshift,k-ksr(rootal))
    # plt.plot(rshift,k-ksr(rootal))
    # plt.title("{0} layer model, a={1}, dk={2:09.4f}, roota={3:.4f}".format(len(a),a, kt-ksr(rootal), rootal))

    # plt.vlines(r,min(k),max(k),linestyle='--',color='red')
    # plt.vlines(0,min(k),max(k),linestyle='--',color='red')
    # plt.hlines(0,0,1,linestyle='--')
    # plt.xlabel('r')
    # plt.ylabel('k(a,r)-k(roota)')
    # plt.show()
    #print('k',k)
    # print('ksr',ksr(rootal))
    #print(rootal,'rootal')
    'helicity returned with root a value append (last element) for later use'
    return np.append(k,rootal)


def lindw(gamma,num_lays,a0):
    'function to create linear a vs r profile for n layers'
    'linear alpha function for intial values'
    alin=np.zeros(num_lays)
    alin[0]=a0
    rlin=np.linspace(0,1,num_lays+1)
    'removing 0 from rlin'
    rlin=rlin[1:]
    for i in np.arange(1,num_lays,1):
        alin[i]=alin[0]+(0.5*gamma*rlin[i-1])

    #print('rlin',rlin)
    #print('alin',alin)
    'ploting script'
    # plt.bar(rlin,alin,align='edge' ,edgecolor='red')#,alpha=0.5)
    # plt.xlim(rlin[0],rlin[-1])
    # rshift=0.5*(rlin[1]-rlin[0])
    # rshift=rshift+rlin
    # ashift=0.5*(alin[1]-alin[0])
    # ashift=alin-ashift
    # plt.scatter(rshift,alin)
    # #print(ashift,'ashift')
    # #print('rshift',rshift)
    # plt.plot(rlin,ashift,color='black')#,alpha=0.5)
    # plt.xlabel('r')
    # plt.ylabel('a')
    # plt.title('Linear cylindrical model alpha vs r')
    # plt.show()
    return alin

'Code to subtract average relaxed alpha from linear alpha profile to obtain helicity tranfer'    
nl,g=10,2.5

'al is alpha array with dim=# of layers , g is the gradent (steepness of alpha), nl is the # of layers'
al=lindw(g,nl,1)
print('al',al)
'defining rl'
rl=np.linspace(0,1,len(al)+1)
'removing 0 from rlin'
rl=rl[1:]
k1=get_dk_n(al,rl)
print('k1',k1)
'last value for k1 is roota'
roota=k1[-1]
print('roota',roota)
'array of roota of dimension of the initial alpha array'
arootar=np.full_like(al,roota,dtype=float) #array full of roota with dim of al
print('arootar',arootar)
k2=get_dk_n(arootar,rl)
print(k2)
'removing root alpha (the last element) from k1'
k1=k1[0:len(al)]
k2=k2[0:len(al)]
print('k1',k1)
print('k2',k2)
's checks whether dk=0, as it should'
s=np.sum(k1-k2)
print('s',s)
rshift=0.5*(rl[1]-rl[0])
rshift=-rshift+rl
'k1-k2 gives the change in helicty in each layer'
'k(alpha)-k(alpha root) for each layer'
plt.scatter(rshift,k1-k2)
plt.plot(rshift,k1-k2)
#plt.title("{0} layer model, a={1}, dk={2:09.4f}, roota={3:.4f}".format(len(a),a, kt-ksr(rootal), rootal))
#plt.savefig('gamma{}nl{}ktran'.format(g,nl))
plt.vlines(rl,min(k1-k2),max(k1-k2),linestyle='--',color='red')
plt.ylim(min(k1-k2),max(k1-k2))
plt.vlines(0,min(k1-k2),max(k1-k2),linestyle='--',color='red')
plt.hlines(0,0,1,linestyle='--',color='green')
plt.xlabel('r')
plt.ylabel('$\delta K$')
plt.show()
    