'working program using for N layer model over arange of r'
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
a=[1,-2,3,5,-8,-5,-6]
r=[]
'for loop to generate array of r values coresp. to n alpha values'
for i in np.linspace(len(a),1,len(a)):
  ra=1/i
  r.append(ra)

print('here is a',a)
print('here is r',r)

def get_dw_n(a,r):
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
 
    # print('ctbr term in k',ctbr)
    # print('k',k)
    # print('total k,',kt)
    # 'if changing (total loop radius) r=!1 multiply arguments in ks by r'
    def bsq(a):
       B0=np.where(a==0,1/(np.pi),(a/(2*pi*j(1,a)))**2)
       return B0
    # plt.plot(np.linspace(-4,4,88888),bsq(np.linspace(-4,4,88888))) #to avoid zero
    # #plt.plot(np.linspace(-4,4,99999),1/(2+j(1,np.linspace(-4,4,99999))**2))
    # plt.plot(np.linspace(-4,4,88888),j(1,np.linspace(-4,4,88888)))
    # plt.axhline(color='r',linestyle=':')
    # plt.axvline(color='r',linestyle=':')
    # # #plt.title('Wsingle with normalisation a2=0')
    # # #plt.ylabel('Ws')
    # plt.ylim([-1,4])
    # # #plt.savefig('ws41')
    # plt.show()

    def ksr(al):
      'k single layer for root finding, with r as total radius of loop, B,L=1'
      'using normalisation'   
      ks1=(2*pi/al)*bsq(al)*((js(0,al)**2+js(1,al)**2)-(2/al)*js(0,al)*js(1,al))
      return ks1-kt
    
    rootal=scipy.optimize.brentq(ksr, -3.83, 3.83)
    # print('alpha root',rootal)
    #plt.plot(np.arange(-3,3,0.003),ksr(np.arange(-3,3,0.003)))
    #plt.show()
    'first case is seperate, L=1'
    w=np.zeros(len(a))
    mu0=1#4*pi*10**-7
    w0=(r[0]**2)*(j(0,a[0]*r[0])**2+j(1,a[0]*r[0])**2)
    w[0]=((pi/mu0)*bn[0]**2)*(w0-(r[0]/abs(a[0]))*j(0,a[0]*r[0])*j(1,a[0]*r[0]))
    # print('bn[0]',bn[0])
    for i in np.arange(1,len(a),1):
      if i<=(len(a)-1):
        F0a2r2=F0(a[i]*r[i],cbr[i])
        F1a2r2=F1(a[i]*r[i],cbr[i])
        F0a2r1=F0(a[i]*r[i-1],cbr[i])
        F1a2r1=F1(a[i]*r[i-1],cbr[i])
        ws=(r[i]**2)*(F0a2r2**2+F1a2r2**2)-r[i]*F0a2r2*F1a2r2/abs(a[i])-(r[i-1]**2)*(F0a2r1**2+F1a2r1**2)
        w[i]=(ws+r[i-1]*F0a2r1*F1a2r1/abs(a[i]))*((pi/mu0)*bn[i]**2) #here whole term multipied by constant
    # print('work',w)
    # print('intial energy',np.sum(w))
    def ws(al1):
        'w single (relaxed state) for root finding,R=1'
        j00=js(0,al1)
        j11=js(1,al1)  
        ws1=np.where(al1==0,0.5*np.pi*bsq(al1),(np.pi)*bsq(al1)*((j00**2+j11**2)-j00*j11/al1))
        return ws1
    # print('final energy',ws(rootal))
    dw=np.sum(w)-ws(rootal)
    # plt.plot(np.linspace(-9,10,99999),ksr(np.linspace(-9,10,99999)))
    # plt.axhline(color='r',linestyle=':')
    # plt.title("2 layer model, a={0}, dw={1:09.4f}, roota={2:.4f}".format(a, dw, rootal))
    # #plt.xlim([alpha1,alpha2])
    # plt.ylim([-0.5,3])
    # plt.ylabel('k(a)-k(a1,a2)')
    # plt.xlabel('a')
    # #plt.savefig('2krootfinda1={}a2={}'.format(alpha1 ,alpha2))
    # plt.show()

    # print('a=',a)
    # print('r=',r)
    # print('rootal=',rootal)
    # print('Change in work=', dw)
    return dw


    'next step to look for similarities between k1'
    'and ksingle then to find roots'  


#print(get_dw_n(a,r))

#function to generate linear alpha profile wrt r
def lindw(gamma,num_lays,a0):
	#specifify number of alpha and r values
	#need to cut r values into discrete equally spaced chunks between 0-1
	#num-lays=num_alpha
	#ignor r[0] as it's 0
    alin=np.zeros(num_lays)
    alin[0]=a0
    rlin=np.linspace(0,1,num_lays+1)
    rlin=rlin[1:]
    for i in np.arange(1,num_lays,1):
	    alin[i]=alin[0]+0.5*gamma*rlin[i-1]
    #ar[num_lays]=ar[num_lays-1]+0.5*gamma*r[num_lays-1]
    #ar=ar[1:]
    
    # print('rlin',rlin)
    # print('alin',alin)
    # plt.scatter(rlin,alin)
    # plt.show()
    return get_dw_n(alin,rlin)
a=lindw(5,5,1)
print(a)
ng = 60
nlays=15
a0=1
ng=np.arange(1,ng,1)
dwrange=np.zeros(len(ng))
for i in range(1,len(ng)):
    dwrange[i]=lindw(ng[i],nlays,a0)
plt.scatter(ng,dwrange)
plt.ylim(0,50)
plt.xlabel('gamma')
plt.ylabel('dw')
plt.title('repeated {} layer relaxation model with a0={}'.format(nlays ,a0))
plt.xlim(0,15)
plt.show()
    


