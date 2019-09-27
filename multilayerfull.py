'working program using for multi-layer model over arange of r'
'program can calculate k,w, dw'
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
from matplotlib.ticker import StrMethodFormatter
scsp=scipy.special
pi=np.pi

'Array of helicity "a", "r" used as inputs'
# a1 =[ 3.9995996,   3.93234086,  3.64907622,  2.79882793,  0.95293295, -1.30292591,
#  -2.88752223, -3.46387526, -2.61891399,  0.37831325]
# #r1=[0.25,0.5,0.75,1]
# r1=np.linspace(0.01,1,len(a1))
# r1=r1[1:]

# print('here is a1',a1)
# print('here is r1',r1)

def getdwn(a,r):
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
    #ida = np.where(b!=0)
    #idb = np.where(b==0)
    #cbr[ida] = c[ida]/b[ida]
    #cbr[idb] = 0

    ida = b!=0
    cbr[ida]=c[ida]/b[ida]
    cbr[~ida] = 0
    print('here is cbr',cbr)

    dis=np.zeros(len(a))
    bn=np.zeros(len(a))
    'psi array will be useful for checks'
    psi=np.zeros(len(a))
    'dis[0]==dis2'
    dis[0]=pi*r[0]*abs(a[1])/2
    bn[0]=dis[0]*b[0]

    ctbdiff=j(1,a[0]*r[0])*(1/a[0]-1/a[1])
    'ctbdiff known as ctb1rat in 2lay code'
    psi[0]=2*pi*(r[1]*b[0]*(F1(a[1]*r[1],cbr[0])/(abs(a[1]))*dis[0])+r[0]*ctbdiff)
    
    # print('psi[0]',psi[0])
    'bn term computed outside of loop to allow for'
    'bn[i-1] is the recursion term, cbr[0]==c2/b2'
    'for loop starts at 1'
    for i in np.arange(1,len(a),1):
      if i<=(len(a)-2):
        'for dis[1]==dis3, use a[2]==a3'
        dis[i]=pi*r[i]*abs(a[i+1])/2
        bn[i]=dis[i]*b[i]*bn[i-1]
        psi[i]=(2*pi*bn[i]/abs(a[i+1]))*(r[1+i]*F1(a[i+1]*r[i+1],cbr[i])-r[i]*F1(a[i+1]*r[i],cbr[i]))
    # print('array of psi',psi)
    'len(a)-1 values for psi and psi0 containts 1st and 2nd layer'    
    psit=np.sum(psi)
    # print('total psi',psit)
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
    'if changing (total loop radius) r=!1 multiply arguments in ks by r'
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
    #print('final energy',ws(rootal))
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
    # print('rootal=',rootal)
    return dw

#print(getdwn(a1,r1))

'Hood et al 2009 field profile'
def aprofile(r,l):
    bz=(1-(l**2)/7+((l**2)/7)*(1-r**2)**7 - l**2*r**2*(1-r**2)**6)**0.5
    alp=2*l*(1-r**2)**2*(1-4*r**2)/bz
    bt=l*r*((1-r**2))**3
    return alp
def varelax(num_lays):
    'function to relax Hood bfields'
    rl=np.linspace(0,1.01,num_lays+1)
    rl=rl[1:]
    # print(rl,'rl')
    y=aprofile(rl,2)
    # print(y,'y')
    #y_tmp = copy(y)
    #y=np.where(y>0,y-0.5,y)
    'removing 0 values from array y and coresp. rl values'
    mask0 = y != 0
    rm=rl[mask0]
    y=y[mask0]
    # print('rm',len(rm))
    # print('y',y)
    # plt.bar(rm,y,width=rl[1]-rl[0],edgecolor='red',align='edge')#,alpha=0.5)
    # plt.xlim(rm[0],rm[-1])
    # rshift=0.5*(rm[1]-rm[0])
    # rshift=rm+rshift
    # yshift=0.5*(y[1]-y[0])
    # yshift=y#-yshift
    # plt.scatter(rshift,y)
    # plt.plot(rshift,y,color='black')#,alpha=0.5)
    # plt.xlabel('r')
    # plt.ylabel('a')
    # plt.xlim(rm[0],rm[-1]+(rm[1]-rm[0]))
    # #plt.savefig('hoodalprofile')
    # plt.show()
    # if len(y) == 0: return 0
    # else: return getdwn(y,rm)
    return getdwn(y,rm)
#print(varelax(10))
def repeat_varelax(nl):
	'function to calculate varelax for different # of layers'
    'nl=number of layers, tl= array of layer number'
    nl=20
    tl=np.linspace(1,nl,nl).astype(np.int32)
    #print(tl,'tl')
    #print(len(tl),'tl')
    dwrange=np.zeros(len(tl))
    for i in range(1,len(tl)):
        dwrange[i]=varelax(tl[i])
    dwrange=np.roll(dwrange,-1)
    dwrange[-1]=varelax(len(tl))
    #mask1=dwrange<0.5
    #dwrange=dwrange[mask1]
    #tl=tl[mask1]
    plt.plot(tl,dwrange)
    plt.scatter(tl,dwrange,color='red')
    #plt.ylim(0,dwrange[3])
    plt.xlabel('Layers')
    plt.ylabel('$\delta W$')
    plt.xticks(tl)
    plt.savefig('hooddwrelax')
    #plt.xticks(np.linspace(0,nl,nl/4+1))
    plt.show()
    'code to generate subplot in higher resoltion'
    # fig, ax = plt.subplots(nrows=1, figsize=(7,7))
    # ax1 = fig.add_axes([0.5,0.15,0.4,0.3]) #co-ordinates from the left corner, in figure relative units [x,y,dx,dy]
    # ax.plot(tl,dwrange)
    # ax.scatter(tl,dwrange,color='red')
    # #plt.ylim(0,2)
    # ax.set_xlabel('Layer')
    # ax.set_ylabel('K')
    # ax.set_xticks(tl)
    # #plt.title('Repeated {} layer relaxation model with a0={}'.format(nlays-1 ,a0))
    # #plt.xlim(0,13)
    # tl1=[3,6,9,12,15,18]
    # ax1.set_xlim(3,20)
    # ax1.set_xticks(tl1)
    # ax1.set_ylim(dwrange[3]-0.01,dwrange[-1]+0.01)
    # ax1.axhline(0,color='green',linestyle='--')
    # ax1.plot(tl,dwrange)
    # ax1.scatter(tl,dwrange,color='red')
    # ax1.set_xlabel('Layer')
    # ax1.set_ylabel('K')
    #plt.savefig('hooddwvsnl')
    plt.show()   

#print(repeat_varelax(20)) 

def lindw(gamma,num_lays,a0):
    'function to create linear a vs r profile for n layers'
    alin=np.zeros(num_lays)
    alin[0]=a0
    rlin=np.linspace(0,0.98,num_lays+1)
    'removing 0 from rlin'
    rlin=rlin[1:]
    for i in np.arange(0,num_lays,1):
        alin[i]=alin[0]+(0.5*gamma*rlin[i-1])

    print('rlin',rlin)
    print('alin',alin)
    # plt.bar(rlin,alin,align='edge' ,edgecolor='red')#,alpha=0.5)
    # plt.xlim(rlin[0],rlin[-1])
    # rshift=0.5*(rlin[1]-rlin[0])
    # rshift=rshift+rlin
    # ashift=0.5*(alin[1]-alin[0])
    # ashift=alin-ashift
    # plt.scatter(rshift,alin)
    #print(ashift,'ashift')
    #print('rshift',rshift)
    # plt.plot(rlin,ashift,color='black')#,alpha=0.5)
    # plt.xlabel('r')
    # plt.ylabel('a')
    #plt.title('Linear cylindrical model alpha vs r')
    # plt.show()
    return getdwn(alin,rlin)

#print(lindw(15,10,0.5))

def lindwrep(nl):
	'function to relax lindw for different number of layers'
    'nl=number of layers, tl= array of layer number'
    nl=20
    tl=np.linspace(1,nl,nl).astype(np.int32)
    #print(tl,'tl')
    #print(len(tl),'tl')
    dwrange=np.zeros(len(tl))
    for i in range(1,len(tl)):
        dwrange[i]=lindw(2,tl[i],1)
    #mask1=dwrange<0.5
    #dwrange=dwrange[mask1]
    #tl=tl[mask1]
    dwrange=np.roll(dwrange,-1)
    dwrange[-1]=lindw(2,len(tl),1)
    plt.plot(tl,dwrange)
    plt.scatter(tl,dwrange,color='red')
    #plt.ylim(0,dwrange[3])
    plt.xlabel('Layer')
    plt.ylabel('$\delta W$')
    plt.xticks(tl)
    plt.savefig('lineardwvslay')
    #plt.xticks(np.linspace(0,nl,nl))
    plt.show()
print(lindwrep(20))

def lindwgradvar(ng,nlays,a0):
	'function to relax lindw with varying gradient'
    'ng=# of gradient variations, nlays= ~ lays for each relax, a0=starting point'
    ng=np.linspace(1,ng,15)
    print(ng,'ng')
    dwrange=np.zeros(len(ng))
    for i in range(0,len(ng)):
        dwrange[i]=lindw(ng[i],nlays,a0)
    plt.plot(ng**2,dwrange)
    print(dwrange,'dwrange')
    plt.scatter(ng**2,dwrange,color='red')
    #plt.ylim(0,2)
    plt.xlabel('$\gamma^{2}$')
    plt.ylabel('$\delta W$')

    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}')) # 2 decimal places    #plt.title('Repeated {} layer relaxation model with a0={}'.format(nlays-1 ,a0))
    #plt.xlim(0,13)
    plt.savefig('gammavsdwsq')
    plt.show()
#print(lindwgradvar(1.8,10,1))

'simple code to calculate dw while keeping other alpha values fixed as a check'
# a2=np.linspace(-8,8,99)
# a2w=np.zeros(len(a2))
# for i in range(0,len(a2)):
#     a2w[i]=getdwn([a2[i],a2[i],1,1],r1)
#  #print(a2w)
# plt.plot(a2,a2w)
# plt.ylim(0,2)
# plt.title('2 layer dw vs a1 when a2=1')
# plt.xlabel('a1')
# plt.ylabel('dw')
# #plt.savefig('2laya2=1a1r')
# plt.show()
