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
a1,a2,a3,a4=3,3,2,-1
'Array of helicity "an" used as inputs'

an=[a1,a2,a3,a4]
def get_helicityn(a):
	'function to get n values of helicty for n layers'
	'rn is array for radii such that a[1]->rn[1]...a[n]->rn[n]'
	'''rn=[]
	'for loop to generate array of r values coresp. to n alpha values'
	for i in np.linspace(len(an),1,len(an)):
		r=1/i
		rn.append(r)'''
	rn=[0.25,0.5,0.75,1]
	print(rn)
	#print(rn[0],rn[2])
	#print(len(an))
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
	'if some args. always +1 then always have diff order'
	'hence some functions always remain same'
	'increasing order F0->G0... ,F1->G1 wrt 3 lay code'
	cn=np.zeros(len(an))
	bn=np.zeros(len(an))
	cbrat=np.zeros(len(an))
	c2prod=j(0,an[0]*rn[0])*j(1,an[1]*rn[0])-sig(0,1)*j(1,an[0]*rn[0])*j(0,an[1]*rn[0])
	b2prod=sig(0,1)*j(1,an[0]*rn[0])*y(0,an[1]*rn[0])-j(0,an[0]*rn[0])*y(1,an[1]*rn[0])
	cn[0]=c2prod
	bn[0]=b2prod
	cbrat[0]=0
	for i in range(len(an)):
		if i<=(len(an)-2):
			cbrat[i]=np.where(bn[i-1]==0,0,cn[i-1]/bn[i-1])
			print(cbrat)
			cprod=F0(an[i]*rn[i],cbrat[i])*j(1,an[i+1]*rn[i])
			cprod=cprod-sig(i,i+1)*F1(an[i]*rn[i],cbrat[i])*j(0,an[i+1]*rn[i])
			print(cprod)
			bprod=sig(i,i+1)*F1(an[i]*rn[i],cbrat[i])*y(0,an[i+1]*rn[i])
			bprod=bprod-F0(an[i]*rn[i],cbrat[i])*y(1,an[i+1]*rn[i])
			print(bprod)
			cn[i]=cprod
			bn[i]=bprod
		'i+1 must not exceed len(an)'
	#print(len)
	return	'here is cn',cn,'here is bn',bn
'next step would be to check values for constants with 2lay and 3 lay code'
'need to fix program for when a vals negative'
print(get_helicityn(an))

