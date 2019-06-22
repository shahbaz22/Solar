import numpy as np
import matplotlib.pyplot as plt
an=[3,2,4,5,6,4,7,6,5]
an1=[]
count=[]
for i in np.arange(0,len(an)):
  if i==0:
    kn=2*4
    an1.append(kn)

  elif i<(len(an)-1):
    kn=2*an[i]
    an1.append(kn)
  else:
  	kn=0
  	an1.append(kn)
print(an1)
print(len(an))

for ind, el in enumerate(an):
  if ind ==0:
    an1.append(2*4)
    count.append(ind)
  elif ind == 2:
    an1.append(el*2)
    count.append(ind)
  else:
    an1.append(0)
    count.append(ind)


print(an1,count)
print(len(an))

#print(np.arange(0,len(an)))