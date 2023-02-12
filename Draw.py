import numpy as np
import matplotlib.pyplot as plt

f=open('../f.txt','r')
Z=[]
for l in f.readlines():
    v=[float(x) for x in l.split()]
    Z.append(v)
f.close()

Z=np.array(Z)

#plt.imshow(Z, extent=[0, 10, 0, 10], origin='lower', cmap='RdGy')
plt.imshow(Z.T, extent=[0, 20, 0, 10], origin='lower', cmap=plt.cm.jet)
#plt.imshow(Z, extent=[0, 10, 0, 10], origin='lower', cmap=plt.cm.gray)
plt.colorbar()
plt.show()
