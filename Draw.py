# https://juejin.cn/post/7034117434659471397 vector field
import numpy as np
import matplotlib.pyplot as plt

f=open('../f.txt','r')
Z=[]
for l in f.readlines():
    v=[float(x) for x in l.split()]
    Z.append(v)
f.close()

Z=np.array(Z)

plt.imshow(Z.T, extent=[0, 10e-2, 0, 5e-2], origin='lower', cmap=plt.cm.jet)
#plt.imshow(Z, extent=[0, 10, 0, 10], origin='lower', cmap=plt.cm.gray)
plt.colorbar()
plt.show()
