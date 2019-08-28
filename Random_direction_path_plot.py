import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import random

random.seed()
dimensions = 3
num_photons = 1
num_paths = 100


def add_column(v):
    size = np.shape(v)
    new = np.zeros((num_photons,dimensions,size[2]+1))
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                new[i,j,k+1] = v[i,j,k]
    return new


v=np.zeros((num_photons,dimensions,num_paths+1))
for k in range(num_paths+1):
    for i in range(dimensions):
        for j in range(num_photons):
            v[j,i,k] = random.random()

v=add_column(v)

fig=plt.figure()

if(dimensions == 3):
    ax=fig.add_subplot(111,projection='3d')
    for i in range(num_photons):
            plt.plot(v[i,0,:],v[i,1,:],v[i,2,:])
elif(dimensions > 3 or dimensions <2):
    print("too many dimensions (max 3, min 2)")
    exit()
else:
    for i in range(num_photons):
        plt.plot(v[i,0,:],v[i,1,:])


plt.show()
