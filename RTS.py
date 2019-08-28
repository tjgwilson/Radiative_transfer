import numpy as np
import math
import random
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

random.seed()
num_photons = 10000
radius = 1
dim = 3

tau_max = 2
tau_min = 2


show_path = False
find_average_scatters = True
scatter_loops = 100

class rand:

    def __init__(self,dimensions):
        self.dimensions=dimensions

    def box_muller():
        u1 = random.random()
        u2 = random.random()
        return (math.sqrt(-2.0*math.log(u1))*math.cos(2.0*math.pi*u2))

    def vector(self):
        vector = np.zeros((self.dimensions))

        for i in range(self.dimensions):
            vector[i] = rand.box_muller()

        norm = math.sqrt(sum(j*j for j in vector))

        for i in range(self.dimensions):
            vector[i] = vector[i] / norm
        return vector

    def polar_vector(self):
        vector = rand.vector(self)
        sph_vect = np.zeros((self.dimensions))

        if(self.dimensions == 3):
            sph_vect[0] = math.sqrt(sum(i*i for i in vector))
            sph_vect[1] = math.atan(vector[1] / vector[0])
            sph_vect[2] = math.acos(vector[2] / sph_vect[0])
        elif(self.dimensions == 2):
            sph_vect[0] = math.sqrt(sum(i*i for i in vector))
            sph_vect[1] = math.atan(vector[1] / vector[0])
        elif (self.dimensions == 1):
            sph_vect[0] = vector[0]
        else:
            print("ERROR: Unable to create polar coordinates with inputed dimensions. Terminating simulations")
            exit()
        return sph_vect

    def optical_depth():
        return (-math.log(1-random.random()))

    def path_length(self,depth,tau_max):
        return ((rand.optical_depth() * depth) / tau_max)






three_dim=rand(dim)
av_scat_array = np.zeros((scatter_loops,2))

"""
if (tau_max < tau_min):
    x =tau_max
    tau_max = tau_min
"""
if(find_average_scatters):
    show_path = False
elif(show_path):
    scatter_loops = 1
    path = np.zeros(())

t = tau_min
t_step = (tau_max - tau_min) / scatter_loops

for k in range(scatter_loops):
    vects = np.zeros((num_photons,dim))
    num_scatters=np.zeros((num_photons))
    if(show_path):
        fig=plt.figure(1)
        ax=fig.add_subplot(111,projection='3d')
    for i in range(num_photons):
        distance = 0.0
        if(show_path):
            xo = 0.0
            yo = 0.0
            zo = 0.0
        while(distance < radius):
            vects[i,:] += three_dim.vector() * three_dim.path_length(radius,t)
            distance = math.sqrt(sum(j*j for j in vects[i,:]))
            num_scatters[i] += 1

            if(show_path):
                xi = [xo,vects[i,0]]
                yi = [yo,vects[i,1]]
                zi = [zo,vects[i,2]]
                xo = vects[i,0]
                yo = vects[i,1]
                zo = vects[i,2]
                ax.plot(xi,yi,zi)

    for i in range(num_photons):
        av_scat_array[k,0] += num_scatters[i]
    av_scat_array[k,0] = av_scat_array[k,0] / num_photons
    av_scat_array[k,1] = t
    t += t_step
    if(show_path):
            ax.set_xlim(-radius,radius)
            ax.set_ylim(-radius,radius)
            ax.set_zlim(-radius,radius)
            ax.scatter(0.0,0.0,0.0,"o",s=50,c='r')
            u,v = np.mgrid[0:2.0*np.pi:20j, 0:np.pi:10j]
            x = radius * np.cos(u)*np.sin(v)
            y = radius * np.sin(u)*np.sin(v)
            z = radius * np.cos(v)
            ax.plot_wireframe(x,y,z,color="k",alpha=0.3)
            plt.show()


print(av_scat_array[:,1])

plot_num = 2
if(t_step > 0.0):
    plot_num =1


if(find_average_scatters):
    fig=plt.figure(2)
    if(t_step <= 0.0):

        fig.add_subplot(plot_num,1,1)
        plt.hist(x=av_scat_array[:,0],bins='auto',density=True)
        plt.axvline(av_scat_array[:,0].mean(), color='r', linestyle='dashed', linewidth=1,label="Average = {}" .format(av_scat_array[:,0].mean()))
        plt.xlabel("Mean number of scatters")
        plt.legend(title='Expected  scatter = {}'.format((tau_max*tau_max)/2))

        fig.add_subplot(plot_num,1,2)
        plt.plot(av_scat_array[:,0],'k',label='Observed mean scatters')
        plt.plot([1,scatter_loops],[(tau_max*tau_max)/2,(tau_max*tau_max)/2],'r',label='Expected mean scatters')
        plt.legend()
        plt.xlabel("Loop No.")
        plt.ylabel("Mean scatters")
    elif(t_step > 0.0):
        fig.add_subplot(plot_num,1,1)
        plt.plot(av_scat_array[:,1],av_scat_array[:,0],'k',label='Observed mean scatters')
        plt.plot(av_scat_array[:,1],(av_scat_array[:,1]*av_scat_array[:,1])/2,'r',label='Expected mean scatters')
        plt.legend()
        plt.xlabel("Max Optical Depth")
        plt.ylabel("Mean scatters")
    plt.show()
