import numpy as np
import random
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

class rand:
    gaussian_random = 0.0
    check = False
    def __init__(self,dimensions):
        self.dimensions=dimensions

    def box_muller():
        if(rand.check == False):
            u1 = random.random()
            u2 = random.random()
            rand.gaussian_random = math.sqrt(-2.0*math.log(u1))*math.sin(2.0*math.pi*u2)
            rand.check = not rand.check
            return (math.sqrt(-2.0*math.log(u1))*math.cos(2.0*math.pi*u2))
        elif(rand.check == True):
            rand.check = not rand.check
            return rand.gaussian_random

    #Uses guassian distributed vectors to produce a N dimensional cartesian vector (normalized) uniformal distributed in all dimensions
    def Ndim_vector(self):
        vector = np.zeros((self.dimensions))

        for i in range(self.dimensions):
            vector[i] = rand.box_muller()

        norm = math.sqrt(sum(j*j for j in vector))

        for i in range(self.dimensions):
            vector[i] = vector[i] / norm
        return vector

    #Uses guassain distributed cartesian vectors to create polar vector in 1,2,3 dimensions.
    def vector_polar(self):
        vector = rand.Ndim_vector(self)
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
            print("ERROR: Unable to create polar coordinates with input dimensions. Terminating simulations")
            exit()
        return sph_vect

    #Creates 3D vector that is not uniformally distriubed in all dimensions (used for testing purposes). Spherical cooridinate output
    def vector_simple(self):
        vector = np.zeros((3))
        vector[0] = 2.0 * np.pi * random.random()
        vector[1] = np.pi * random.random()
        if(vector[0] >= np.pi):
            vector[0] = (-np.pi)+(vector[0] - np.pi)
        vector[2] = 1.0

        return vector

    #uses cummulative probability theory to create a uniformally distriubted 3D random vector (normalized). Spherical cooridinate output
    def vector_cpf(self):
        vector = np.zeros((3))
        vector[0] = 2.0 * np.pi * random.random()
        vector[1] = math.acos(1.0 - (2.0*random.random()))
        if(vector[0] >= np.pi):
            vector[0] = (-np.pi)+(vector[0] - np.pi)
        vector[2] = 1.0

        return vector

    #Plots 3D plot of normalized vectors for specified type of RNG
    def plot_vecs_example(self, num_vecs, type):

        fig=plt.figure()
        ax=fig.add_subplot(111,projection='3d')
        origin = np.zeros((3))
        point = np.zeros((3))

        if(self.dimensions != 3):
            print("Incorrect number of dimensions")
            exit()

        for i in range(num_vecs):
            if(type == "simple"):
                vec = rand.vector_simple(self)
                point[0] = vec[2] * math.sin(vec[1]) * math.cos(vec[0])
                point[1] = vec[2] * math.sin(vec[1]) * math.sin(vec[0])
                point[2] = vec[2] * math.cos(vec[1])
            elif(type == "gaussian"):
                point = rand.Ndim_vector(self)
            elif(type == "cpf"):
                vec = rand.vector_cpf(self)
                point[0] = vec[2] * math.sin(vec[1]) * math.cos(vec[0])
                point[1] = vec[2] * math.sin(vec[1]) * math.sin(vec[0])
                point[2] = vec[2] * math.cos(vec[1])

            c = np.c_[origin,point]
            plt.plot(c[0,:],c[1,:],c[2,:])
        return fig

    #Some routines to visualize distribution of random numbers from different RNG's to see how good they are.
    def test_random(self, num_photons, type):
        data = np.zeros((num_photons,2))
        num_bins = 200
        threshold = 1

        if(type == "gaussian"):
            title = "Gaussian Produced Vectors"
            for i in range(num_photons):
                vect = rand.Ndim_vector(self)
                data[i,0] = math.atan2(vect[1],vect[0])
                data[i,1] = math.acos(vect[2])
        elif(type == "cpf"):
            title = "CPF Produced Vectors"
            for i in range(num_photons):
                vect = rand.vector_cpf(self)
                if(vect[2] == 1):
                    data[i,0] = vect[0]
                    data[i,1] = vect[1]
                else:
                    exit()
        elif(type == "simple"):
                title = "Incorrectly Produced Vectors"
                for i in range(num_photons):
                    vect = rand.vector_simple(self)
                    if(vect[2] == 1):
                        data[i,0] = vect[0]
                        data[i,1] = vect[1]
                    else:
                        exit()

        dphi = 2.0*np.pi / num_bins
        dtheta = np.pi / num_bins

        binsphi = np.linspace((-np.pi),(np.pi),num_bins,endpoint=True)
        binstheta = np.linspace(0.0,(np.pi),num_bins,endpoint=True)
        h, binsphi, binstheta = np.histogram2d(data[:,0],data[:,1],bins=[binsphi,binstheta])

        for i in range(num_bins-1):
            for j in range(num_bins-1):
                h[i,j] = h[i,j] / ((dtheta*dphi*math.sin(binstheta[i])))


        count = 0
        rand_vec = np.zeros((num_bins**2,3))
        for i in range(num_bins-1):
            for j in range(num_bins-1):
                rand_vec[count,0] = binsphi[i]
                rand_vec[count,1] = binsphi[j]
                rand_vec[count,2] = h[i,j]
                count += 1


        fig = plt.figure()
        ax=fig.add_subplot(223)
        im = ax.imshow(h, interpolation='none', origin='low',norm=LogNorm(), extent=[binsphi[0], binsphi[-1], binstheta[0], binstheta[-1]])
        plt.title("Density normalized to grid size")
        plt.xlabel("Phi [rad]")
        plt.ylabel("Theta [rad]")
        fig.colorbar(im,orientation="vertical")
        ax=fig.add_subplot(221)
        hst1=ax.hist(x=data[:,0],bins=num_bins,density=True)
        plt.title("phi")
        ax=fig.add_subplot(222)
        hst2=ax.hist(x=data[:,1],bins=num_bins,density=True)
        plt.plot(binstheta,0.5*np.sin(binstheta),'r',label="sin(theta)")
        plt.legend()
        plt.title("theta")

        fig.suptitle(title)
        return fig


############################
############################
three  = rand(3)
#a=three.test_random(100000,"simple")
#b=three.test_random(100000,"gaussian")
#c=three.test_random(100000,"cpf")

#three.plot_vecs_example(100, "gaussian")
three.plot_vecs_example(5, "cpf")
plt.show()
############################
############################
