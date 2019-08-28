import numpy as np
import matplotlib.pyplot as plt
import math
import random
from threading import Thread

random.seed()
x=0
y=1
z=2

class RNG:
    gaussian_random = 0.0
    check = False
    def box_muller(self):
        if(rand.check == False):
            u1 = random.random()
            u2 = random.random()
            rand.gaussian_random = math.sqrt(-2.0*math.log(u1))*math.sin(2.0*math.pi*u2)
            rand.check = not rand.check
            return (math.sqrt(-2.0*math.log(u1))*math.cos(2.0*math.pi*u2))
        elif(rand.check == True):
            rand.check = not rand.check
            return rand.gaussian_random

    def tau(self):
        return (-math.log(random.random()))

class gridData:
    def __init__(self):
        self.gridSize = 1.
        self.cellSize = 0.1
        self.radius = 0.8   #radius of sphere of matter

        self.N = int(self.gridSize/self.cellSize)
        n = self.N

        self.density = np.zeros((n,n,n))
        self.temp = np.zeros((n,n,n))
        self.vel = np.zeros((n,n,n,3))
        self.albedo = np.ones((n,n,n))
        self.tauMax = np.ones((n,n,n))*30.
        self.nLeaving = np.zeros((n,n,n))
        self.nEntering = np.zeros((n,n,n))
        self.totalPath = np.zeros((n,n,n))

    def findGridCell(self,point):
        ds = self.cellSize
        i=0
        j=0
        k=0

        for i in range(self.N):
            if(point[x] >= float(i*ds)):
                if(point[x] < float((i+1)*ds)):
                    break

        for j in range(self.N):
            if(point[y] >= float(j*ds)):
                if(point[y] < float((j+1)*ds)):
                    break

        for k in range(self.N):
            if(point[z] >= float(k*ds)):
                if(point[z] < float((k+1)*ds)):
                    break
        return (i,j,k)

    def setupGrid(grid):
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    continue


class packet:
    def __init__(self):
        self.wavelength = 500.e-9
        self.currentCell = []
        self.pos = np.zeros((3))
        self.dir = np.zeros((3))
        self.dtau = 0.
        self.L = 0.
        self.scatters = 0


def emmitPacket(rand,grid,point):
    p = packet()

    for i in range(3):
        p.dir[i] = rand.box_muller()
    norm = math.sqrt(sum(j*j for j in p.dir))
    if (norm == 0.0):
         norm = 1.e-20
    for i in range(3):
        p.dir[i] = p.dir[i] / norm

    for i in range(3):
        p.pos[i] = point[i]

    p.dtau = rand.tau()
    p.L = p.dtau * grid.radius / grid.tauMax[grid.findGridCell(point)]
    return p

def modulus(v):
    return math.sqrt(v[x]**2+v[y]**2+v[z]**2)

def movePacket(grid,p):
    check = True
    currentCell = grid.findGridCell(p.pos)
    vec = p.L * p.dir
    endPos = vec+p.pos

    endCell = grid.findGridCell(endPos)
    s = np.zeros((3))
    if(currentCell != endCell):
        face = 0
        dx = endPos[x] - p.pos[x]
        dy = endPos[y] - p.pos[y]
        dz = endPos[z] - p.pos[z]

        if(dx > 0.0 and (abs(dx) > abs(dy) and abs(dx) > abs(dz))):
            face = (x+1)
            t = (grid.cellSize*float(currentCell[x]+1) - p.pos[x]) / vec[x]
            s[x] = grid.cellSize*float(currentCell[x]+1)
            s[y] = p.pos[y] + (t*vec[y])
            s[z] = p.pos[z] + (t*vec[z])
            # s = (grid.cellSize*(currentCell[x]+1),vec[y]+p.pos[y],vec[z]+p.pos[z])

        elif(dx < 0.0 and (abs(dx) > abs(dy) and abs(dx) > abs(dz))):
            face = -(x+1)
            t = (grid.cellSize*float(currentCell[x]) - p.pos[x]) / vec[x]
            s[x] = grid.cellSize*float(currentCell[x])
            s[y] = p.pos[y] + (t*vec[y])
            s[z] = p.pos[z] + (t*vec[z])
            # s = (grid.cellSize*(currentCell[x]),vec[y]+p.pos[y],vec[z]+p.pos[z])

        elif(dy > 0.0 and (abs(dy) > abs(dx) and abs(dy) > abs(dz))):
            face = (y+1)
            t = (grid.cellSize*float(currentCell[y]+1) - p.pos[y]) / vec[y]
            s[y] = grid.cellSize*float(currentCell[y]+1)
            s[x] = p.pos[x] + (t*vec[x])
            s[z] = p.pos[z] + (t*vec[z])
            # s = (vec[x]+p.pos[x],grid.cellSize*(currentCell[y]+1),vec[z]+p.pos[z])

        elif(dy < 0.0 and (abs(dy) > abs(dx) and abs(dy) > abs(dz))):
            face = -(y+1)
            t = (grid.cellSize*float(currentCell[y]) - p.pos[y]) / vec[y]
            s[y] = grid.cellSize*float(currentCell[y])
            s[x] = p.pos[x] + (t*vec[x])
            s[z] = p.pos[z] + (t*vec[z])
            # s = (vec[x]+p.pos[x],grid.cellSize*(currentCell[y]),vec[z]+p.pos[z])

        elif(dz > 0.0 and (abs(dz) > abs(dx) and abs(dz) > abs(dy))):
            face = (z+1)
            t = (grid.cellSize*float(currentCell[z]+1) - p.pos[z]) / vec[z]
            s[z] = grid.cellSize*float(currentCell[z]+1)
            s[x] = p.pos[x] + (t*vec[x])
            s[y] = p.pos[y] + (t*vec[y])
            # s = (vec[x]+p.pos[x],vec[y]+p.pos[y],grid.cellSize*(currentCell[z]+1))

        elif(dz < 0.0 and (abs(dz) > abs(dx) and abs(dz) > abs(dy))):
            face = -(z+1)
            t = (grid.cellSize*float(currentCell[z]) - p.pos[z]) / vec[z]
            s[z] = grid.cellSize*float(currentCell[z])
            s[x] = p.pos[x] + (t*vec[x])
            s[y] = p.pos[y] + (t*vec[y])
            # s = (vec[x]+p.pos[x],vec[y]+p.pos[y],grid.cellSize*(currentCell[z]))

        p.L = modulus(s-p.pos)
        p.pos = s
        tempCell = np.zeros((3))
        for i in range(3):
            tempCell[i] = currentCell[i]

        tempCell[abs(face)-1] += (face/abs(face))
        p.currentCell = tempCell
        p.scatters += 1

        grid.totalPath[currentCell] += p.L
        grid.nLeaving[currentCell] += 1
        grid.nEntering[currentCell] += 1
    else:
        p.pos += vec
        p.L = modulus(vec)
        p.currentCell = currentCell
        p.scatters += 1
        grid.totalPath += p.L


    for i in range(3):
        if(int(p.currentCell[i]) > (grid.N-1)):
            check = False
            break
        elif(int(p.currentCell[i]) < 0):
            check = False
            break
    return check

def scatterPacket(grid,p):
    while(1):
        if(not movePacket(grid,p)):
            print("packet left grid")
            print("num or scatters =",p.scatters)
            break
        currentCell = (int(p.currentCell[x]),int(p.currentCell[y]),int(p.currentCell[z]))
        a = grid.albedo[currentCell]
        print(currentCell)
        if(random.random() <= a):
            print("scattered")
            for i in range(3):
                p.dir[i] = rand.box_muller()
            norm = math.sqrt(sum(j*j for j in p.dir))
            p.dtau = rand.tau()
        else:
            print("absorbed")




rand = RNG()
grid = gridData()
startPoint = ((0.55,0.55,0.55))
p = emmitPacket(rand,grid,startPoint)
scatterPacket(grid,p)
