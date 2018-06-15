#Import
import numpy as np
import matplotlib.pylab as plt
import random

#Class containing methods to solve the Reaction Diffusion PDE
class PDE(object):

#Init method. Takes time and space intervals, the dimension of the system and initial value of phi
    def __init__(self,dt,dx,dimension,k,sigma,tol,vo):
        self.dt, self.dx = dt, dx
        self.dimension = dimension
        #Order array is the concentration array
        self.order_array = np.full((self.dimension,self.dimension),0.5)
        self.k, self.sigma = float(k),float(sigma)
        self.tol = tol/(self.dimension**2)
        self.vo = vo
    
    #Instance method which imposes periodic boundary conditions on a specific index
    def PBC(self,index):
        return index%(self.dimension)
    
    #Instance method to add random noise to each index of the order parameter array
    def Noise(self):
        for i in range(self.dimension):
            for j in range(self.dimension):
                self.order_array[i,j] = self.order_array[i,j] + random.uniform(-0.1,0.1)


        #Instance method which returns the magnitude of an index from the centre of the lattice
    def Centre_Position(self,i,j,sq):
        x_centre = float(self.dimension/2)
        y_centre = float(self.dimension/2)
        magsq = (i - x_centre)**2 + (j - y_centre)**2
        if sq == True: return magsq
        else: return np.sqrt(magsq)

    #Finite Difference algorithm to update the order parameter array according to Euler Algorithm
    def Update_Order_Array(self,oldarray):
        const = (self.dt)/(self.dx**2.)
        for i in range(self.dimension):
            for j in range(self.dimension):
                #Euler Algorithm
                self.order_array[i,j] = oldarray[i,j] + const*(np.exp(-1*self.Centre_Position(i,j,sq=True)/self.sigma**2) - self.k*oldarray[i,j]) + const*(oldarray[self.PBC(i+1),j]+oldarray[self.PBC(i-1),j]+oldarray[i,self.PBC(j+1)]+oldarray[i,self.PBC(j-1)] - 4*oldarray[i,j]) - self.dt*(self.Adv_Term(i,j))

    #Sweep method which updates the order parameter for the whole lattice. Returns True if successive arrays are within a tolerance ie steady state reached
    def Sweep(self):
        old_array = np.copy(self.order_array)
        self.Update_Order_Array(old_array)
        return np.allclose(old_array,self.order_array,0,self.tol)

    #Instance method to calculate the discretised Laplacian of a given array for a given index
    def Laplacian_2D(self,array,i,j):
        wrtx = (array[self.PBC(i+1),j] + array[self.PBC(i-1),j] -2*array[i,j])/(self.dx**2)
        wrty = (array[i,self.PBC(j+1)] + array[i,self.PBC(j-1)] -2*array[i,j])/(self.dx**2)
        return [wrtx,wrty]

    #Instance method to calculate the discretised Grad squared of a given array and index
    def Grad_Sq(self,array,i,j):
        wrtx = (array[self.PBC(i+1),j] - array[self.PBC(i-1),j])/(2*self.dx)
        wrty = (array[i,self.PBC(j+1)] - array[i,self.PBC(j-1)])/(2*self.dx)
        dot_prod = PDE.Dot_Product([wrtx,wrty],[wrtx,wrty])
        return dot_prod
        
    #Method which calculates the advection term in the PDE. If self.vo is zero then this term does not affect the solution
    def Adv_Term(self,i,j):
        wrtx = (self.order_array[self.PBC(i+1),j] - self.order_array[self.PBC(i-1),j])/(2*self.dx)
        vx = -self.vo*np.sin(2*j*np.pi/self.dimension)
        return vx*wrtx

    #Method to compute the dot product of two vectors (lists)
    def Dot_Product(list1,list2):
        dotlist = []
        for i in range(len(list1)):
            dotlist.append(list1[i]*list2[i])
            return sum(dotlist)

#Method to write the indices, concentration value and distance from centre to a text file. Also returns this data in lists.
    def WriteData(self,name):
        dist_list = []
        conc_list = []
        f = open(str(name)+'.txt','w')
        for i in range(self.dimension):
            for j in range(self.dimension):
                f.write(str(i)+" "+str(j)+" "+str(self.order_array[i,j])+" "+str(self.Centre_Position(i,j,sq=False))+" "+"\n")
                    
                dist_list.append(self.Centre_Position(i,j,sq=False))
                conc_list.append(self.order_array[i,j])
        f.close()
        return np.array(dist_list),np.array(conc_list)

    #Method to create the contour plot of the concentration array
    def PlotArray(self):
        x = y = np.linspace(0,self.dimension-1,self.dimension)
        Y,X = np.meshgrid(x,y)
        plt.contourf(X,Y,self.order_array[:,:],30)
        plt.title("Contour plot of Steady State Concentration Array")
        plt.xlabel("Lattice Position X")
        plt.ylabel("Lattice Position Y")
        plt.colorbar()
        plt.show()









