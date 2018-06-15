import numpy as np
import sys
import matplotlib.pylab as plt
from matplotlib.animation import FuncAnimation
from Adv_Class import PDE

#Define space and time increments that adhere to Von-Neumann stability.
dx = 1
dt = 0.1

#Create instance of PDE class and apply random noise to concentration array
A = PDE(dt,dx,int(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]))
A.Noise()

#If user selects the task to be visual simulation of the system
if str(sys.argv[6]) == 'viz':
    
    #Animation function
    def UpdatePlot(*args):
        image = ax.imshow(A.order_array)
        A.Sweep()
        return image,
    
    #Create animation
    fig,ax = plt.subplots()
    image = ax.imshow(A.order_array)
    ani = FuncAnimation(fig,UpdatePlot,blit=True)
    plt.title("2D Diffusion Reaction Equation Numerical Solution")
    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    plt.show()

#If user selects the task to be obtaining data relating to the convergence of the solution
elif str(sys.argv[6]) == 'data':
    
    #Loops to calculate the number of iterations to convergence
    n = 0
    steady_state = False
    while steady_state == False:
        n +=1
        if n%500 == 0:
            print(n)
        steady_state = A.Sweep()
    print("Iterations to convergence: " +str(n))

    #Write data to text file and plot the converged array
    r_list,conc_list = A.WriteData('SS_Concentration')
    A.PlotArray()

    #Plot the concentration of the converged array as a function of the distance from the centre
    sorted_r,sorted_conc = zip(*sorted(zip(r_list,conc_list)))
    plt.plot(r_list,conc_list)
    plt.title("Steady State plot of Concentration vs Dist from Centre")
    plt.xlabel("|r|")
    plt.ylabel("Phi")
    plt.show()



