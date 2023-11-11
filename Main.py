
import numpy as np
import os
import time 
import math as mt
import matplotlib.pyplot as plt
from Modules.Momentum import Momentum
from Modules.Pressure import Pressure_adjust
from Modules.pressuredemo import Pressure_adjust1
from Update.Pressure_update import add
from Modules.convergence import convergence
#from Plots.Pressure_plot import plot_pressure
from Plots.Plot import Graph_PV,plot
from Plots.Velocity_plot import Graph_V,plot_V
from Plots.Pressure_plot import  Graph_P,plot_P 


"""
const paramers through-out the analysis
All units are taken in SI system (Kilogram,Meter,Second)

"""
Total_time=1
dt = 0.00001                              # Time step size
#n = int(Total_time/dt)                   # Number of time steps
n=1
rho = 996.550                          # Density kg/m3
g = 9.8                                # Gravitational acceleration m/s2
Grid_points=500                   #Total number of grid points (Velocity faces)
Inletmassflux=900                      #kg/m2.s
P_atm=1                                #atm
dia=0.0154                             #m
d_vis=0.000854                         #Ns/m2

p_exit= P_atm* 101325                  #N/m2
A= (mt.pi)*(dia**2)/4                  #m2
u_inlet= Inletmassflux/(rho)           #m/s
length=2                               #m
dx = length/Grid_points                # Spatial grid size


"""
Variables evolve-with time
"""
# Initialize arrays for velocity (u) and pressure (p) at staggered grid locations

u_n = np.zeros(Grid_points)                            # Velocity at faces
p_star = np.zeros(Grid_points)                         # pressure at cell centers
A_n = np.full((2 * Grid_points + 1), A)                # Area at n_th step
perimeter= (mt.pi)*dia                                 #perimeter at n_th step 
p_s=np.full((2 * Grid_points + 1), perimeter)          #perimeter at n_th step
dp= 128*d_vis*length*A*u_inlet/((mt.pi)*dia**4)        #pressure change for end points(Change the length to get the desired pressure drop)
Area=np.full((2 * Grid_points + 1), A)


"""
Initialize the pressure 
"""
p_star[len(p_star)-1]=p_exit
for i in range(len(p_star)-2,-1,-1):
    p_star[i]=p_star[i+1]+((dp*dx)/length)




"""
Initialize the velocity 
"""
u_n[0]=u_inlet
for i in range(1,len(u_n)):
    u_n[i]=0.01    

graph=Graph_PV()
#graph_v=Graph_V(u_n)
#graph_p=Graph_P(p_star)
x=np.arange(0,len(u_n))
#plot_V(u_n,x,graph_v)
#plot_P(p_star,x,graph_p)

# Main Algo
def unsteady_1D_flow(A_n,Area,u_n,p_star,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit):
    
    for t in range(0,n):                                                                                        
        u_star=Momentum(p_star,u_n,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
        if np.array_equal(u_star,u_n):
            print("Yes, the arrays are the same.")
        else:
            print("No, the arrays are not the same.")
            #Get u_star for guessed pressure
        converge=convergence(u_star,Grid_points,rho,Area,A_n,dx,dt)
        converge1=convergence(u_n,Grid_points,rho,Area,A_n,dx,dt)                             #Check the convergence
        #graph=Graph_PV(u_star,p_star)
        x=np.arange(0,len(u_star))
        #plot(u_star,p_star,x,graph)
        #graph_v=Graph_V(u_star)
        #graph_p=Graph_P(p_star)
        start_time = time.time()
        print(converge,converge1)
        i=0
        while(i<0):
            p_add=Pressure_adjust(p_star,u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
            p_star=add(p_star,p_add,Grid_points)                                                #update the pressure
            u_star=Momentum(p_star,u_n,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
            converge=convergence(u_star,Grid_points,rho,Area,A_n,dx,dt)
            #debug
            plot(u_star,p_star,x,graph)
            #plot_V(u_star,x,graph_v)
            #plot_P(p_star,x,graph_p)
            plt.pause(0.01)
            #debug
            i=i+1


            elasp_time = time.time() - start_time
            if elasp_time >4:
                break

        # Save the current figure with a unique file name
        # # Unique file name based on iteration
        #plt.savefig(file_name)
        #plt.clf()  # Clear the current figure to prepare for the next iteration
        plt.show()    
      
        u_n=u_star

    return u_star,p_star



"""
Final result 
"""

velocity, Pressure =unsteady_1D_flow(A_n,A,u_n,p_star,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit)
