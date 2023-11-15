
import numpy as np
import os
import time 
import math as mt
import matplotlib.pyplot as plt
from Modules.Momentum import Momentum
from Modules.Pressure import Pressure_adjust
from Update.Pressure_update import add
from Modules.convergence import convergence
from Modules.Pressure import Pressure_adjust
from Plots.Pressure_plot import plot_pressure
from Plots.Plot import Graph_PV,plot
from Plots.Velocity_plot import Graph_V,plot_V
from Plots.Pressure_plot import  Graph_P,plot_P 


"""
const paramers through-out the analysis
All units are taken in SI system (Kilogram,Meter,Second)

"""
Total_time=1
#dt = 0.001                              # Time step size
#n = int(Total_time/dt)                  # Number of time steps
n=1
rho = 996.550                          # Density kg/m3
g = 9.8                                # Gravitational acceleration m/s2
Grid_points=200                     #Total number of grid points (Cell centres including extra cell at the end)
Inletmassflux=900                      #kg/m2.s
P_atm=1                                #atm
dia=0.0154                             #m
d_vis=0.000854                         #Ns/m2
p_exit= P_atm* 101325                  #N/m2
A= (mt.pi)*(dia**2)/4                  #m2
u_inlet= Inletmassflux/(rho)           #m/s
length=2                               #m
dx = length/Grid_points 
dt=0.0001               # Spatial grid size
CFL=u_inlet*dt/dx
minimum_dt=dx/u_inlet     #<0.0044


"""
Variables evolve-with time
"""
A_n = np.full((2 * Grid_points + 1), A)                # Area at n_th step
perimeter= (mt.pi)*dia                                 #perimeter at n_th step 
p_s=np.full((2 * Grid_points + 1), perimeter)          #perimeter at n_th step
dp= 128*d_vis*length*A*u_inlet/((mt.pi)*dia**4)        #pressure change for end points(Change the length to get the desired pressure drop)
print(dp,"dp")

"""
Initialize the velocity 
"""
u_n = np.zeros(Grid_points+1)  
u_n[0]=u_inlet
for i in range(1,len(u_n)):
    u_n[i]=0.01    


# Main Algo
def unsteady_1D_flow(A_n,A,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit):
    p_star = np.zeros(Grid_points)
    u_star = np.zeros(Grid_points+1) 
    x2=np.arange(0,len(u_star))
    """
    Initialize the pressure 
    """
    p_star[len(p_star)-1]=p_exit
    for i in range(len(p_star)-2,-1,-1):
        p_star[i]=p_star[i+1]+((dp*dx)/length)
    
    u_star=u_n
    Area=np.full((2 * Grid_points + 1), A)

    for t in range(0,n):
        start_time=time.time()
        i=0
        while True:
            print(i)                                                                          
            u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)  
            converge=convergence(u_star,Grid_points,rho,Area,A_n,dx,dt)
            p_add=Pressure_adjust(u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
            p_star=add(p_add,p_star,Grid_points)
            elasp_time = time.time() - start_time
            if elasp_time >5:
                break
        u_n=u_star
    plt.show()    
    return u_star,p_star

"""
Final result 
"""

velocity, Pressure =unsteady_1D_flow(A_n,A,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit)
