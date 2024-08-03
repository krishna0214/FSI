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
#from Plots.Pressure_plot import plot_pressure
from Plots.Plot import Graph_PV,plot
from Plots.Velocity_plot import Graph_V,plot_V
from Plots.Velocity_w_Time import Graph_V_w_time,plot_V_time,plot_P_time
from Plots.Pressure_plot import  Graph_P,plot_P 
from Plots.Area_profile import Graph_Area_profile,plot_Area
from Plots.Time_averaged_shear import Graph_shear,plot_shear
from Modules.Correct_velocity import Correct_velocity
from Modules.Correct_pressure import Correct_pressure
from Area.linear_area import stenoisis_Area5
from Boundary_conditions.velocity_p_time import velocity
from Modules.friction_factor import friction_factors_shear
import pandas as pd
"""
const paramers through-out the analysis
All units are taken in SI system (Kilogram,Meter,Second)

"""
Total_time=1.0
#dt = 0.001                        # Time step size                 # Number of time steps
rho=1060                         # Density kg/m3
g = 9.8                                # Gravitational acceleration m/s2
Grid_points=400          #Total number of grid points (Cell centres including extra cell at the end)
#x=110
P_atm=1
d_vis=0.004
#u_inlet=0.2                    #Ns/m2                              
dia=0.0039
p_exit= P_atm* 101325                 #N/m2
#A= (mt.pi)*(dia**2)/4                  #m2
#u_inlet= Inletmassflux/(rho) 
length=0.07                              #m
dx = length/Grid_points 
RE=500
u_inlet=0.2
minimum_dt=dx/u_inlet     #<0.0044
dt=minimum_dt                 # Spatial grid size
n = int(Total_time/dt)  
#n=1
CFL=u_inlet*dt/dx

print(dia,"dia")
"""
Variables evolve-with time
"""


r_x,A_n,p_s,x5=stenoisis_Area5(dia,Grid_points,length)
#r_x,A_n,p_s=Const_area_profile(dia ,length,Grid_points)
#r_x,A_n,p_s=Linear_area_profile(dia ,length,Grid_points)

#graph_v=Graph_V()
#x2=np.arange(len(A_n))
#plot_V(A_n*10**4,x2*100*dx/2,graph_v)
#plt.pause(2) 
#plt.show()

"""
Initialize the velocity 
"""
print(A_n[0],"A_n[0]")
u_n = np.zeros(Grid_points+1) 
u_n[0]=u_inlet
for i in range(1,len(u_n)):
    u_n[i]=0.0001    



'''
Get r_x at some length
for i in range(0,29):
    L1=i/4
    n11=int(L1*2*Grid_points/(length*100))
    print(A_n[n11]*10**4,"A_n at " , L1)

'''
for i in range(0,29):
    L1=i/4
    n11=int(L1*2*Grid_points/(length*100))
    print(r_x[n11]*10**4,"r_n at " , L1)

data = """
0	0.2
3.50E-03	0.206348
7.00E-03	0.21291
1.05E-02	0.21966
1.40E-02	0.224755
1.75E-02	0.232139
2.10E-02	0.242742
2.45E-02	0.293371
2.80E-02	0.42439
3.15E-02	0.628085
3.50E-02	0.762172
3.85E-02	0.665265
4.20E-02	0.486604
4.55E-02	0.361606
4.90E-02	0.320664
5.25E-02	0.329156
5.60E-02	0.342336
5.95E-02	0.352277
6.30E-02	0.366753
6.65E-02	0.382526
7.00E-02	0.399339
"""

# Split the data into lines and then split each line into columns
lines = data.strip().split('\n')
values = [line.split() for line in lines]

# Convert values to numpy arrays
values = np.array(values)

# Convert values to float
values = values.astype(float)

# Separate the values into two numpy arrays
axial_values = values[:, 0]
velocity_values = values[:, 1]

# Print the arrays
print("Axial Values:", axial_values)
print("Velocity Values:", velocity_values)


#print(p_exit*760/101325,"p_exit")
#Area=np.array(A_n)
#print(A_n.shape,"A_n")
#print(Area.shape,"Area")
# Main Algo
def unsteady_1D_flow(A_n,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit):
    #len of p_star is u_star-1
    u_star = np.zeros(Grid_points+1)
    #p_star = np.zeros(Grid_points)
    #p_star[len(p_star)-1]=p_exit 
    x2=np.arange(len(u_star))
    x1=np.arange(len(u_star)-1)

    """
    Initialize the pressure 
    """
    u_star=u_n
    
    
    graph_v=Graph_V()
    #u_inlet=u_inlet*(mt.cos(mt.pi*Time[0]))
    #graph_v=Graph_V()
    for t in range(0,1):
        print(t,"t")
        Area=np.array(A_n)
        #graph=Graph_PV()
        #graph_v=Graph_V()
        #graph_p=Graph_P()
        #graph=Graph_PV()
        #graph_v=Graph_V()
        #graph_p=Graph_P()
        i=0
        #u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit) 
        while True:
            i=i+1
            print(i)
            converge=convergence(u_star,Grid_points,rho,Area,A_n,dx,dt)
            print(converge,"converge")
            if converge<10**-3:
                break                                                                         
            p_add=Pressure_adjust(u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
            relaxation_factor=0.8
            p_add=relaxation_factor*p_add
            #p_star=add(p_add,p_star,relaxation_factor)
            u_star=Correct_velocity(u_star,p_add,Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)
            #p_star=Correct_pressure(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s)
            #u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
            #plot_V(u_star,x2,graph_v)
        #print(df,"df")
            #plt.pause(1) 
        u_n=u_star
        #print(u_n,"u_n")
        p_star=Correct_pressure(p_exit,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s)
        #print(u_star,"p_star")
        #x=np.linspace(0,len(Area),len(Area))\
        data = np.column_stack((u_star, x2))

    # Save data to CSV file
        np.savetxt("velocity.csv", data, delimiter=",", header="u_star,x2", comments="")    
        x2=x2*dx
        plot_V(u_star,x2,axial_values,velocity_values,graph_v)
        #plot_V(u_star,x2,graph_v)
        #print(df,"df")
        plt.pause(2) 
    
    
    plt.show()    
    return u_star,p_star

"""
Final result 
"""
velocity, Pressure =unsteady_1D_flow(A_n,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit)





