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
from Area.linear_area import Linear_area_profile, Const_area_profile,Area,stenoisis_Area5
from Boundary_conditions.velocity_p_time import velocity
from Modules.friction_factor import friction_factors_shear
import pandas as pd
"""
const paramers through-out the analysis
All units are taken in SI system (Kilogram,Meter,Second)

"""
Total_time=1.0


#dt = 0.001                        # Time step size                 # Number of time steps
#rho = 1000
rho=1060                         # Density kg/m3
g = 9.8                                # Gravitational acceleration m/s2
Grid_points=200            #Total number of grid points (Cell centres including extra cell at the end)
#Inletmassflux=900                      #kg/m2.s
P_atm=1 
#d_vis=1 
d_vis=0.0035
u_inlet=0.2                      #Ns/m2                              
#dia= (RE*d_vis)/(rho*u_inlet)
dia=0.0039

print(dia,"dia")    #m
p_exit= P_atm* 101325                  #N/m2
#A= (mt.pi)*(dia**2)/4                  #m2
#u_inlet= Inletmassflux/(rho) 
length=0.07                              #m
dx = length/Grid_points 
minimum_dt=dx/u_inlet     #<0.0044
dt=minimum_dt                 # Spatial grid size
n = int(Total_time/dt)  
#n=1
CFL=u_inlet*dt/dx


"""
Variables evolve-with time
"""
a1=(2*length)/(12.5*dia)

r_x,A_n,p_s,x5=stenoisis_Area5(dia,Grid_points,length)

"""
Initialize the velocity 
"""
print(A_n[0],"A_n[0]")
u_n = np.zeros(Grid_points+1) 
u_n[0]=u_inlet
for i in range(1,len(u_n)):
    u_n[i]=0.0001    

# Main Algo
def unsteady_1D_flow(A_n,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit):
    p_star = np.zeros(Grid_points)
    u_star = np.zeros(Grid_points+1) 
    x1=np.arange(len(p_star))
    x2=np.arange(len(u_star))

    """
    Initialize the pressure 
    """
    p_star[len(p_star)-1]=p_exit
    for i in range(len(p_star)-2,-1,-1):
        #p_star[i]=p_star[i+1]+((dp*dx)/length)
        p_star[i]=p_star[i+1]
    
    u_star=u_n
    Area=A_n
    graph_v=Graph_V()
    Time=np.linspace(0,0.8,n) 
    #u_inlet=u_inlet*(mt.cos(mt.pi*Time[0]))
    #graph_v=Graph_V()
    u_in=u_inlet
    for t in range(0,1):
        print(t,"t")
        #graph=Graph_PV()
        #graph_v=Graph_V()
        #graph_p=Graph_P()
        print(Time[t],"time")
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
            if converge<10**-4:
                break                                                                         
            p_add=Pressure_adjust(u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
            relaxation_factor=0.8
            p_add=relaxation_factor*p_add
            #p_star=add(p_add,p_star,relaxation_factor)
            u_star=Correct_velocity(u_star,p_add,Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)
            p_star=Correct_pressure(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s)
            #u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
        
        u_n=u_star
       
    #y=np.linspace(0,length,100)
    xx=np.linspace(0,length,Grid_points+1)
    friction_factor,shear,Re,A_re,D=friction_factors_shear(rho,d_vis,A_n,u_star,dx)
    xx1=np.linspace(0,length,len(shear))
    #print(shear,"shear")
    #print(xx,"xx")
    #u_star1=np.zeros(len(y))
    #for i in range(0,len(y)):
        #u_star1[i]=u_star[int(y[i]/dx)]
    # Load DataFrame from the pickle file
    #df= pd.read_pickle(r'E:\BTP\Linux\Fluid_flow\Unsteady-Fluid-Flow-Model\Results\results.pkl')
    #print(df.index)
    #df.loc[:, ('Jawahar', 'Stenosis', '70%', 'velocity')] = u_star
    #df.loc[:, ('Jawahar', 'Stenosis', '70%', 'dist')] = xx
    #df.to_pickle(r'E:\BTP\Linux\Fluid_flow\Unsteady-Fluid-Flow-Model\Results\results.pkl')
    plot_V((u_n/D)*(A_re*1),xx1,graph_v)
    #print(df,"df")
    plt.pause(2)    
    plt.show()    
    return u_star,p_star

"""
Final result 
"""
velocity, Pressure =unsteady_1D_flow(A_n,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit)





