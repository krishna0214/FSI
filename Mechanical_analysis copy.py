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
from Boundary_conditions.velocity_p_time import velocity,Pressure
from Modules.friction_factor import friction_factors_shear
import pandas as pd
from Modules.mechanical_model import area
from Interpolate.inter import interpolate_Area

"""
const paramers through-out the analysis
All units are taken in SI system (Kilogram,Meter,Second)

"""
Total_time=1.0
#dt = 0.001                        # Time step size                 # Number of time steps
rho=1000                         # Density kg/m3
g = 9.8                                # Gravitational acceleration m/s2
Grid_points=250          #Total number of grid points (Cell centres including extra cell at the end)
d_vis=0.0035
dt=0.0001
n=int(Total_time/dt)
length=0.126                             #m
dx = length/Grid_points 




"""
Initialize the velocity 
"""


u_n = np.zeros(Grid_points+1) 
A_n=np.zeros(2*Grid_points+1)


def unsteady_1D_flow(A_n,u_n,Grid_points,rho,dx,dt,d_vis,n):
    graph=Graph_PV()
    Time=np.linspace(0,1,n)
    #len of p_star is u_star-1
    u_star = np.zeros(Grid_points+1) 
    p_star=np.ones(Grid_points)
    Area=np.zeros(2*Grid_points+1)  
    p_s=np.zeros(2*Grid_points+1) 
    A_n=np.ones(2*Grid_points+1)*(mt.pi*(0.003)**2)
    #A_n=np.ones(2*Grid_points+1)*0.22*(10**-4)
    u_n[0]=velocity(A_n[0],Time[0])
    
    for i in range(1,len(u_n)):
        u_n[i]=0.0001  
    Area=A_n
    u_star=u_n
    start_time=time.time()
     
    for t in range(0,5): 
        print(t,'t')
        print(Time[t],"time")
        p_s=2*(np.sqrt(np.pi*A_n))
        u_inlet=velocity(Area[0],Time[t])
        if u_inlet<=0:
            u_inlet=0.000001
        u_n[0]=u_inlet
        #print(u_inlet,"u_inlet")
        p_exit=Pressure(Time[t])*(10**3)
        p_star=p_star*p_exit
        #print(u_n[0])
        u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit) 
        i=0
        while True:
            i=i+1
            print(i)
            p_star=Correct_pressure(p_exit,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s)
            Area=area(p_star/100)
            Area=interpolate_Area(Area)
            converge=convergence(u_star,Grid_points,rho,Area,A_n,dx,dt)
            print(converge,"converge")
            if converge<10**-3:
                break 
            elif i>20:
                break                                                                        
            p_add=Pressure_adjust(u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
            relaxation_factor=1
            p_add=relaxation_factor*p_add
            #p_star=add(p_add,p_star,relaxation_factor)
            u_star=Correct_velocity(u_star,p_add,Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)
            #p_star=Correct_pressure(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s)
            #u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
        
        
        #print(p_star[0])
        x1=np.arange(len(u_star))
        x2=np.arange(len(p_star))
        x3=np.arange(len(Area))
        plot(u_n,p_star,Area,x1,x2,x3,graph)
        plt.pause(1) 
        
        
        u_n=u_star
        A_n=Area
        #print(A_n[0])  
    plt.show()    
    return u_star,p_star

"""
Final result 
"""
velocity,pressure =unsteady_1D_flow(A_n,u_n,Grid_points,rho,dx,dt,d_vis,n) 