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
from Area.linear_area import Linear_area_profile, Const_area_profile,Area,stenoisis_Area
from Boundary_conditions.velocity_p_time import velocity
from Modules.friction_factor import friction_factors_shear
"""
const paramers through-out the analysis
All units are taken in SI system (Kilogram,Meter,Second)

"""
Total_time=1.0


#dt = 0.001                        # Time step size                 # Number of time steps
rho = 1060                         # Density kg/m3
g = 9.8                                # Gravitational acceleration m/s2
Grid_points=500              #Total number of grid points (Cell centres including extra cell at the end)
#Inletmassflux=900                      #kg/m2.s
P_atm=1                                #atm
dia=0.02                             #m
d_vis=0.004                         #Ns/m2
p_exit= P_atm* 101325                  #N/m2
#A= (mt.pi)*(dia**2)/4                  #m2
#u_inlet= Inletmassflux/(rho) 
u_inlet=0.042          #m/s
length=0.07                              #m
dx = length/Grid_points 
dt=0.0040             # Spatial grid size
n = int(Total_time/dt)  
#n=1
CFL=u_inlet*dt/dx
minimum_dt=dx/u_inlet     #<0.0044


"""
Variables evolve-with time
"""
a1=(2*length)/(12.5*dia)
l1=(length/2)-(length/a1)
l2=(length/2)+(length/a1)
#r_x,A_n,p_s=Linear_area_profile(dia ,length,Grid_points)  #Area and perimeter at n_th step  
#r_x,A_n,p_s=Area(dia,Grid_points,length,l1,l2)
r_x,A_n,p_s=Const_area_profile(dia,length,Grid_points)
#r_x,A_n,p_s=stenoisis_Area(dia,Grid_points,length,l1,l2)
print(A_n[0],"A_n[0]")
n1=int(Grid_points/2)
print(A_n[n1],"A_n[Grid_points/2]")
print(((A_n[0]-A_n[n1])/(A_n[0]))*100,"A_n[Grid_points/2]-A_n[0] it is the reduction ")

"""
Initialize the velocity 
"""
print(A_n[0],"A_n[0]")
n1=int(Grid_points/2)
print(A_n[n1],"A_n[Grid_points/2]")
print(A_n[n1]-A_n[0],"A_n[Grid_points/2]-A_n[0] it is the reduction ")
u_n = np.zeros(Grid_points+1) 
u_n[0]=u_inlet
for i in range(1,len(u_n)):
    u_n[i]=0.001    

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
    u_temp=np.zeros((n,len(u_star)))
    shear_temp=np.zeros((n,len(u_star)))
    pressure_temp=np.zeros((n,len(p_star)))
    #Area=np.full((2 * Grid_points + 1), A)
    Area=A_n
    #graph=Graph_PV()
    #graph_p=Graph_P()
    #graph_v=Graph_V()
    Time=np.linspace(0,0.8,n) 
    #u_inlet=u_inlet*(mt.cos(mt.pi*Time[0]))
    graph_v=Graph_V()
    u_in=u_inlet
    for t in range(0,n):
        print(t,"t")
        #graph=Graph_PV()
        #graph_v=Graph_V()
        #graph_p=Graph_P()
        print(t,"t")
        print(Time[t],"time")
        #graph=Graph_PV()
        #graph_v=Graph_V()
        #graph_p=Graph_P()
        #for i in range(1,len(u_n)):
        #u_n[i]=u_inlet
        start_time=time.time()
        u_inlet=velocity(A_n[0],Time[t])
        if u_inlet<=0:
            u_inlet=0.001
        u_n[0]=u_inlet
        print(u_inlet,"u_inlet")
        #plot_V(u_star,x2,graph_v)
        #plt.pause(0.1)
        #plot_P(x1,p_star,graph_p)
        #plt.pause(0.1)
        #plot(p_star,u_star,x1,x2,graph)
        #plt.pause(0.5)
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
            p_star=Correct_pressure(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s)
            #u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
            #plot(p_star,u_star,x1,x2,graph)
            #plt.pause(1) 
            #relaxation_factor=10**-5
            #p_star=add(p_add,p_star,relaxation_factor)
            #plot_P(x1,p_star,graph_p)
            #plt.pause(0.5)
            #u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
            #plot_P(x1,p_star,graph_p)
            #plt.pause(0.5)
            #u_star=u_update
            #plot_V(shear[1],x2,graph_v)
            #plt.pause(0.1)
            
            # Save the current figure with a unique file name
            # # Unique file name based on iteration
            #plt.savefig(file_name)
        
            #plt.clf()  # Clear the current figure to prepare for the next iteration
         
        #print(p_update[0],"p_update")
        #print(p_star[0],"p_star")
        #plot_P(x1,p_update,graph_p)
        #plt.pause(0.05)
        u_temp[t,:]=u_star
        pressure_temp[t,:]=p_star
        #shear_temp[t,:]=shear[1]
        #print(graph[1],"graph")
        #plot_V(u_star,x2,graph_v)
        #plt.pause(0.1)
        u_n=u_star
        #plot(p_star,u_star,x1,x2,graph)
        #plt.pause(0.1)
    #graph1=Graph_V_w_time()    
    #plot_V_time(u_temp,x2,Time,graph1,"velocity")
    #plt.pause(0.5)
    #graph2=Graph_V_w_time()
    #plot_P_time(pressure_temp,x1,Time,graph2,"pressure")
    #plt.pause(0.5) 
    shear=friction_factors_shear(rho, d_vis, A_n, u_temp)
    graph3=Graph_V_w_time()
    plot_V_time(shear[1],x2,Time,graph3,"shear")
    plt.pause(0.5)
    #np.savetxt('Velocity_dia.dat',u_star, delimiter=',')
    #np.savetxt('Pressure_dia.dat',p_star, delimiter=',')
    #np.savetxt('distanct_dia.dat',x2, delimiter=',')
    #np.savetxt('shear_dia.dat',shear[1], delimiter=',')
    #np.savetxt('u_temp_dia.dat',u_temp, delimiter=',')
    #np.savetxt('pressure_temp_dia.dat',pressure_temp, delimiter=',')

    time_avg_shear=((np.sum(shear[1],axis=0))/(Total_time))*(dt)
    time_avg_velocity=((np.sum(u_temp,axis=0))/(Total_time))*(dt)
    time_avg_pressure=((np.sum(pressure_temp,axis=0))/(Total_time))*dt
    graph_shear=Graph_shear()
    avg_shear=time_avg_shear.reshape(len(x2))
    print(avg_shear[0],"avg_shear")
    print(avg_shear.shape,"avg_shear.shape",x2.shape,"x2.shape")
    plot_shear(time_avg_shear,x2,graph_shear)
    #plot_V(time_avg_velocity,x2,graph_v)
    plt.pause(0.5)
    #np.savetxt('time_avg_velocity.dat',time_avg_velocity, delimiter=',')
    #np.savetxt('time_avg_pressure.dat',time_avg_pressure, delimiter=',')
    #np.savetxt('time_average_shear.dat',time_avg_shear, delimiter=',')
    plt.show()    
    return u_star,p_star

"""
Final result 
"""
velocity, Pressure =unsteady_1D_flow(A_n,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit)
