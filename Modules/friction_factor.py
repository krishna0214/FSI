
import numpy as np
import math as mt

"""
def friction_factors_shear(rho, d_vis, A_n, u_n):
    friction_factor=np.zeros(len(u_n))
    A_reshaped = A_n[::2]
    A_reshaped=A_reshaped*(np.ones_like(u_n))
    Re = np.abs((rho * u_n * np.sqrt(4 * A_reshaped / (mt.pi))) / d_vis)
    friction_factor = np.where(Re < 2000, 64 / Re, 0.3164 * Re ** -0.25)
    shear = (1/2) * rho * (u_n**2) * friction_factor        
    return friction_factor,shear 


    def friction_factors_shear(rho, d_vis, A_n, u_n):
    shear=np.ones_like(u_n)
    for i in range (0,len(u_n)):
        Re=(rho*u_n[i]*(mt.sqrt(4*A_n[2*i]/(mt.pi))))/d_vis
        if Re < 2000:
            f=64/Re
        else :
            f=0.3164*Re**-0.25 
        shear[i]=(f*rho*(u_n[i]**2))/(8)   
    return shear 

     friction_factor=np.zeros((u_n.shape[0],u_n.shape[1]))
   
"""
def friction_factors_shear1(rho, d_vis, A_n, u_n):
    friction_factor=np.zeros((u_n.shape[0],u_n.shape[1]))
    A_reshaped = A_n[::2]
    A_reshaped=A_reshaped*(np.ones_like(u_n))
    u_n = np.where(u_n == 0, 0.001, u_n)
    Re = np.abs((rho * u_n * np.sqrt(4 * A_reshaped / (mt.pi))) / d_vis)
    friction_factor = np.where(Re < 2000, 64 / Re, 0.3164 * Re ** -0.25)
    shear = (1/8) * rho * (u_n**2) * friction_factor        
    return friction_factor,shear

def friction_factors_shear(rho,d_vis,A_n,u_n,dx):
    A_reshaped = A_n[::2]
    #u_interpolated=np.array(0,len(u_n)-1)
    #friction_factor_inter=np.array(0,len(u_n)-1)
    #A_interpolated=np.array(0,len(u_n)-1)
    #A_difference=np.array(0,len(u_n)-1)
    #shear=np.array(0,len(u_n)-1)
    #u_n = np.where(u_n == 0, 0.001, u_n)
    D=np.sqrt(4 * A_reshaped / (mt.pi))
    #D_inter=np.array(0,len(u_n)-1)
    #D1=np.sqrt(4 * A_n / (mt.pi))
    Re = np.abs((rho * u_n * D) / d_vis)
    for i in range(0,len(Re)):
        if Re[i]<2000:
            friction_factor=64/Re
        elif Re[i]>2000:
            friction_factor= 0.3164 * Re ** -0.25   
    #friction_factor = np.where(Re < 2000, 64 / Re, 0.3164 * Re ** -0.25)
    #for i in range(0,len(u_n)-1):
        #u_interpolated[i]=(u_n[i]+u_n[i+1])/2
        #D_inter=(D[i]+D[i+1])/2
        #friction_factor_inter=(friction_factor[i]+friction_factor[i+1])/2
        #A_interpolated=(A_reshaped[i]+A_reshaped[i+1])/2
        #A_difference[i]=(A_reshaped[i]-A_reshaped[i+1])/2
        #N=(A_difference[i])/(mt.pi*D_inter[i]*dx)
        #shear[i] = (1/8) * rho * (u_interpolated[i]**2) * friction_factor_inter[i]*N
        
    N=(A_reshaped)/(mt.pi*D*dx)
    shear = (1/2) * rho * (u_n**2) * friction_factor*N
    #shear=8*(d_vis)*(A_n[0]*u_n[0])/(A_n*D1**2)    
    return friction_factor,shear,Re,A_reshaped,D