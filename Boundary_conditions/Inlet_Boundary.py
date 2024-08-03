
import math as mt
import numpy as np

"""
    Inlet boundary equation co-efficinets 
    """
def Inlet_Boundary_equation(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit):

    Re_inlet=(rho*u_n[1]*(mt.sqrt(4*A_n[2]/(mt.pi))))/d_vis
    if Re_inlet < 2000:
        f_inlet=64/Re_inlet
    else :
        f_inlet=0.3164*Re_inlet**-0.25 
    
    
    lamda_1= (rho*Area[0])/dx
    lamda_2=(rho*Area[2])/dx
    mu_in=(rho*Area[1]-rho*A_n[0])/dt
    a_2=((rho*Area[1]/dt)+(f_inlet*rho*u_n[1]*p_s[2]/8))
    c_1=(Area[1]/dx)
    c_2=(Area[3]/dx)
    z1=(a_2/lamda_2)*(lamda_1*u_inlet-lamda_2*u_star[1]-mu_in)
    return c_1,c_2,z1