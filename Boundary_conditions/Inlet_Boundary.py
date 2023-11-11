
import math as mt
import numpy as np

"""
    Inlet boundary equation co-efficinets 
    """
def Inlet_Boundary_equation(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit):

    Re_inlet=(rho*u_n[2]*(mt.sqrt(4*A_n[2]/(mt.pi))))/d_vis
    if Re_inlet < 2000:
        f_inlet=64/Re_inlet
    else :
        f_inlet=0.3164*Re_inlet**-0.25 

    mu1=((rho*Area[1])-(rho*A_n[1]))/dt
    mu2=((rho*Area[3])-(rho*A_n[3]))/dt
    u2= (((rho*Area[0])/dx)*u_inlet-(mu1))/((rho*Area[2])/dx)
    u3=(((rho*Area[0])/dx)*u_inlet-(mu1+mu2))/((rho*Area[4])/dx)
    a1=((rho*Area[2]/dt)+(f_inlet*rho*u_n[2]*p_s[2]/8))
    b0=(rho*Area[0]/(2*dx))
    b2=(rho*Area[4]/(2*dx))
    c_1=(Area[1]/dx)
    c_2=(Area[3]/dx)
    z1=a1*(u2-u_star[1])-b0*(u_inlet-u_star[0])-b2*(u3-u_star[2])

    return c_1,c_2,z1