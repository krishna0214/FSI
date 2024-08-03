import numpy as np
from Co_efficients.momentum_coefficients import get_abcd
from Boundary_conditions.Outlet_Boundary import Outlet_Boundary_const


def Correct_velocity(u_star,p_add,Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s):
    a,b,c,d=get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)
    l1_extra_cell,mu_extra_cell,l2_extra_cell=Outlet_Boundary_const(rho,A_n,Area,dt,dx,Grid_points)
    u_update=np.zeros((len(u_star)))
    u_update[0]=u_star[0]
    for i in range(1,len(u_star)):
        if i<len(u_star)-1:
            u_update[i]=u_star[i]+ ((c[i-1]/a[i])*p_add[i-1])-((c[i]/a[i])*p_add[i])
        else:
            u_update[i]=(l1_extra_cell*u_star[i-1]-mu_extra_cell)/l2_extra_cell
    return u_update