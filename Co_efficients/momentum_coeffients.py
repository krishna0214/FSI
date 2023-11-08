import numpy as np
import math as mt


def get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s):
    b=np.ones((Grid_points))
    a=np.ones((Grid_points))
    c=np.ones((Grid_points-1))
    d=np.ones((Grid_points))
    
    for i in range (0,Grid_points):
        Re=(rho*u_n[i]*(mt.sqrt(4*A_n[i]/(mt.pi))))/d_vis
        if Re < 2000:
            f=64/Re
        else :
            f=0.3164*Re**-0.25 
        a[i]=((rho*Area[2*i]/dt)+(f*rho*u_n[i]*p_s[2*i]/8))
        b[i]=((rho*Area[2*i]/(2*dx)))
        d[i]=((rho*A_n[2*i]*u_n[i]/dt))  
        if i<Grid_points-1:
            c[i]=(Area[(2*i)+1]/dx)
        if i==Grid_points-1:
            mu=((rho*Area[2*i+1])-(rho*A_n[2*i+1]))/dt



    return a,b,c,d,mu  


