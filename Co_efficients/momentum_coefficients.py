import numpy as np
import math as mt


def get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s):
    b=np.ones((Grid_points+1))
    a=np.ones((Grid_points+1))
    c=np.ones((Grid_points))
    d=np.ones((Grid_points+1))
    
    for i in range (0,Grid_points+1):
        Re=(rho*u_n[i]*(mt.sqrt(4*A_n[i]/(mt.pi))))/d_vis
        if Re < 2000:
            f=64/Re
        else :
            f=0.3164*Re**-0.25 
    
        a[i]=(rho*Area[2*i]/dt)+(f*rho*u_n[i]*p_s[2*i]/8)
        b[i]=((rho*Area[2*i]/(2*dx)))
        d[i]=((rho*A_n[2*i]*u_n[i]/dt))
        if i<Grid_points:
            c[i]=(Area[(2*i)+1]/dx)

    return a,b,c,d 


