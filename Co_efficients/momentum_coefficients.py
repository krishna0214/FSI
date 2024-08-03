import numpy as np
import math as mt


def get_abcd1(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s):
    b=np.ones((Grid_points+1))
    a=np.ones((Grid_points+1))
    c=np.ones((Grid_points))
    d=np.ones((Grid_points+1))
    
    for i in range (0,Grid_points+1):
        if u_n[i]==0:
            u_n[i]=0.00001
            Re=np.abs((rho*u_n[i]*(mt.sqrt(4*A_n[2*i]/(mt.pi))))/d_vis)
            if Re < 2000:
                f=64/Re
            else :
                f=0.3164*Re**-0.25 
        else:
            Re=np.abs((rho*u_n[i]*(mt.sqrt(4*A_n[2*i]/(mt.pi))))/d_vis)
            if Re < 2000:
                f=64/Re
            else :
                f=0.3164*Re**-0.25

        a[i]=(rho*Area[2*i]/dt)+(f*rho*u_n[i]*p_s[2*i]/8)
        b[i]=((rho*Area[2*i]/(2*dx)))
        d[i]=((rho*A_n[2*i]*u_n[i]/dt))
        if i<Grid_points:
            c[i]=(Area[(2*i)+1]/dx)#b=0.01*b
    #a=2*a

    return a,b,c,d 


def get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s):
    b=np.ones((len(u_n)))
    a=np.ones(len(u_n))
    c=np.ones((len(u_n)))
    d=np.ones(len(u_n))
    f=np.ones(len(u_n))
    
    for i in range (0,len(u_n)):
        if u_n[i]==0:
            u_n[i]=0.000001
            Re=np.abs((rho*u_n[i]*(mt.sqrt(4*A_n[2*i]/(mt.pi))))/d_vis)
            if Re < 2000:
                f[i]=64/Re
            else :
                f[i]=0.3164*Re**-0.25 
        else:
            Re=np.abs((rho*u_n[i]*(mt.sqrt(4*A_n[2*i]/(mt.pi))))/d_vis)
            if Re < 2000:
                f[i]=64/Re
            else :
                f[i]=0.3164*Re**-0.25

      
        a[i]=(rho*Area[2*i]/dt)+(f[i]*rho*u_n[i]*p_s[2*i]/8)
        #a[i]=(f[i]*rho*u_n[i]*p_s[2*i]/8)
        b[i]=((rho*Area[2*i]*u_n[i])/(2*dx))
        d[i]=((rho*A_n[2*i]*u_n[i])/(dt))
        #d[i]=0
        if i<len(u_n)-1:
            c[i]=(Area[(2*i)+1]/dx)#b=0.01*b
    #a=2*a

    return a,b,c,d 



def get_abcd_steady(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s):
    b=np.ones((len(u_n)))
    a=np.ones(len(u_n))
    c=np.ones((len(u_n)))
    d=np.ones(len(u_n))
    f=np.ones(len(u_n))
    
    for i in range (0,len(u_n)):
        if u_n[i]==0:
            u_n[i]=0.000001
            Re=np.abs((rho*u_n[i]*(mt.sqrt(4*A_n[2*i]/(mt.pi))))/d_vis)
            if Re < 2000:
                f[i]=64/Re
            else :
                f[i]=0.3164*Re**-0.25 
        else:
            Re=np.abs((rho*u_n[i]*(mt.sqrt(4*A_n[2*i]/(mt.pi))))/d_vis)
            if Re < 2000:
                f[i]=64/Re
            else :
                f[i]=0.3164*Re**-0.25

      
        a[i]=(f[i]*rho*u_n[i]*p_s[2*i]/8)
        #a[i]=(f[i]*rho*u_n[i]*p_s[2*i]/8)
        b[i]=((rho*Area[2*i]*u_n[i])/(2*dx))
        #d[i]=0
        if i<len(u_n)-1:
            c[i]=((Area[(2*i)]+Area[2*i+1]+Area[2*i+2])/3*dx)#b=0.01*b
    #a=2*a

    return a,b,c