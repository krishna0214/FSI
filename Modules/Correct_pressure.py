import numpy as np
from Co_efficients.momentum_coefficients import get_abcd
from Matrix_solver.LU_decompose import solve_lu
from Update.pressure_combine import pressure_Combine

def u_hat(Grid_points,u_star,b):
    u_nb=np.zeros(len(u_star))
    u_nb[0]=0
    u_nb[len(u_nb)-1]=0
    for i in range(1,len(u_star)-1):
        u_nb[i]=(u_star[i-1]*b[i-1]-u_star[i+1]*b[i+1])
    return u_nb

"""
def solve_pressure(a, c, d, Grid_points,u_hat,p_exit):
    n=Grid_points-1                                             
    A = np.zeros((n, n+1))
    Rhs=np.zeros((n,1))                           
    for i in range(1,n+1):
        A[i-1, i-1] = -1*c[i-1]
        A[i-1,i] = c[i]
        const = u_hat[i-1]+d[i]-a[i]
        Rhs[i-1,0]=const
    #print(A,Rhs) 
    Rhs[0,0]= Rhs[0,0]-(A[n-1,n-1]*p_exit)
    #print(Rhs[0,0],Rhs[n-3,0])     
    coefficient_matrix= A[:,0:n]
    #print(coefficient_matrix,Rhs)       
    p_star1=solve_lu(coefficient_matrix,Rhs)
    #print(solve_linear_system(coefficient_matrix,Rhs))   
    return p_star1



def Correct_pressure(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit):
   
    a_values,b_values,c_values,d_values=get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)       #solve co-efficients at all grid points
    #Imposed continuity on Extra cell at outlet boundary 
    u_hat1=u_hat(Grid_points,u_star,b_values)
    p_star1= solve_pressure(a_values, c_values, d_values, Grid_points, u_hat1,p_exit)# Solve momentum equation to find intermediate velocity (u_star)
    p_star=pressure_Combine(p_star,p_star1,Grid_points)
    #print(u_n,'at combine')
    #print(u_star)
    #651
    return p_star

"""

def Correct_pressure(p_exit,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s):
    a,b,c,d=get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)
    #u_hat1=u_hat(Grid_points,u_star,b)
    n=len(u_star)-1
    p_updated=np.ones(n)*1000
    p_updated[len(p_updated)-1]=p_exit
    for i in range(len(p_updated)-2,-1,-1):
        #p_updated[i]=(c[i+1]*p_updated[i+1])+(d[i+1])+(a[i+2]*u_star[i+2])-(a[i]*u_star[i])+((a[i+1]*u_star[i+1]))/c[i]
        #p_updated[i]=((c[i+1]*p_updated[i+1])-(a[i]*u_star[i]))
        c1=Area[i+1]/dx
        c2=Area[i]/dx
        p_updated[i]=((c1*(p_updated[i+1]))-(d[i+1])+(a[i+1]*u_star[i+1])+(a[i+2]*u_star[i+2])-(a[i]*u_star[i]))/c2
    #print(p_updated,"p_updated")    
    return p_updated
