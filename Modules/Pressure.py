import numpy as np
import math as mt 
import matplotlib.pyplot as plt
from Matrix_solver.LU_decompose import solve_lu,lu_decomposition
from Co_efficients.mass_coefficients import get_lamda_mu
from Co_efficients.momentum_coefficients import get_abcd
from Co_efficients.pressure_new import get_pressure
from Boundary_conditions.Inlet_Boundary import Inlet_Boundary_equation
from Plots.Plot import Graph_PV,plot

def correct_pressure(Grid_points,alpha,beta,gama,z,p_exit,c_1,c_2,z1):
    """
    Calculate P' vector using the provided inputs and matrix inversion.

    Parameters:
    alpha (numpy.ndarray):
    beta (numpy.ndarray): 
    gama (numpy.ndarray):
    z(numpy.ndarray):

    Returns:
    p' (numpy.ndarray): Calculated p'(It is the value of p that has to be added) vector.
    """
    
    n=Grid_points
    A1 = np.zeros((n-1, n))
    Rhs1=np.zeros((n-1,1))

    A1[0,0] =c_1
    A1[0,1] =-1*c_2
    Rhs1[0,0]=z1
    
    for i in range(1,n-1):
        A1[i,i-1] = gama[i-1]
        A1[i,i] = (-1*alpha[i])
        A1[i,i+1]= (beta[i+1])
        const = (-1*z[i])
        Rhs1[i,0]=const
        

    #print(A1,Rhs1)    
    Rhs1[n-2,0]=Rhs1[n-2,0]-(A1[n-2,n-1]*p_exit)
    co_efficientmatrix1=A1[:,0:n-1] 
    #print(co_efficientmatrix1,Rhs1)
    p_=solve_lu(co_efficientmatrix1,Rhs1)    
    return p_


def Pressure_adjust(u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit):

    """    
    Pressure Adjustment part 
    """    
    alpha,beta,gama,z=get_pressure(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis) 
    #print(u_star,"afteralphabeta")
    #print(alpha,beta,gama,z,"adjust1")
   
    """
    The equation is c_1p1'- c_2p2'= z1
    At inlet_cells using u_inlet

    """
    c_1,c_2,z1=Inlet_Boundary_equation(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
    add_p1=correct_pressure(Grid_points,alpha,beta,gama,z,p_exit,c_1,c_2,z1)#Solve pressure correction equation 
    #print(add_p1)
    return add_p1