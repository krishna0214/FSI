import numpy as np
import math as mt 
from Matrix_solver.LU_decompose import solve_lu,lu_decomposition
from Co_efficients.mass_coefficients import get_lamda_mu
from Co_efficients.momentum_coefficients import get_abcd
from Co_efficients.pressure_coefficients import get_pressure_coefficients
from Co_efficients.pressure_new import get_pressure

def correct_pressure(Grid_points,alpha,beta,gama,z,p_exit,c_1,c_2,z1,beta_boundary):
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
    #populate the coefficent matrix
    A1[0,0] =c_1
    A1[0,1] = -1*c_2
    const = (z1)
    Rhs1[0,0]=const

    for i in range(1,n-1):
        if i<n-1:
            A1[i,i-1] = gama[i-1]
            A1[i,i] = (-1*alpha[i-1])
            A1[i,i+1]= (beta[i-1])
            const = (-1*z[i-1])
            Rhs1[i,0]=const
        if i==n-1:
            A1[i,i-1] = gama[i-1]
            A1[i,i] = (-1*alpha[i-1])
            A1[i,i+1]= beta_boundary
            const = (-1*z[i-1])
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
    alpha,beta,gama,z=get_pressure_coefficients(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis) 
    #print(alpha,beta,gama,z)
    #print(alpha.size)

    """
    Inlet boundary equation co-efficinets 
    """

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


    """
    The equation is c_1p1'- c_2p2'= z1
    """
    Re_boundary=(rho*u_n[Grid_points-1]*(mt.sqrt(4*A_n[2*Grid_points-2]/(mt.pi))))/d_vis
    if Re_boundary < 2000:
        f_boundary=64/Re_boundary
    else :
        f_boundary=0.3164*Re_inlet**-0.25 

    x_1B_exit=(rho*Area[2*Grid_points-2]/dx)
    x_2B_exit=(rho*Area[2*Grid_points-2]/dt)
    x_3B_exit=(f_boundary*rho*u_n[Grid_points-1]*p_s[2*Grid_points-2]/8)
    
    beta_boundary=(x_1B_exit/(x_2B_exit+x_3B_exit))*(Area[2*Grid_points-1]/dx)
    
    add_p1=correct_pressure(Grid_points,alpha,beta,gama,z,p_exit,c_1,c_2,z1,beta_boundary)#Solve pressure correction equation 
    #print(add_p1)
    return add_p1