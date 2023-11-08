import numpy as np
from Co_efficients.momentum_coeffients import get_abcd
from Matrix_solver.LU_decompose import solve_lu
from Update.velocity_update import Combine


def solve_momentum(a, b, c, d, p_star,Grid_points,mu,u_entry,p_exit):
    """
    Calculate u_star vector using the provided inputs and matrix inversion.

    Parameters:
    a (numpy.ndarray): Coefficients a_{i+1/2}^{n+1}.
    b (numpy.ndarray): Coefficients b_{i-1/2}^{n+1} & b_{i+3/2}^{n+1}.
    c (numpy.ndarray): Coefficients c_i^{n+1}.
    d (numpy.ndarray): Coefficients d_{i+1/2}^{n}.
    p_star (numpy.ndarray): Guessed pressure vector.

    Returns:
    u_star (numpy.ndarray): Calculated u_star vector.
    """
    n =Grid_points
    if n>2:
                                                       # Initialize the coefficient matrix A and the right-hand side vector Rhs
        A = np.zeros((n-1, n))
        Rhs=np.zeros((n-1,1))
                                                       # Populate the coefficient matrix and the right-hand side vector
        for i in range(1,n):
            if i<n-1:
                A[i-1, i-1] =b[i-1]
                A[i-1,i] = (-1*a[i])
                A[i-1,i+1]=b[i+1]
                const = ((c[i] * p_star[i]) - (c[i-1] * p_star[i-1]) - d[i])
                Rhs[i-1,0]=const
            else:
                A[i-1, i-1] =b[i-1]
                A[i-1,i] = (-1*(a[i]-b[i]))

                const = ((c[i-1] * p_exit) - (c[i-1] * p_star[i-1]) - d[i]+ mu/2)
                Rhs[i-1,0]=const

        #print(A,Rhs) 
        Rhs[0,0]= Rhs[0,0]-(A[0,0]*u_entry)
        #print(Rhs[0,0],Rhs[n-3,0])     
        coefficient_matrix= A[:,1:n]
        #print(coefficient_matrix,Rhs)       
        u_star1=solve_lu(coefficient_matrix,Rhs)
        #print(solve_linear_system(coefficient_matrix,Rhs))

    return u_star1


def Momentum(p_star,u_n,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit):
    """    
    #Mometum equation part 
    """
    a_values,b_values,c_values,d_values,mu=get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)       #solve co-efficients at all grid points
    #print(u_n)
                                                                                                        #define Boubdary conditions for velocity
    u_star1= solve_momentum(a_values, b_values, c_values, d_values, p_star,Grid_points,mu,u_inlet,p_exit)# Solve momentum equation to find intermediate velocity (u_star)
    #print(u_n)
    #print(u_star,u_star1)
    u_temp= u_n.copy()
    u_star=Combine(u_star1,u_temp,Grid_points)
    #print(u_n,'at combine')
    #print(u_star)
    return u_star