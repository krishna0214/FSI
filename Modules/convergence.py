import numpy as np
from Co_efficients.mass_coeffcinets import get_lamda_m_u

def convergence(u_star,Grid_points,rho,Area,A_n,dx,dt):
    lamda,m_u=get_lamda_m_u(Grid_points,rho,Area,A_n,dx,dt)
    residue=np.zeros((Grid_points-1))
    for i in range(0,len(residue)):
        residue[i]=(lamda[i+1]*u_star[i+1])-(lamda[i]*u_star[i])+(m_u[i])
    mean= np.mean(np.abs(residue))
    #squared_differences = (residue - mean_with_mod) ** 2
    #mean_squared_error = squared_differences.mean()
    #rmse = np.sqrt(mean_squared_error)
    return mean