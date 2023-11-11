import numpy as np
import math as mt 
from Co_efficients.momentum_coefficients import get_abcd
from Co_efficients.mass_coefficients import get_lamda_mu


#pressure correction coefficients:
def get_pressure(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis):

    alpha=np.zeros(Grid_points)
    beta=np.zeros(Grid_points)
    gama=np.zeros(Grid_points)
    z= np.zeros(Grid_points)

    a,b,c,d=get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)
    lamda,mu=get_lamda_mu(Grid_points,rho,Area,A_n,dx,dt)

    for i in range (0,Grid_points):
        alpha[i]=((lamda[i]/a[i])+(lamda[i+1]/a[i+1]))*(c[i])
        beta[i]=(c[i])*(lamda[i]/a[i])
        gama[i]=(c[i])*(lamda[i+1]/a[i+1])
        z[i]=mu[i]+(lamda[i+1]*u_star[i+1])+(lamda[i]*u_star[i])       
    return alpha,gama,beta,z