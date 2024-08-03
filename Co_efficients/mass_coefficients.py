import numpy as np
import math as mt

def get_lamda_mu(Grid_points,rho,Area,A_n,dx,dt):
    lamda=np.zeros((Grid_points+1))
    mu=np.ones((Grid_points))
    for i in range (0,len(lamda)):
        lamda[i]=(rho*Area[2*i])/dx
        
    for j in range (0,len(mu)):
        mu[j]=(rho*(Area[2*j+1]-A_n[2*j+1]))/dt    
    return lamda, mu