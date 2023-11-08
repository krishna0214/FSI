import numpy as np
def get_lamda_m_u(Grid_points,rho,Area,A_n,dx,dt):
    lamda=np.zeros((Grid_points+1))
    m_u=np.zeros((Grid_points))
    for i in range (0,len(lamda)):
        lamda[i]=(rho*Area[2*i])/dx
    for j in range (0,len(m_u)):
        m_u[j]=(rho*(Area[2*j+1]-A_n[2*j+1]))/dt    
    
    return lamda, m_u