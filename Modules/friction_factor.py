
import numpy as np
import math as mt

def frictionfactors_shear(rho, d_vis, A_n, u_n):
    friction_factor=np.zeros((u_n.shape[0],u_n.shape[1]))
    A_reshaped = A_n[::2]
    A_reshaped=A_reshaped*(np.ones_like(u_n))
    Re = np.abs((rho * u_n * np.sqrt(4 * A_reshaped / (mt.pi))) / d_vis)
    friction_factor = np.where(Re < 2000, 64 / Re, 0.3164 * Re ** -0.25)
    shear = (1/2) * rho * (u_n**2) * friction_factor        
    return friction_factor,shear 

