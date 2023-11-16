
import numpy as np

def pressure_Combine(p_star,p_star1,Grid_points):
    p_star1=p_star1.reshape(1,Grid_points-1)
    p_star[0:Grid_points-1]=p_star1[:, :]
    return p_star