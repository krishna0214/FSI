
import numpy as np

def Combine(u_star,u_star1,Grid_points):
    u_star1=u_star1.reshape(1,Grid_points)
    u_star[1:Grid_points+1]=u_star1[:, :]
    return u_star
     
