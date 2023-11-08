import numpy as np
import math as mt 


#pressure correction coefficients:
def get_pressure_coeffcients(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis):
    alpha=np.zeros(Grid_points-1)
    beta=np.zeros(Grid_points-1)
    gama=np.zeros(Grid_points-1)
    z=np.zeros(Grid_points-1)
    if (Grid_points>1): 
        for i in range(1,(2*Grid_points)-2,2):
            Re_left=(rho*u_n[(i-1)//2]*(mt.sqrt(4*A_n[i-1]/(mt.pi))))/d_vis
            Re_right=(rho*u_n[(i+1)//2]*(mt.sqrt(4*A_n[i+1]/(mt.pi))))/d_vis

            if Re_left < 2000:
                f_left=64/Re_left
            else :
                f_left=0.3164*Re_left**-0.25

            if Re_right < 2000:
                f_right=64/Re_right
            else :
                f_right=0.3164*Re_right**-0.25 

            a_left=(rho*Area[i-1]/dx)
            b_left=(rho*Area[i-1]/dt)
            c_left=(f_left*rho*u_n[(i-1)//2]*p_s[i-1]/8) //2
            a_right=(rho*Area[i+1]/dx)
            b_right=(rho*Area[i+1]/dt)
            c_right=(f_right*rho*u_n[(i+1)//2]*p_s[i+1]/8)

            beta1_left=(a_left/(b_left+c_left))
            gama1_right=(a_right/(b_right+c_right))

            alpha[(i-1)//2]=(beta1_left+gama1_right)*(Area[i]/dx)
            beta[(i-1)//2]=beta1_left*(Area[i]/dx)
            gama[(i-1)//2]=gama1_right*(Area[i]/dx)
            z[(i-1)//2]=((a_left*(u_star[(i-1)//2]))-(a_right*(u_star[(i+1)//2]))+(((rho*A_n[i])-(rho*Area[i]))/dt))

            
    return alpha,gama,beta,z