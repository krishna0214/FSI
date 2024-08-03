        


def Outlet_Boundary_const(rho,A_n,Area,dt,dx,Grid_points):
    i=Grid_points-1
    l1=rho*Area[2*i]/dx
    mu_extra_cell=((rho*Area[2*i+1])-(rho*A_n[2*i+1]))/dt
    l2=rho*Area[2*(i+1)]/dx
    return l1,mu_extra_cell,l2
