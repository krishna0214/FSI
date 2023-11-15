import numpy as np
import os
import time 
import math as mt
import matplotlib.pyplot as plt


"""
const paramers through-out the analysis
All units are taken in SI system (Kilogram,Meter,Second)

"""
Total_time=1
dt=0.001                            # Time step size
#n = int(Total_time/dt) 
n=1                # Number of time steps
rho = 996.550                          # Density kg/m3
g = 9.8                                # Gravitational acceleration m/s2
Grid_points=100                    #Total number of grid points (Cell centres including extra cell at the end)
Inletmassflux=900                      #kg/m2.s
P_atm=1                                #atm
dia=0.0154                             #m
d_vis=0.000854                         #Ns/m2
p_exit= P_atm* 101325                  #N/m2
A= (mt.pi)*(dia**2)/4                  #m2
u_inlet= Inletmassflux/(rho)           #m/s
length=2                               #m
dx = length/Grid_points 
dt=0.001               # Spatial grid size
CFL=u_inlet*dt/dx
minimum_dt=dx/u_inlet     #<0.0044

"""
Variables evolve-with time
"""
# Initialize arrays for velocity (u) and pressure (p) at staggered grid locations
A_n = np.full((2 * Grid_points + 1), A)                # Area at n_th step
perimeter= (mt.pi)*dia                                 #perimeter at n_th step 
p_s=np.full((2 * Grid_points + 1), perimeter)          #perimeter at n_th step
dp= 128*d_vis*length*A*u_inlet/((mt.pi)*dia**4)        #pressure change for end points(Change the length to get the desired pressure drop)
print(dp,"dp")

"""
Initialize the velocity 
"""
u_n = np.zeros(Grid_points+1)  
u_n[0]=u_inlet
for i in range(1,len(u_n)):
    u_n[i]=0.01    

#Boundary_conditions
def Inlet_Boundary_equation(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit):
    Re_inlet=(rho*u_n[2]*(mt.sqrt(4*A_n[2]/(mt.pi))))/d_vis
    if Re_inlet < 2000:
        f_inlet=64/Re_inlet
    else :
        f_inlet=0.3164*Re_inlet**-0.25 
    mu1=((rho*Area[1])-(rho*A_n[1]))/dt
    mu2=((rho*Area[3])-(rho*A_n[3]))/dt
    u2= (((rho*Area[0])/dx)*u_inlet-(mu1))/((rho*Area[2])/dx)
    u3=(((rho*Area[0])/dx)*u_inlet-(mu1+mu2))/((rho*Area[4])/dx)
    a1=((rho*Area[2]/dt)+(f_inlet*rho*u_n[2]*p_s[2]/8))
    b0=(rho*Area[0]/(2*dx))
    b2=(rho*Area[4]/(2*dx))
    c_1=(Area[1]/dx)
    c_2=(Area[3]/dx)
    z1=a1*(u2-u_star[1])-b0*(u_inlet-u_star[0])-b2*(u3-u_star[2])
    return c_1,c_2,z1
def Outlet_Boundary_const(rho,A_n,Area,dt,dx,Grid_points):
    i=Grid_points-1
    l1=rho*Area[2*i]/dx
    mu_extra_cell=((rho*Area[2*i+1])-(rho*A_n[2*i+1]))/dt
    l2=rho*Area[2*(i+1)]/dx
    return l1,mu_extra_cell,l2    
#Co-effiecinet matrixs
def get_lamda_mu(Grid_points,rho,Area,A_n,dx,dt):
    lamda=np.zeros((Grid_points+1))
    mu=np.ones((Grid_points))
    for i in range (0,len(lamda)):
        lamda[i]=(rho*Area[2*i])/dx   
    for j in range (0,len(mu)):
        mu[j]=(rho*(Area[2*j+1]-A_n[2*j+1]))/dt    
    return lamda, mu
def get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s):
    b=np.ones((Grid_points+1))
    a=np.ones((Grid_points+1))
    c=np.ones((Grid_points))
    d=np.ones((Grid_points+1))
    for i in range (0,Grid_points+1):
        Re=(rho*u_n[i]*(mt.sqrt(4*A_n[i]/(mt.pi))))/d_vis
        if Re < 2000:
            f=64/Re
        else :
            f=0.3164*Re**-0.25 
        a[i]=((rho*Area[2*i]/dt)+(f*rho*u_n[i]*p_s[2*i]/8))
        b[i]=((rho*Area[2*i]/(2*dx)))
        d[i]=((rho*A_n[2*i]*u_n[i]/dt))
        if i<Grid_points:
            c[i]=(Area[(2*i)+1]/dx)
    return a,b,c,d 
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
        z[i]=mu[i]-(lamda[i+1]*u_star[i+1])+(lamda[i]*u_star[i])       
    return alpha,gama,beta,z    
def convergence(u_star,Grid_points,rho,Area,A_n,dx,dt):
    lamda,m_u=get_lamda_mu(Grid_points,rho,Area,A_n,dx,dt)
    residue=np.zeros((Grid_points-1))
    for i in range(0,len(residue)):
        residue[i]=(lamda[i+1]*u_star[i+1])-(lamda[i]*u_star[i])+(m_u[i])
    mean= np.mean(np.abs(residue))
    return mean
def lu_decomposition(A):
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))  
    for i in range(n):
        L[i, i] = 1.0
        for j in range(i, n):
            U[i, j] = A[i, j] - np.dot(L[i, :i], U[:i, j])
        for j in range(i+1, n):
            L[j, i] = (A[j, i] - np.dot(L[j, :i], U[:i, i])) / U[i, i]
    return L, U
def solve_lu(A, b):
    L, U = lu_decomposition(A)
    n = A.shape[0]
    y = np.zeros((n, 1))
    x = np.zeros((n, 1))

    # Solve Ly = b using forward substitution
    for i in range(n):
        y[i, 0] = b[i, 0] - np.dot(L[i, :i], y[:i, 0])

    # Solve Ux = y using backward substitution
    for i in range(n - 1, -1, -1):
        x[i, 0] = (y[i, 0] - np.dot(U[i, i+1:], x[i+1:, 0])) / U[i, i]
    return x

def solve_momentum(a, b, c, d, p_star,Grid_points,u_entry,p_exit,l1_extra_cell,mu_extra_cell,l2_extra_cell):
    n =Grid_points
    if n>=2:
        A = np.zeros((n, n+1))
        Rhs=np.zeros((n,1))
        for i in range(1,n+1):
            if i<n:
                A[i-1, i-1] =b[i-1]
                A[i-1,i] = (-1*a[i])
                A[i-1,i+1]=(-1*b[i+1])
                const = ((c[i] * p_star[i]) - (c[i-1] * p_star[i-1]) - d[i])
                Rhs[i-1,0]=const
            else:
                A[i-1, i-1] = l1_extra_cell
                A[i-1,i] = -1*l2_extra_cell
                Rhs[i-1,0]=mu_extra_cell
        Rhs[0,0]= Rhs[0,0]-(A[0,0]*u_entry)
        coefficient_matrix= A[:,1:n+1]
        u_star1=solve_lu(coefficient_matrix,Rhs)
    return u_star1
def Combine(u_star,u_star1,Grid_points):
    u_star1=u_star1.reshape(1,Grid_points)
    u_star[1:Grid_points+1]=u_star1[:, :]
    return u_star
def Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit):
    a_values,b_values,c_values,d_values=get_abcd(Grid_points,rho,dt,dx,d_vis,A_n,Area,u_n,p_s)    
    l1_extra_cell,mu_extra_cell,l2_extra_cell=Outlet_Boundary_const(rho,A_n,Area,dt,dx,Grid_points)
    u_star1= solve_momentum(a_values, b_values, c_values, d_values, p_star,Grid_points,u_inlet,p_exit,l1_extra_cell,mu_extra_cell,l2_extra_cell)
    u_star=Combine(u_star,u_star1,Grid_points)
    return u_star
def correct_pressure(Grid_points,alpha,beta,gama,z,p_exit,c_1,c_2,z1):
    n=Grid_points
    A1 = np.zeros((n-1, n))
    Rhs1=np.zeros((n-1,1))
    A1[0,0] =c_1
    A1[0,1] =-1*c_2
    Rhs1[0,0]=z1
    for i in range(1,n-1):
        A1[i,i-1] = gama[i-1]
        A1[i,i] = (-1*alpha[i])
        A1[i,i+1]= (beta[i+1])
        const = (-1*z[i])
        Rhs1[i,0]=const
    Rhs1[n-2,0]=Rhs1[n-2,0]-(A1[n-2,n-1]*p_exit)
    co_efficientmatrix1=A1[:,0:n-1] 
    p_=solve_lu(co_efficientmatrix1,Rhs1)    
    return p_
def Pressure_adjust(u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit):
    alpha,beta,gama,z=get_pressure(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis) 
    c_1,c_2,z1=Inlet_Boundary_equation(Grid_points,u_n,u_star,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
    add_p1=correct_pressure(Grid_points,alpha,beta,gama,z,p_exit,c_1,c_2,z1)#Solve pressure correction equation 
    return add_p1   
def add(p_add,p_star,Grid_points):
    relaxation_factor=0.001
    p_add=p_add.reshape(Grid_points-1)
    p_star[0:Grid_points-1] += relaxation_factor*p_add[0:Grid_points-1]
    return p_star
# Main Algo
def unsteady_1D_flow(A_n,A,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit):
    p_star = np.zeros(Grid_points)
    u_star = np.zeros(Grid_points+1) 
    """
    Initialize the pressure 
    """
    p_star[len(p_star)-1]=p_exit
    for i in range(len(p_star)-2,-1,-1):
        p_star[i]=p_star[i+1]+((dp*dx)/length)
    
    u_star=u_n
    Area=np.full((2 * Grid_points + 1), A)

    for t in range(0,n):
        start_time=time.time()
        def Graph_P():
            fig=plt.figure()
            ax2 = fig.add_subplot(111)
            return fig, ax2 
        def plot_P(x,y,graph,labelx,labely):
            graph[1].plot(x,y,alpha=1) 
            graph[1].set_xlabel(labelx)
            graph[1].set_ylabel(labely)
        graph_p=Graph_P()
        x2=np.arange(0,len(p_star),1)
        while True:                                                                        
            u_star=Momentum(p_star,u_n,u_star,Grid_points,rho,dt,dx,d_vis,A_n,Area,p_s,u_inlet,p_exit)
            converge=convergence(u_star,Grid_points,rho,Area,A_n,dx,dt)
            plot_P(x2,p_star,graph_p,"Along x","p_star")
            plt.pause(0.01)
            p_add=Pressure_adjust(u_star,Grid_points,u_n,Area,A_n,p_s,rho,dx,dt,d_vis,u_inlet,p_exit)
            p_star=add(p_add,p_star,Grid_points)
            elasp_time = time.time() - start_time 
            if elasp_time >10:
                break 
        u_n=u_star
    plt.show()    
    return u_star,p_star

"""
Final result 
"""
velocity, Pressure =unsteady_1D_flow(A_n,A,u_n,p_s,Grid_points,rho,dx,dt,d_vis,n,u_inlet,p_exit)
