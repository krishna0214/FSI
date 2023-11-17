import math as mt
import numpy as np
from Plots.Velocity_plot import Graph_V,plot_V
import matplotlib.pyplot as plt


def Area(dia,Grid_points,length,l1,l2):
    d_in=dia
    r_in=dia/2
    H=dia/100
    W=(l2-l1)
    dx=length/(2*Grid_points)
    Area_in=mt.pi*(d_in**2)/4
    p_in=mt.pi*d_in
    r_x1=np.full((2*Grid_points+1),d_in/2)
    r_x2=np.full((2*Grid_points+1),d_in/2)
    A_n=np.full((2*Grid_points+1),Area_in)
    p_s=np.full((2*Grid_points+1),p_in)
    n1=l1
    n2=l2
    for i in range (n1,n2):
        center=((n2+n1)*length)/(2*Grid_points)
        z=(center)-((length*i)/(2*Grid_points))
        r_x1[i]= r_in+(H*(mt.exp((-z**2)/(2*(W**2)))))
        r_x2[i]= -1*(r_in(H*(mt.exp((-z**2)/(2*(W**2))))))
        A_n[i]=mt.pi*(r_x1[i]**2)
        p_s[i]=mt.pi*2*r_x2[i]
    return r_x1,r_x2   




    




dia=0.0154
length=2
Grid_points=200
H=dia/100
W=length/4
dx=length/Grid_points
r_effective=np.zeros((2*Grid_points+1))

for i in range(0,2*Grid_points+1):
    z= (length/2)-(i*dx/2)
    r_effective[i]=(dia/2)+(H*mt.exp(-z**2/(2*(W**2))))

graph_v=Graph_V()
x2=np.arange(len(r_effective))
plot_V(r_effective,x2,graph_v)



H=dia/100
W=length/4
A_n = np.zeros(2 * Grid_points + 1)  
p_s=np.zeros(2*Grid_points+1)              # Area at n_th step
r_effective=np.zeros((2*Grid_points+1))
for i in range(0,2*Grid_points+1):
    z= (length/2)-(i*dx/2)
    r_effective[i]=(dia/2)+(H*mt.exp(-z**2/(2*(W**2))))
    A= (mt.pi)*(r_effective[i]**2)
    A_n[i]=A
    perimeter= (mt.pi)*2*r_effective[i]
    p_s[i]=perimeter
