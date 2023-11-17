import math as mt
import numpy as np
from Plots.Velocity_plot import Graph_V,plot_V
import matplotlib.pyplot as plt

dia=0.0154
length=2
Grid_points=200
H=0.5*dia
W=length/4
dx=length/Grid_points
r_effective=np.full((2*Grid_points+1),dia/2)
#A_n = np.zeros(2*Grid_points + 1)  
#p_s=np.zeros(2*Grid_points+1)


A= (mt.pi)*r_effective**2
perimeter= (mt.pi)*2*r_effective
A_n=np.full((2*Grid_points+1),A)
p_s=np.full((2*Grid_points+1),perimeter)

l1=length/4
l2=3*length/4

for i in range(int(Grid_points/4),int(3*Grid_points/4)):
    z= ((l2-l1)/2)-(i*dx/2)
    r_effective[i]=(dia/2)+(H*mt.exp(-z**2/(2*(W**2))))
    A_n[i]=(mt.pi)*(r_effective[i]**2)
    p_s[i]=(mt.pi)*2*r_effective[i]


graph_v=Graph_V()
x2=np.arange(len(r_effective))
plot_V(A_n,x2,graph_v)
plt.show()






