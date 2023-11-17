import math as mt
import numpy as np
from Plots.Velocity_plot import Graph_V,plot_V
import matplotlib.pyplot as plt

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