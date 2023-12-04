import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def Graph_V_w_time():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_ylim(0, 1.0)
    return fig, ax

    
def plot_V_time(u_star,x,time,graph,colour):
    #graph[1].plot(x, time, u_star)
    #x, time = np.meshgrid(x, time)
     # Replace this with your actual value of t1
    Time=time*np.ones_like(x)
    graph[1].plot(x, Time, u_star,linewidth=1)
    #graph[1].plot(u_star, x, time, marker='o')  # Plotting V vs Time in 3D
    graph[1].set_xlabel('x',fontsize=12)
    graph[1].set_ylabel('time',fontsize=12)
    graph[1].set_zlabel('u_star',fontsize=12)
    # Set title
    graph[1].set_title('3D Line Plot', fontsize=16)
    #plt.show()