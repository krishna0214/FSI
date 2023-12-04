import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def Graph_V_w_time():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_ylim(0, 1.0)
    return fig, ax

    
#def plot_V_time(u_star,x,time,graph,colour):
    #graph[1].plot(x, time, u_star)
    #x, time = np.meshgrid(x, time)
     # Replace this with your actual value of t1
    Time=time*np.ones_like(x)
    graph[1].plot(x, Time, u_star,cmap='plasma',linewidth=1)
    #graph[1].plot(u_star, x, time, marker='o')  # Plotting V vs Time in 3D
    graph[1].set_xlabel('x',fontsize=12)
    graph[1].set_ylabel('time',fontsize=12)
    graph[1].set_zlabel('u_star',fontsize=12)
    # Set title
    graph[1].set_title('3D Line Plot', fontsize=16)
    #plt.show()

#def plot_V_time(u_star, x, time, graph, colour):
    Time = time * np.ones_like(x)
    
    # Creating meshgrid for surface plot
    X, Y = np.meshgrid(x, Time)
    
    # Replace this line with your actual value of u_star
    Z = u_star * np.ones_like(X)

    # Plotting the surface with colormap
    graph[1].plot_surface(X, Y, Z, cmap='plasma', edgecolor='none')
    
    # Set labels and title
    graph[1].set_xlabel('x', fontsize=12)
    graph[1].set_ylabel('time', fontsize=12)
    graph[1].set_zlabel('u_star', fontsize=12)
    graph[1].set_title('3D Surface Plot', fontsize=16)
    
    # Adding colorbar for reference
    #fig.colorbar(surf, ax=graph[1], shrink=0.5, aspect=5)


def plot_V_time(u_star, x, time,graph):
    X, T = np.meshgrid(x, time)  # Create meshgrid for X and time
    # Ensure u_star has the same shape as X and T
  # Repeat u_star along time axis
    # Plotting the surface with colormap
    #graph[1].plot_surface(X, T,u_star)
    surf = graph[1].plot_surface(X, T, u_star, cmap='viridis', edgecolor='none')
    graph[1].figure.colorbar(surf, ax=graph[1], shrink=0.5, aspect=5)
    #graph[1].plot_wireframe(X, T, U_star, color='blue')
    # Creating the surface by plotting lines with varying heights
    #for i in range(len(x)):
        #graph[1].plot([x[i], x[i]], [time[0], time[0]], [0, U_star[i, 0]], color='blue')

    # Set labels and title
    graph[1].set_xlabel('x', fontsize=12)
    graph[1].set_ylabel('time', fontsize=12)
    graph[1].set_zlabel('u_star', fontsize=12)
    graph[1].set_title('3D Surface Plot', fontsize=16)

    # Adding colorbar for reference
    #graph[0].colorbar(surf, ax=graph[1], shrink=0.5, aspect=5)

