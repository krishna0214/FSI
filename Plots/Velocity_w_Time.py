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


def plot_V_time(u_star, x, time,graph,title):
    X, T = np.meshgrid(x, time)  # Create meshgrid for X and time
    # Ensure u_star has the same shape as X and T
  # Repeat u_star along time axis
    # Plotting the surface with colormap
    #graph[1].plot_surface(X, T,u_star)
    #surf = graph[1].plot_surface(X, T, u_star, cmap='viridis', edgecolor='none')
    min=np.min(u_star)
    max=np.max(u_star)
    graph[1].set_box_aspect([1, 1, 0.5])
    surf = graph[1].plot_surface(X, T, u_star, cmap='coolwarm', edgecolor='none')
    surf.set_clim(min, max)
    graph[1].view_init(elev=30, azim=45)
    graph[1].figure.colorbar(surf, ax=graph[1], shrink=0.5, aspect=5)
    #graph[1].plot_wireframe(X, T, U_star, color='blue')
    # Creating the surface by plotting lines with varying heights
    #for i in range(len(x)):
        #graph[1].plot([x[i], x[i]], [time[0], time[0]], [0, U_star[i, 0]], color='blue')
    graph[1].set_zlim(min, max) 
    # Set labels and title
    graph[1].set_xlabel('x', fontsize=12)
    graph[1].set_ylabel('time', fontsize=12)
    graph[1].set_zlabel(title, fontsize=12)
    graph[1].set_title(title, fontsize=16)

    # Adding colorbar for reference
    #graph[0].colorbar(surf, ax=graph[1], shrink=0.5, aspect=5)

def plot_P_time(u_star, x, time, graph, title):
    X, T = np.meshgrid(x, time)  # Create meshgrid for X and time
    
    # Determine the minimum and maximum values of u_star for proper scaling
    z_min = 100000
    z_max = np.max(u_star)

    # Plotting the surface with colormap and adjusted z-axis limits
    surf = graph[1].plot_surface(X, T, u_star, cmap='coolwarm', edgecolor='none')
    surf.set_clim(z_min, z_max)  # Set the color limits based on min and max values

    graph[1].view_init(elev=30, azim=45)
    graph[1].figure.colorbar(surf, ax=graph[1], shrink=0.5, aspect=5)

    # Set z-axis limits for proper scaling
    graph[1].set_zlim(z_min, z_max)  # Set the z-axis limits based on min and max values
    # Assuming 'ax' is your 3D Axes
    graph[1].xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))  # Set the x pane color to transparent
    graph[1].yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))  # Set the y pane color to transparent
    graph[1].zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))  # Set the z pane color to transparent

    # Set labels and title
    graph[1].set_xlabel('x', fontsize=12)
    graph[1].set_ylabel('time', fontsize=12)
    graph[1].set_zlabel(title, fontsize=12)
    graph[1].set_title(title, fontsize=16)
