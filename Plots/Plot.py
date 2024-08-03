import numpy as np
import matplotlib.pyplot as plt



def Graph_PV():
    fig=plt.figure()
    ax1=fig.add_subplot(311)
    ax2=fig.add_subplot(312)
    ax3=fig.add_subplot(313)
    ax1.set_xlabel('Grid Points')
    ax1.set_ylabel('velocity')
    ax2.set_xlabel('Grid Points')
    ax2.set_ylabel('pressure')
    ax3.set_xlabel('Grid Points')
    ax3.set_ylabel('Area')
    
    plt.tight_layout()
    return fig, ax1, ax2, ax3
    
def plot(u_star,p_star,Area,x1,x2,x3,graph):
    
    graph[1].plot(x1, u_star,alpha=0.8)
    #graph[1].legend()
    graph[2].plot(x2, p_star,alpha=0.8)
    #graph[2].legend()
    graph[3].plot(x3, Area,alpha=0.8)



def Graph_P():
    fig=plt.figure()
    ax2 = fig.add_subplot(111)
    return fig, ax2 
def plot_P(x,y,graph,labelx,labely):
    graph[1].plot(x,y,alpha=1) 
    graph[1].set_xlabel(labelx)
    graph[1].set_ylabel(labely)

