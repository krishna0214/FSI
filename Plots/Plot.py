import numpy as np
import matplotlib.pyplot as plt



def Graph_PV():
    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    ax1.set_xlabel('Grid Points')
    ax1.set_ylabel('p_star')
    ax2.set_xlabel('Grid Points')
    ax2.set_ylabel('u_star')
    plt.tight_layout()
    return fig, ax1, ax2
    
def plot(p_star,u_star,x1,x2,graph):
    
    graph[1].plot(x1, p_star,alpha=0.8)
    #graph[1].legend()
    graph[2].plot(x2, u_star,alpha=0.8)
    #graph[2].legend()



def Graph_P():
    fig=plt.figure()
    ax2 = fig.add_subplot(111)
    return fig, ax2 
def plot_P(x,y,graph,labelx,labely):
    graph[1].plot(x,y,alpha=1) 
    graph[1].set_xlabel(labelx)
    graph[1].set_ylabel(labely)

