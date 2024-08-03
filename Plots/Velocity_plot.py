import numpy as np
import matplotlib.pyplot as plt


def Graph_V():
    fig=plt.figure(figsize=(5, 4),dpi=200)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('x')
    ax1.set_ylabel('velocity in m/s')
    ax1.legend()
    plt.tight_layout()
    return fig, ax1
    
def plot_V(u_star,x,x2,velocity,graph):   
    graph[1].plot(x, u_star, label='Velocity',alpha=0.8,linewidth=1.5)
    graph[1].scatter(x2,velocity,color='red',alpha=1.0)