import numpy as np
import time
import threading
import matplotlib.pyplot as plt



def Graph_PV(u_star,p_star):
    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    ax1.set_xlabel('Grid Points')
    ax1.set_ylabel('u_star')
    ax1.legend()
    ax2.set_xlabel('Grid Points')
    ax2.set_ylabel('p_star')
    ax2.legend()
    plt.tight_layout()
    return fig, ax1, ax2
    
def plot(u_star,p_star,x,graph):
    graph[1].plot(x, u_star, label='u_star',alpha=0.8)
    graph[2].plot(x, p_star, label='p_star',alpha=0.8)