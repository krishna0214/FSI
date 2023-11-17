import numpy as np
import matplotlib.pyplot as plt


def Graph_V():
    fig=plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Grid Points')
    ax1.set_ylabel('u_star')
    ax1.legend()
    plt.tight_layout()
    return fig, ax1
    
def plot_V(u_star,x,graph):
    graph[1].plot(x, u_star, label='u_star',alpha=0.8)