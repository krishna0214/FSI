import numpy as np
import matplotlib.pyplot as plt


def Graph_shear():
    fig=plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('x')
    ax1.set_ylabel('Shear_avg')
    ax1.legend()
    plt.tight_layout()
    return fig, ax1
    
def plot_shear(shear,x,graph):
    graph[1].plot(x, shear, label='u_star',alpha=0.8)