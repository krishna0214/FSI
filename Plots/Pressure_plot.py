import numpy as np
import matplotlib.pyplot as plt


def Graph_P():
    fig=plt.figure(figsize=(8, 6),dpi=200)
    ax2 = fig.add_subplot(111)
    ax2.set_xlabel('Grid Points')
    ax2.set_ylabel('p_star')
    ax2.legend()
    plt.tight_layout()
    return fig, ax2
    
def plot_P(x,p_star,graph):
    graph[1].plot(x, p_star, label='p_star',alpha=1)

    