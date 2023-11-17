import numpy as np
import matplotlib.pyplot as plt


def Graph_Area_profile():
    fig=plt.figure()
    ax1=fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3=fig.add_subplot(313)
    ax1.set_xlabel('Along z')
    ax1.set_ylabel('Radius')
    ax1.legend()
    ax2.set_xlabel('Along z')
    ax2.set_ylabel('Area')
    ax2.legend()
    ax3.set_xlabel('Along z')
    ax3.set_ylabel('Perimeter')
    ax3.legend()
    plt.tight_layout()
    return fig, ax1,ax2,ax3
    
def plot_Area(x,r1,r2,Area,Perimeter,graph):
    graph[1].plot(x, r1, label='Radius',alpha=0.8)
    graph[1].plot(x, r2, label='Radius',alpha=0.8)
    graph[2].plot(x, Area, label='Area',alpha=0.8)
    graph[3].plot(x, Perimeter, label='Perimeter',alpha=0.8)