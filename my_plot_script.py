import numpy as np
import os
import time
import math as mt
import matplotlib.pyplot as plt
from ipywidgets import interact, widgets
from Modules.Momentum import Momentum
from Modules.Pressure import Pressure_adjust
from Update.Pressure_update import add
from Modules.convergence import convergence
from Plots.Plot import Graph_PV, plot
from Plots.Velocity_plot import Graph_V, plot_V
from Plots.Pressure_plot import Graph_P, plot_P

# ... (your existing code)

# Define a function to update the plot based on dt and dx
def update_plot(dt, dx):
    Total_time = 1
    n = 1
    rho = 996.550
    g = 9.8
    Grid_points = 100
    Inletmassflux = 900
    P_atm = 1
    dia = 0.0154
    d_vis = 0.000854
    p_exit = P_atm * 101325
    A = (mt.pi) * (dia**2) / 4
    u_inlet = Inletmassflux / (rho)
    length = 2


    u_n = np.zeros(Grid_points)
    p_star = np.zeros(Grid_points)
    A_n = np.full((2 * Grid_points + 1), A)
    perimeter = (mt.pi) * dia
    p_s = np.full((2 * Grid_points + 1), perimeter)
    dp = 128 * d_vis * length * A * u_inlet / ((mt.pi) * dia**4)

    p_star[len(p_star) - 1] = p_exit
    for i in range(len(p_star) - 2, -1, -1):
        p_star[i] = p_star[i + 1] + ((dp * dx) / length)

    u_n[0] = u_inlet
    for i in range(1, len(u_n)):
        u_n[i] = 0.01

    def unsteady_1D_flow(A_n, A, u_n, p_star, p_s, Grid_points, rho, dx, dt, d_vis, n, u_inlet, p_exit):
        Area = np.full((2 * Grid_points + 1), A)
        for t in range(0, n):
            start_time = time.time()
            i = 0
            while True:
                u_star = Momentum(p_star, u_n, Grid_points, rho, dt, dx, d_vis, A_n, Area, p_s, u_inlet, p_exit)
                converge1 = convergence(u_n, Grid_points, rho, Area, A_n, dx, dt)
                converge = convergence(u_star, Grid_points, rho, Area, A_n, dx, dt)
                if i > 0:
                    break
                p_add = Pressure_adjust(p_star, u_star, Grid_points, u_n, Area, A_n, p_s, rho, dx, dt, d_vis, u_inlet,p_exit)
                p_star = add(p_star, p_add, Grid_points)
                i = i + 1
                elasp_time = time.time() - start_time
                if elasp_time > 1000:
                    break
            u_n = u_star
        return u_star, p_star

    velocity, Pressure = unsteady_1D_flow(A_n, A, u_n, p_star, p_s, Grid_points, rho, dx, dt, d_vis, n, u_inlet, p_exit)

    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    ax1.set_xlabel('Grid Points')
    ax1.set_ylabel('u_star')
    ax2.set_xlabel('Grid Points')
    ax2.set_ylabel('p_star')
    ax1.plot(velocity, label='Velocity')
    ax2.plot(Pressure, label='Pressure')
    plt.tight_layout()
    plt.show()

# Create interactive widgets for dt and dx
dt_widget = widgets.FloatSlider(value=0.0001, min=0.0001, max=0.001, step=0.0001, description='dt:')
dx_widget = widgets.FloatSlider(value=0.02, min=0.01, max=0.1, step=0.01, description='dx:')

# Use interact to connect the widgets with the update_plot function
interact(update_plot, dt=dt_widget, dx=dx_widget)
