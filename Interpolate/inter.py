import numpy as np
from scipy.interpolate import interp1d

def interpolate_pressure(p_star):
    n = len(p_star)
    x = np.arange(n)  # Generate x values for interpolation
    f = interp1d(x, p_star, kind='quadratic')  # Linear interpolation function
    x_interpolated = np.linspace(0, n - 1, 2 * n - 1)  # Generate interpolated x values
    interpolated_pressure = f(x_interpolated)  # Interpolate pressure values
    pressure_inlet_Boundary=(2*interpolated_pressure[0]-interpolated_pressure[1])
    pressure_outlet_Boundary=(2*interpolated_pressure[-1]-interpolated_pressure[-2])
    interpolated_pressure = np.insert(interpolated_pressure, 0, pressure_inlet_Boundary)
    interpolated_pressure = np.append(interpolated_pressure, pressure_outlet_Boundary)
    return interpolated_pressure


def interpolate_Area(Area):
    n = len(Area)
    x = np.arange(n)  # Generate x values for interpolation
    f = interp1d(x, Area,kind='quadratic')  # Linear interpolation function
    x_interpolated = np.linspace(0, n - 1, 2 * n - 1)  # Generate interpolated x values
    interpolated_Area = f(x_interpolated)  # Interpolate area values
    Area_inlet_Boundary = 2 * interpolated_Area[0] - interpolated_Area[1]
    Area_outlet_Boundary = 2 * interpolated_Area[-1] - interpolated_Area[-2]
    interpolated_Area = np.insert(interpolated_Area, 0, Area_inlet_Boundary)
    interpolated_Area = np.append(interpolated_Area, Area_outlet_Boundary)
    return interpolated_Area
