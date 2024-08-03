
import numpy as np
from math import pi

def stenoisis_Area5(dia, Grid_points, Length):
  L = Length
  n = 2*Grid_points
  dx = L / n
  r_i = dia/2
  r_cs = 0.00166
  r_o = 0.00138
  x_i = 0
  x_o = L
  x_is = 0.02
  x_cs = 0.035
  x_os = 0.05
  r = np.zeros(n + 1)
  x = np.zeros(n + 1)
  x_star = np.zeros(n + 1)
  theta = np.zeros(n + 1)
  b = 0.4
  for i in range(n + 1):
    x[i] = dx * i
    r[i] = r_i + (x[i] - x_i) * (r_o - r_i) / (x_o - x_i)
    if x[i] == x_is:
      r_is = r[i]
    if x[i] == x_os:
      r_os = r[i]
  for i in range(n + 1):
    if (x[i] >= x_is) and (x[i] <= x_os):
      x_star[i] = (x[i] - x_cs) / (x_os - x_is)
      theta[i] = x_star[i] * 2 * pi
      r[i] = r[i] * (1 - b * 0.5 * (1 + np.cos(theta[i])))

  Area = pi * r**2
  p_s = 2 * pi * r
  return r, Area, p_s, x

