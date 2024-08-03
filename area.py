import csv
import os
import numpy as np
import time 
import math as mt
import matplotlib.pyplot as plt
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



Grid_points = 400
length = 0.07
dia = 0.0039
r_x, A_n, p_s, x5 = stenoisis_Area5(dia, Grid_points, length)

dx = (100 * length) / (2 * Grid_points)
L1 = 2
L2 = 5

n = 25  # No of segments
length_of_stenoised_region = 3
size = length_of_stenoised_region / (n - 2)
L1 = 2  # Starting point
L21 = 2
data = []  # List to hold the data for CSV

# Add A_n[0] at the beginning
data.append([A_n[0] * 10 ** 4, 0, L1])

for i in range(0, n - 1):
    n11 = int((L1 + (i * size)) * 2 * Grid_points / (length * 100))
    data.append([A_n[n11] * 10 ** 4, L1 + (i * size), size])
    print(A_n[n11] * 10 ** 4, "A_n at ", (L1 + (i * size)))

# Add A_n[len(A_n)-1] at the end
data.append([A_n[-1] * 10 ** 4, length * 100, L21])

# Writing data to CSV
save_directory = r'D:\Simvascular\S_40\Models'
os.makedirs(save_directory, exist_ok=True)
csv_file_path = os.path.join(save_directory, 'output.csv')
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['A_n', 'L1+(i*size)', 'length_segment'])
    writer.writerows(data)

print("CSV file has been generated successfully!")
print("Done")




































"""
import numpy as np
import os
import time 
import math as mt
import matplotlib.pyplot as plt
from Area.linear_area import stenoisis_Area5
import csv
"""


"""


Grid_points = 400
length = 0.07
dia = 0.0039
r_x, A_n, p_s, x5 = stenoisis_Area5(dia, Grid_points, length)

dx = (100 * length) / (2 * Grid_points)
L1 = 2
L2 = 5
n1 = int(L1 / dx)
n2 = int(L2 / dx)

# Calculate the total number of points needed
total_points_needed = 50

# Calculate the number of points to be generated between n1 and n2
points_between_n1_n2 = total_points_needed - 2  # Subtract 2 for n1 and n2

# Generate equally spaced points between n1 and n2
x_points = np.linspace(n1+1, n2-1, points_between_n1_n2)

# Round the points to integers
x_points = np.round(x_points).astype(int)

# Add n1 and n2 to the list of points
x_points = np.concatenate(([n1], x_points, [n2]))

# Print the total number of points
print(len(x_points), "Total number of grid points between n1 and n2:")

# Print the points
print("Grid points:", x_points)

Length_indices= (x_points*dx)+0.005
# Access A_n values at x_points indices
area_values = (A_n[x_points]*(10**4))

# Print the area values
print("Area values at the generated points:", area_values)
print("Length indices:", Length_indices)


# Write to CSV file
with open('area_length_indices.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Area', 'Length_indices'])
    for area, length_index in zip(area_values, Length_indices):
        writer.writerow([area, length_index])

print("Data saved to 'area_length_indices.csv'")
print(A_n[0]*10**4, "A_n[0]")
print(A_n[len(A_n)-1]*10**4, "A_n[len(A_n)-1]")

"""

