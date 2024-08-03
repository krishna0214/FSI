import numpy as np
from stl import mesh

def stenosis_Area(dia, Grid_points, Length):
    L = Length
    n = 2 * Grid_points
    dx = L / n
    r_i = dia / 2
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
    b = 0.2
    
    for i in range(n + 1):
        x[i] = dx * i
        r[i] = r_i + (x[i] - x_i) * (r_o - r_i) / (x_o - x_i)
        
        if x[i] == x_is:
            r_is = r[i]
        
        if x[i] == x_os:
            r_os = r[i]
    
    for i in range(n + 1):
        if x_is <= x[i] <= x_os:
            x_star[i] = (x[i] - x_cs) / (x_os - x_is)
            theta[i] = x_star[i] * 2 * np.pi
            r[i] = r[i] * (1 - b * 0.5 * (1 + np.cos(theta[i])))
    
    Area = np.pi * r**2
    p_s = 2 * np.pi * r
    return r, Area, p_s, x

# Parameters
dia = 0.02  # Diameter
Grid_points = 50  # Number of grid points
Length = 1.0  # Length of the tube

# Calculate area profile
r, Area, _, x = stenosis_Area(dia, Grid_points, Length)

# Create a 3D mesh
vertices = np.column_stack((x, np.zeros_like(x), Area))
num_vertices = len(x)
num_faces = num_vertices - 1

faces = np.zeros((num_faces, 3), dtype=int)
for i in range(num_faces):
    faces[i][0] = i
    faces[i][1] = i + 1
    faces[i][2] = i + num_vertices

faces = np.tile(faces, (2, 1))

# Create the mesh object
stl_mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        stl_mesh.vectors[i][j] = vertices[f[j], :]

# Save the mesh to an STL file
stl_mesh.save('area_profile.stl')
