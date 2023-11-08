def add(p_star,p_add,Grid_points):
    p_add=p_add.reshape(Grid_points-1)
    p_star[0:Grid_points-1] += p_add[0:Grid_points-1]
    return p_star