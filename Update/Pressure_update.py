import numpy as np

def add(p_add,p_star,relaxation_factor):
    p_star[0:len(p_star)] += relaxation_factor*p_add[0:len(p_add)]
    return p_star