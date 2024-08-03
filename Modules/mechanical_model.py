import numpy as np
import math as mt
#from Interpolate.inter import interpolate_Area



Ao=0.22*(10**-4)
r_d=3*10**-3
A_d=mt.pi*(r_d**2)
E=700*10**3
h=0.3*10**-3
beta=4/3*(mt.sqrt(mt.pi))*E*h
p_d=10900
rho=1000


def area(Pressure):
    A=((np.sqrt(A_d))+((Pressure-p_d)*A_d/beta))**2
    #A=interpolate_Area(A)
    #c=np.sqrt(beta/(2*rho*A_d))*(A)**0.25
    #print(c,"c")
    return A


"""
x=area(np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7]))
print(x)"""