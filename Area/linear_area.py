import numpy as np
import math as mt




def Linear_area_profile(dia ,length,Grid_points):
    d_in=dia
    d_out=dia/2
    Area_in=mt.pi*(d_in**2)/4
    p_in=mt.pi*d_in
    r_x=np.full((2*Grid_points+1),d_in/2)
    A_n=np.full((2*Grid_points+1),Area_in)
    p_s=np.full((2*Grid_points+1),p_in)
    slope=(d_out-d_in)/(2*length)
    for i in range(1,len(A_n)):
        r_x[i]=r_x[i-1]+((slope*(length/Grid_points))/2)
        A_n[i]=mt.pi*(r_x[i]**2)
        p_s[i]=mt.pi*2*r_x[i]
    return r_x,A_n,p_s
    


def Const_area_profile(dia ,length,Grid_points):
    d_in=dia
    Area_in=mt.pi*(d_in**2)/4
    p_in=mt.pi*d_in
    r_x=np.full((2*Grid_points+1),d_in/2)
    A_n=np.full((2*Grid_points+1),Area_in)
    p_s=np.full((2*Grid_points+1),p_in)
    return r_x,A_n,p_s


def Area(dia,Grid_points,length,l1,l2):
    d_in=dia
    r_in=dia/2
    H=dia/10
    W=(l2-l1)/10
    dx=length/(2*Grid_points)
    Area_in=mt.pi*(d_in**2)/4
    p_in=mt.pi*d_in
    r_x1=np.full((2*Grid_points+1),d_in/2)
    r_x2=np.full((2*Grid_points+1),-1*d_in/2)
    A_n=np.full((2*Grid_points+1),Area_in)
    p_s=np.full((2*Grid_points+1),p_in)
    n1=l1/dx
    n2=l2/dx
    n1=int(n1)
    n2=int(n2)

    for i in range (n1,n2):
        #center=((n2+n1)*length)/(2*Grid_points)
        #z=(center)-((length*i)/(2*Grid_points))
        z=(((n1+n2)/2)-i)*(dx)
        r_x1[i]= r_in+(H*(mt.exp((-z**2)/(2*(W**2)))))
        r_x2[i]= -1*(r_in+(H*(mt.exp((-z**2)/(2*(W**2))))))
        A_n[i]=mt.pi*(r_x1[i]**2)
        p_s[i]=mt.pi*2*r_x1[i]
    return r_x1,r_x2,A_n,p_s   

