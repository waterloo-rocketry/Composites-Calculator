import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi

Q0=np.array([
			[155.7*1e7,3.02*1e7,0*1e7],
			[3.02*1e7,12.16*1e7,0*1e7],
			[0*1e7,0*1e7,4.40*1e7]
	])

#Function returns transformed stiffness matrix Qbar
def transform_Q(Q0,ply_angle):
    m=cos(ply_angle)
    n=sin(ply_angle)
    T1=[[m**2,n**2,2*m*n],[n**2,m**2,-2*m*n],[-m*n,m*n,m**2-n**2]]
    T2=[[m**2,n**2,m*n],[n**2,m**2,-m*n],[-2*m*n,2*m*n,m**2-n**2]]
    T1_inv=np.linalg.inv(T1)
    Qbar=np.linalg.multi_dot([T1_inv,Q0,T2])
   
    return Qbar