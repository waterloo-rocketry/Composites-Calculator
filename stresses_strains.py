import numpy as np
from math import cos, sin, pi
from Q_and_Qbar import transform_Q

#Function returns the stress in primary material directions as a column vector
def transform_T1(principalStress,ply_angle):
    m=cos(ply_angle)
    n=sin(ply_angle)
    T1=[[m**2,n**2,2*m*n],[n**2,m**2,-2*m*n],[-m*n,m*n,m**2-n**2]]
    primaryStress=np.dot(T1,principalStress)
    
    return primaryStress


#Calculate 2 values of strain and stress for every height element, excluding 
#the first and last heights which only result in 1 value of strain and stress because they're on the surface
#There's an opportunity to remove the inner loop but when I try to make height a column vector the dimensions are mismatched
def global_ply(Q0, stack, estrain_glob, kstrain_glob,height):
	plystrain_glob=np.zeros(shape=(3,2*len(height)-2))
	plystress_glob=np.zeros(shape=(3,2*len(height)-2))

	for k in range(len(height)-1):
	    for i in range(3):
	        plystrain_glob[i,2*k]=estrain_glob[i]+height[k]*kstrain_glob[i]
	        plystrain_glob[i,2*k+1]=estrain_glob[i]+height[k+1]*kstrain_glob[i]
	    Qbar=transform_Q(Q0,stack[k])
	    plystress_glob[:,2*k]=np.dot(Qbar,plystrain_glob[:,2*k])
	    Qbar=transform_Q(Q0,stack[k])
	    plystress_glob[:,2*k+1]=np.dot(Qbar,plystrain_glob[:,2*k+1])

	return (plystrain_glob, plystress_glob)

def local_ply(plystrain_glob, plystress_glob, height, stack):
	plystrain_loc=np.zeros(shape=(3,2*len(height)-2))
	plystress_loc=np.zeros(shape=(3,2*len(height)-2))

	for k in range(len(height)-1):
	    plystrain_loc[:,2*k]=transform_T1(plystrain_glob[:,2*k],stack[k])
	    plystrain_loc[:,2*k+1]=transform_T1(plystrain_glob[:,2*k+1],stack[k])
	    plystress_loc[:,2*k]=transform_T1(plystress_glob[:,2*k],stack[k])
	    plystress_loc[:,2*k+1]=transform_T1(plystress_glob[:,2*k+1],stack[k])

	return (plystrain_loc, plystress_loc)