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

##Constants (will probably be moved back to main)
#Laminate parameters
stack=[0,pi/2,pi/2,0] #Stacking sequence from top down
thickness=[0.15*1e-3,0.15*1e-3,0.15*1e-3,0.15*1e-3] #Thickness of each ply from top down
#Unit test len(stack)=len(thickness)?

#Input forces N and moments M per unit length (N/m, N)
forces=np.array([50400,1809,0]) #Nx, Ny, Nxy
moments=np.array([0,0,0]) #Mx, My, Mxy
forces.shape=[3,1]
moments.shape=[3,1]
#set force column vectors
force_moment=np.concatenate((forces,moments))

#Configure stacking sequence and create array of ply heights about the midplane
midplane=np.sum(thickness)/2 #Geometric midplane
#Array element right before the geometric midplane
#This is the reason unsymmetric laminates will break the calculator
mid=int(len(thickness)/2)
        
if (len(thickness)%2==1): #Odd number of plies
    height=np.zeros(len(thickness)+2)
    stack.insert(mid+int(1),stack[mid]) #Duplicate the middle ply angle
    thickness[mid]=thickness[mid]/2 #Halve the thickness of the middle ply
    thickness.insert(mid+1,thickness[mid]) #Duplicate the middle ply thickness
    print("Odd \n",height,", ",stack,"\n")
else: #Even number of plies
    height=np.zeros(len(thickness)+1)
    print("Even \n",height,", ",stack,"\n")

height = [np.sum(thickness[0:i])-midplane for i in range(len(height))]
print("(midplane, height array) ",midplane,", ",height,'\n')