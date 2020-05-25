import numpy as np
import matplotlib.pyplot as plt
from math import cos
from math import sin
from math import pi

np.set_printoptions(precision=4)

#Stiffness matrix Q
Q0=[[155.7*1e7,3.02*1e7,0*1e7],[3.02*1e7,12.16*1e7,0*1e7],[0*1e7,0*1e7,4.40*1e7]]
print (np.array(Q0))
#Angle between local and global coordinate system for each ply
theta=0

#Laminate stacking sequence and thicknesses from top down
stack=[0,pi/2,pi/2,0]
thickness=[-0.15*1e-3,-0.15*1e-3,0,0.15*1e-3,0.15*1e-3]
print (thickness)
height=np.zeros(len(thickness))
midplane=int(len(thickness)/2)

#Forces N and moments M per unit length (N/m, N)
forces=np.array([5000,0,0])
moments=np.array([0,0,0])
forces.shape=[3,1]
moments.shape=[3,1]

#Create array storing the heights of k plies about the midplane
for i in range(len(thickness)):
    if (i<midplane):
        height[i]=np.sum(thickness[i:midplane])
    else:
        height[i]=np.sum(thickness[midplane:i+1])

#Transformation matrix T1 inverse and T2 (local <--> global coordinate system)
#Function returns transformed stiffness matrix Qbar
def transform(Q0,ply_angle):
    m=cos(ply_angle)
    n=sin(ply_angle)
    T1=[[m**2,n**2,2*m*n],[n**2,m**2,-2*m*n],[-m*n,m*n,m**2-n**2]]
    T2=[[m**2,n**2,m*n],[n**2,m**2,-m*n],[-2*m*n,2*m*n,m**2-n**2]]
    T1_inv=np.linalg.inv(T1)
    Qbar=np.linalg.multi_dot([T1_inv,Q0,T2])
   
    return Qbar

#Calculate A, B, D matrices by integrating transformed stiffness matrix over all plies based on laminate 
#constitutive equations
#Units for A, B, D matrices respectively are hN/m, hN, hNm

A=np.zeros([3,3])
B=np.zeros([3,3])
B_inv=B
D=np.zeros([3,3])

for i in range(len(Q0)):
    for j in range(len(Q0[0])):
        for k in range(len(stack)):
            Qbar=transform(Q0,stack[k])
            ply_kA=height[k+1]-height[k]
            ply_kB=(height[k+1]**2-height[k]**2)/2
            ply_kD=(height[k+1]**3-height[k]**3)/3
            A[i,j]+=ply_kA*Qbar[i,j]
            B[i,j]+=ply_kB*Qbar[i,j]
            D[i,j]+=ply_kD*Qbar[i,j]

print("\n","A \n",A,"\n B \n",B,"\n D \n",D,"\n")            
A_inv=np.linalg.inv(A)

#Set ABD matrix, and force column vectors
ABD=np.concatenate((np.concatenate((A,B),axis=0),np.concatenate((B,D),axis=0)),axis=1)
print ("\n ABD stiffness matrix in hectanewton/meter, hectanewton, and hectanewton meter respectively: \n",ABD,"\n")
force_moment=np.concatenate((forces,moments))

#Calculate global laminate strains and radii of curvature at midplane
midstrain=np.dot(np.linalg.inv(ABD),force_moment)
print ("\n Strains (ex, ey, exy) and radii of curvature (kx, ky, kxy): \n",midstrain)

#Calculate strain in each ply k and plot
estrain=midstrain[:3]
kstrain=midstrain[3:6]
estrain.shape=[3,1]
kstrain.shape=[3,1]

plystrain=np.zeros(shape=(3,len(height)))
plystrain1=np.zeros(shape=(3,len(height)))

#Sanity check to ensure array dimensions are correct
print ("plystrain shape",plystrain.shape)
print ("height shape",height.shape)
print("estrain shape",estrain.shape)
print("kstarin shape",kstrain.shape)

for i in range(3):
    for j in range(len(height)):
        plystrain[i,j]=estrain[i,0]+np.dot(height[j],kstrain[i,0])
print("\n Strain for each ply k: \n",plystrain,"\n")
        
for i in range(3):
    plt.plot(plystrain[i,:5],height)
    plt.gca().invert_yaxis()
    plt.xlabel("Laminate Global Strains x10^-2")
    plt.ylabel("Height Above and Below Midplane (mm)")
    plt.grid(axis="both")
    plt.title("Strain Through Laminate Thickness ")
    
print ("forces \n",forces,"\n moments \n",moments,"\n strain \n",midstrain)
