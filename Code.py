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

## The Tsai-Wu failure Criteria for laminates.
## If the expresion given is greater than 1 for any of the
## layers of the laminate, the laminate will fail.
## The expression takes as input the maximum stresses
## The laminate can take as well as the stresses it is curretly
## experiencing to determine if the laminate will fail.
## [Note, I'll define all the variable names once I can
## get an explanation as to what they all really mean]

def tsai_wu(F1t, F1c, F2t, F2c, F12, s1, s2, s12):
    f1 = 1/F1t - 1/F1c
    f11 = 1/(F1t * F1c)
    f2 = 1/F2t - 1/F2c
    f22 = 1/(F2t * F2c)
    f12 = -1/2 * (f11 * f22)**(0.5)
    f66 = 1/(F12 ** 2)

    return f1*s1 + f2*s2 + f11*s1*s1 + f22*s2*s2 + f66*s12*s12 + 2*f12*s1*s2

## The following code assumes that an arrays called
## plystresses will be created in the following form
## plystresses = [
##          [s1, s2, s12],
##          [s1, s2, s12],
##          [s1, s2, s12],
##          ...
## ]
## where each element plystresses[i] of the array is the sigma 1, 2 and 12 of
## for the ply at stack[i]

## maxstresses = [
##          [F1t, F1c, F2t, F2c, F12],
##          [F1t, F1c, F2t, F2c, F12],
##          [F1t, F1c, F2t, F2c, F12],
##          ...
## ]
## Where each element of maxstresses[i] of the array is the
## maximum stress that the ply at stack[i]

##len(plystresses) = len(maxstresses)

plystresses = np.zeros(shape=(len(height),3))
##maxstresses can be or can not be a numpy array, it dosen't matter
maxstresses = np.array([[1,1,1,1,1] for _ in range(len(height))])

has_not_failed = True

for i in range(len(plystresses)):
    tsai_wu_value = tsai_wu(
            maxstresses[i][0],
            maxstresses[i][1],
            maxstresses[i][2],
            maxstresses[i][3],
            maxstresses[i][4],
            plystresses[i][0],
            plystresses[i][2],
            plystresses[i][3]
        )

    if (tsai_wu_value >= 1):
        print(tsai_wu_value)
        has_not_failed = False
        print("There was a failure in layer " + str(i+1))
