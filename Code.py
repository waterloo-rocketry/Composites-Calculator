#Only works for symmetric ply thicknesses for now
#Reference: Example 1 page 17 Chapter 9 Lecture Slides for strain calculation
#Reference: Example page 19 Chapter 8 Lecture Slides for Primary and Principle stress plots
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi
from Q_Qbar_and_constants import *
from ABD_matricies import *
from stresses_strains import *

np.set_printoptions(precision=4)

#Transformation matrix T1 inverse (primary strain --> principal strain)
#Transformation T2 (principal stress --> primary stress)

#Function returns the strain in primary material directions as a column vector
def transform_T2(principalStrain,ply_angle):
    m=cos(ply_angle)
    n=sin(ply_angle)
    T2=[[m**2,n**2,m*n],[n**2,m**2,-m*n],[-2*m*n,2*m*n,m**2-n**2]]
    primaryStrain=np.dot(T2,principalStrain)
    
    return primaryStrain

A,B,D = ABD_tripple(Q0,stack,height)

print("\n","A \n",A,"\n B \n",B,"\n D \n",D,"\n")            
A_inv=np.linalg.inv(A)

ABD = ABD_matrix(A,B,D)
print ("\n ABD stiffness matrix in hectanewton/meter, hectanewton, and hectanewton meter respectively: \n",ABD,"\n")

#Calculate principal laminate strains and radii of curvature at midplane
midstrain=np.dot(np.linalg.inv(ABD),force_moment)
print ("\n Strains (ex, ey, exy) and radii of curvature (kx, ky, kxy): \n",midstrain)

#Assign variables for principal strains and stresses in each ply k
#Calculate principal and primary strains and stresses in each ply k
#Strain varies linearly
estrain_glob=midstrain[:3]
kstrain_glob=midstrain[3:6]

plystrain_glob, plystress_glob = global_ply(Q0,stack,estrain_glob,kstrain_glob,height)
plystrain_loc, plystress_loc = local_ply(plystrain_glob, plystress_glob, height, stack)

print("\n Principal Strain for each ply boundary: \n",plystrain_glob,"\n")
print("\n Principal Stress for each ply boundary: \n",plystress_glob,"\n")
print("\n Primary Strain for each ply boundary: \n",plystrain_loc,"\n")
print("\n Primary Stress for each ply boundary: \n",plystress_loc,"\n")
    
plot_height=np.zeros(2*len(height)-2) #x values for plotting purposes mwahahahah >:]  
for k in range(len(height)-1):
    plot_height[2*k]=height[k]
    plot_height[2*k+1]=height[k+1]
