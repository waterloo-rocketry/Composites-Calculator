import numpy as np

# Laminate parameters
stack = [0, np.pi/2, np.pi/2, 0]  # Stacking sequence from top down
# Thickness of each ply from top down
thickness = [0.15*1e-3, 0.15*1e-3, 0.15*1e-3, 0.15*1e-3]
ply_num = len(thickness)

# Input forces N and moments M per unit length (N/m, N)
forces = np.array([50400, 1809, 0])  # Nx, Ny, Nxy
moments = np.array([0, 0, 0])  # Mx, My, Mxy
# set force column vectors
force_moment = np.concatenate((forces, moments))

# Configure stacking sequence and create array of ply heights about the midplane
midplane = np.sum(thickness)/2  # Geometric midplane
# Array element right before the geometric midplane
# This is the reason unsymmetric laminates will break the calculator
mid = int(len(thickness)/2)

if (len(thickness) % 2 == 1):  # Odd number of plies
    height = np.zeros(len(thickness)+2)
    stack.insert(mid+int(1), stack[mid])  # Duplicate the middle ply angle
    thickness[mid] = thickness[mid]/2  # Halve the thickness of the middle ply
    # Duplicate the middle ply thickness
    thickness.insert(mid+1, thickness[mid])
else:  # Even number of plies
    height = np.zeros(len(thickness)+1)

height = [np.sum(thickness[0:i])-midplane for i in range(len(height))]
