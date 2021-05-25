import numpy as np
from Q_and_Qbar import transform_Q

# Calculate A, B, D matrices by integrating transformed stiffness matrix over all plies based on laminate
# constitutive equations
# Units for A, B, D matrices respectively are hN/m, hN, hNm


def ABD_triple(Q0, stack, height):
    A = np.zeros([3, 3])
    B = np.zeros([3, 3])
    D = np.zeros([3, 3])

    for i in range(len(Q0)):
        for j in range(len(Q0[0])):
            for k in range(len(stack)):
                Qbar = transform_Q(Q0, stack[k])
                ply_kA = height[k+1]-height[k]
                ply_kB = (height[k+1]**2-height[k]**2)/2
                ply_kD = (height[k+1]**3-height[k]**3)/3
                A[i, j] += ply_kA*Qbar[i, j]
                B[i, j] += ply_kB*Qbar[i, j]
                D[i, j] += ply_kD*Qbar[i, j]

    return (A, B, D)

# Set ABD matrix
# I think I made a mistake here. The extension-twisting coupling matrix should be B transpose


def ABD_matrix(A, B, D):
    return np.concatenate(
        (
            np.concatenate((A, B), axis=0),
            np.concatenate((B, D), axis=0)
        ),
        axis=1
    )
