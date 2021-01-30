"""
The Tsai-Wu failure Criteria for laminates.
If the expresion given is greater than or equal to
1 for any of the layers of the laminate, the laminate will fail.
The expression takes as input the maximum stresses
The laminate can take as well as the stresses it is curretly
experiencing to determine if the laminate will fail.
As input, it takes in the maximum tensile (t)
and compressive (c) stresses the ply can take
in the 1 and 2 directions and the maximum
sheet (F12) the ply can take. In addition, it takes
the current stresses (plystress1, 2 and 12) the the
ply is under.
"""

def tsai_wu(F1t, F1c, F2t, F2c, F12, plystress1, plystress2, plystress12):
    f1 = 1/F1t - 1/F1c
    f11 = 1/(F1t * F1c)
    f2 = 1/F2t - 1/F2c
    f22 = 1/(F2t * F2c)
    f12 = -1/2 * (f11 * f22)**(0.5)
    f66 = 1/(F12 ** 2)

    return (f1*plystress1 +
            f2*plystress2 +
            f11*plystress1*plystress1 +
            f22*plystress2*plystress2 +
            f66*plystress12*plystress12 +
            2*f12*plystress1*plystress2)

"""
The following code assumes that an arrays called
plystress will be created in the following form
plystress = [
         [s1, s1, s1 ...],
         [s2, s2, s2 ...],
         [s12, s12, s12 ...]
]
where each element plystress[i] of the array is the sigma 1, 2 and 12 of
for the ply at stack[i]
maxstress = [
         [F1t, F1t, F1t ...],
         [F1c, F1c, F1c ...],
         [F2t, F2t, F2t ...],
         [F2c, F2c, F2c ...]
         [F12, F12, F12 ...]
]
Where each element of maxstress[i] of the array is the
maximum stress that the ply at stack[i]
len(plystress) = len(maxstress)
plystress = np.zeros(shape=(3,len(height)))
maxstress = np.ones(shape=(5,len(height)))
Failure Criterion. This function
loops through 2 arrays of the form
described above and applies the tsai-wu
failure criterion to each on individually.
The function outputs
(boolean, [ints], [floats])
if there are no failures, boolean is True,
otherwise, it is false. The lists contain
all of the incidices and recripricols of
the tsai-wu values for plys that have failed
"""

def failure_criterion(plystress, maxstress):
    has_not_failed = True
    ply_failure_indices = []
    ply_failure_reciprocal_tsai_wu_values = []

    for i in range(len(plystress[0])):
        tsai_wu_value = tsai_wu(
                maxstress[0][i],
                maxstress[1][i],
                maxstress[2][i],
                maxstress[3][i],
                maxstress[4][i],
                plystress[0][i],
                plystress[1][i],
                plystress[2][i]
            )

        if (tsai_wu_value >= 1):
            ply_failure_indices.append(i)
            has_not_failed = False
            ply_failure_reciprocal_tsai_wu_values.append(1/tsai_wu_value)

    return (
            has_not_failed,
            ply_failure_indices,
            ply_failure_reciprocal_tsai_wu_values
        )