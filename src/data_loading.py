import csv

import numpy as np

from layer import Layer

def load_Q_0():
    # format of the q matrix file is:
    # Q11,Q12,Q13
    # Q21,Q22,Q23
    # Q31,Q32,Q33
    Q = np.genfromtxt('../data/q.csv', delimiter=',')
    return Q

def calculate_Q_0():
    # format of the material props file is
    # E1,E2,G12,v12
    with open('../data/material_properties.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            [E1, E2, G12, v12] = row
    E1 = float(E1)
    E2 = float(E2)
    G12 = float(G12)
    v12 = float(v12)

    Q11 = E1 / (1 - v12 ** 2)
    Q12 = v12 * E1 / (1 - v12 ** 2)
    Q22 = E2 / (1 - v12 ** 2)
    Q66 = G12
    Q = np.array([[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]])
    return Q

def get_layers(Q_0, layers_file):
    # format of the layers file is:
    # angle1(pi radians),thickness1(m), F1t1, F1c1, F2t1, F2c1, F121
    # ...
    # angleN(pi radians),thicknessN(m), F1tN, F1cN, F2tN, F2cN, F12N
    layers = []
    with open(layers_file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            layers.append(
                [float(row[0]), float(row[1]),{ 'F1t':float(row[2]), 'F1c':float(row[3]), 'F2t':float(row[4]), 'F2c':float(row[5]),
                  'F12':float(row[6])}, Q_0])

    return layers




def load_forces(forces_file):
    # format of the forces file is:
    # forceX(N),forceY(N),forceZ(N)
    # momentX(N/m),momentY(N/m),momentZ(N/m)
    with open(forces_file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        force = np.array(reader.__next__())
        moment = np.array(reader.__next__())
    return force, moment