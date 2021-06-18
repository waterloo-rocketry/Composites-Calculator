import csv

import numpy as np

class Material:
    def __init__(self, material_data):

        self.F1t = float(material_data[0])
        self.F1c = float(material_data[1])
        self.F2t = float(material_data[2])
        self.F2c = float(material_data[3])
        self.F12 = float(material_data[4])
        mode, self.file = material_data[5:]

        if(mode=='load'):
            self.Q_0 = self.load_Q_0()
        elif(mode=='calc'):
            self.Q_0 = self.calculate_Q_0()

    def load_Q_0(self):
        Q = np.genfromtxt(self.file, delimiter=',')
        return Q

    def calculate_Q_0(self):
        with open(self.file) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            [E1, E2, G12, v12] = next(reader)
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




