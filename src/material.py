import csv

import numpy as np

class Material:
    def __init__(self, material_data):

        self.F1t = float(material_data['F1t'])
        self.F1c = float(material_data['F1c'])
        self.F2t = float(material_data['F2t'])
        self.F2c = float(material_data['F2c'])
        self.F12 = float(material_data['F12'])

        if(material_data['Q']):
            self.Q_0 = np.array(material_data['Q'])
        elif(material_data['properties']):
            self.Q_0 = self.calculate_Q_0(material_data['properties'])


    def calculate_Q_0(self, properties):

        E1 = float(properties['E1'])
        E2 = float(properties['E2'])
        G12 = float(properties['G12'])
        v12 = float(properties['v12'])

        Q11 = E1 / (1 - v12 ** 2)
        Q12 = v12 * E1 / (1 - v12 ** 2)
        Q22 = E2 / (1 - v12 ** 2)
        Q66 = G12
        Q = np.array([[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]])
        return Q




