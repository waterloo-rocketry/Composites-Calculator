import numpy as np


class Layer:
    def __init__(self, angle, thickness, max_forces, Q_0):
        self.angle = angle
        self.thickness = thickness
        self.max_forces = max_forces
        self.Q_0 = Q_0

        self.global_ply_strain = None
        self.global_ply_stress = None


    def get_global_values(self, estrain, kstrain):
        self.global_ply_strain = estrain + self.height*kstrain
        self.global_ply_stress = np.dot(self.Q_bar, self.global_ply_strain)

    def get_local_values(self):
        self.local_ply_strain = np.dot(self.T1, self.global_ply_strain)
        self.local_ply_stress = np.dot(self.T2, self.global_ply_stress)

    def set_Q_bar(self):
        m = np.cos(self.angle)
        n = np.sin(self.angle)
        self.T1 = np.array(
            [[m ** 2, n ** 2, 2 * m * n], [n ** 2, m ** 2, -2 * m * n], [-m * n, m * n, m ** 2 - n ** 2]])
        self.T2 = np.array(
            [[m ** 2, n ** 2, m * n], [n ** 2, m ** 2, -m * n], [-2 * m * n, 2 * m * n, m ** 2 - n ** 2]])
        Q_bar = np.matmul(np.matmul(np.linalg.inv(self.T1), self.Q_0), self.T2)
        self.Q_bar = Q_bar

    def set_heights(self, h_prev):
        self.h_prev = h_prev
        self.height = h_prev + self.thickness

    def tsai_wu(self):
        plystress1 = self.local_ply_stress[0]
        plystress2 = self.local_ply_stress[1]
        plystress12 = self.local_ply_stress[2]
        f1 = 1 / self.max_forces['F1t'] - 1 / self.max_forces['F1c']
        f11 = 1 / (self.max_forces['F1t']  * self.max_forces['F1c'] )
        f2 = 1 / self.max_forces['F2t']  - 1 / self.max_forces['F2c']
        f22 = 1 / (self.max_forces['F2t'] * self.max_forces['F2c'])
        f12 = -1 / 2 * (f11 * f22) ** (0.5)
        f66 = 1 / (self.max_forces['F12'] ** 2)

        tsai_wu_criterion = (f1 * plystress1 +
                f2 * plystress2 +
                f11 * plystress1 * plystress1 +
                f22 * plystress2 * plystress2 +
                f66 * plystress12 * plystress12 +
                2 * f12 * plystress1 * plystress2)

        return tsai_wu_criterion