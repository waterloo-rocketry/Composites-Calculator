import numpy as np


class Layer:
    def __init__(self, angle, thickness, F1t, F1c, F2t, F2c, F12):
        self.angle = angle
        self.thickness = thickness
        self.F1t = F1t
        self.F1c = F1c
        self.F2t = F2t
        self.F2c = F2c
        self.F12 = F12
        self.global_ply_strain=[[],[]]
        self.global_ply_stress=[[],[]]


    def get_global_values(self, estrain, kstrain):
        self.global_ply_strain = np.zeros(shape=(2,3))
        self.global_ply_stress = np.zeros(shape=(2,3))
        for i in range(3):
            self.global_ply_strain[0][i]=estrain[i] + self.height * kstrain[i]

        self.global_ply_stress[0] = np.dot(self.Q_bar, self.global_ply_strain[0, :])
        self.global_ply_stress[1] = np.dot(self.Q_bar, self.global_ply_strain[1, :])

    def get_local_values(self):
        self.local_ply_strain = np.zeros(shape=(2, 3))
        self.local_ply_stress = np.zeros(shape=(2, 3))

        self.local_ply_strain[0] = np.dot(self.global_ply_strain[0], self.T1)
        self.local_ply_strain[1] = np.dot(self.global_ply_strain[1], self.T1)
        self.local_ply_stress[0] = np.dot(self.global_ply_stress[0], self.T1)
        self.local_ply_stress[1] = np.dot(self.global_ply_stress[1], self.T1)

    def set_Q_bar(self, Q_0):
        m = np.cos(self.angle)
        n = np.sin(self.angle)
        self.T1 = np.array(
            [[m ** 2, n ** 2, 2 * m * n], [n ** 2, m ** 2, -2 * m * n], [-m * n, m * n, m ** 2 - n ** 2]])
        self.T2 = np.array(
            [[m ** 2, n ** 2, m * n], [n ** 2, m ** 2, -m * n], [-2 * m * n, 2 * m * n, m ** 2 - n ** 2]])
        Q_bar = np.matmul(np.matmul(np.linalg.inv(self.T1), Q_0), self.T2)
        self.Q_bar = Q_bar

    def set_heights(self, h_prev):
        self.h_prev = h_prev
        self.height = h_prev + self.thickness

    def tsai_wu(self):
        plystress1 = self.local_ply_stress[0][0]
        plystress2 = self.local_ply_stress[0][1]
        plystress12 = self.local_ply_stress[0][2]
        f1 = 1 / self.F1t - 1 / self.F1c
        f11 = 1 / (self.F1t * self.F1c)
        f2 = 1 / self.F2t - 1 / self.F2c
        f22 = 1 / (self.F2t * self.F2c)
        f12 = -1 / 2 * (f11 * f22) ** (0.5)
        f66 = 1 / (self.F12 ** 2)

        return (f1 * plystress1 +
                f2 * plystress2 +
                f11 * plystress1 * plystress1 +
                f22 * plystress2 * plystress2 +
                f66 * plystress12 * plystress12 +
                2 * f12 * plystress1 * plystress2)