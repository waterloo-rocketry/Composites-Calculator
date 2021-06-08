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

    def tsai_wu(self, plystress1, plystress2, plystress12):
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