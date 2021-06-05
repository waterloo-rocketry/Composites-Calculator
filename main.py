import argparse
import csv
import math
import numpy as np


class Stack:
    def __init__(self):
        self.Q_0 = Stack.get_Q_0(check_flags())

        self.midplane, self.layers = Stack.get_layers(self.Q_0)

        # format of the forces file is:
        # forceX(N),forceY(N),forceZ(N)
        # momentX(N/m),momentY(N/m),momentZ(N/m)
        with open('./data/forces.csv') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            self.force = np.array(reader.__next__())
            self.moment = np.array(reader.__next__())

    def get_ABD(self):
        A = np.zeros([3, 3])
        B = np.zeros([3, 3])
        D = np.zeros([3, 3])

        for i in range(3):
            for j in range(3):
                A_temp = 0
                B_temp = 0
                D_temp = 0
                for k in range(len(self.layers)):
                    A_temp = A_temp + self.layers[k].Q_bar[i][j] * (self.layers[k].height - self.layers[k].h_prev)
                    B_temp = B_temp + self.layers[k].Q_bar[i][j] * (
                            self.layers[k].height ** 2 - self.layers[k].h_prev ** 2)
                    D_temp = D_temp + self.layers[k].Q_bar[i][j] * (
                            self.layers[k].height ** 3 - self.layers[k].h_prev ** 3)
                print(A_temp)
                A[i][j] = A_temp
                B[i][j] = 1 / 2 * B_temp
                D[i][j] = 1 / 3 * D_temp

        return np.concatenate(
            (
                np.concatenate((A, B), axis=0),
                np.concatenate((np.transpose(B), D), axis=0)
            ),
            axis=1)

    def get_global(self, ABD):
        force_moment = np.append(self.force.astype(float), self.moment.astype(float))
        midstrain = np.dot(np.linalg.inv(ABD), force_moment)
        estrain_glob = midstrain[:3]
        kstrain_glob = midstrain[3:6]
        global_ply_strain = np.zeros(shape=(3, 2 * len(self.layers)))
        global_ply_stress = np.zeros(shape=(3, 2 * len(self.layers)))

        for k in range(len(self.layers)):
            for i in range(3):
                global_ply_strain[i, 2 * k] = estrain_glob[i] + \
                                              self.layers[math.ceil(k / 2)].height * kstrain_glob[i]
            global_ply_stress[:, 2 * k] = np.dot(self.layers[k].Q_bar, global_ply_strain[:, 2 * k])
            global_ply_stress[:, 2 * k + 1] = np.dot(self.layers[k].Q_bar, global_ply_strain[:, 2 * k + 1])

        return (global_ply_stress, global_ply_strain)

    @classmethod
    def get_Q_0(cls, calculateQ):
        if (not calculateQ):
            # format of the q matrix file is:
            # Q11,Q12,Q13
            # Q21,Q22,Q23
            # Q31,Q32,Q33
            Q = np.genfromtxt('data/q.csv', delimiter=',')
        else:
            # format of the material props file is
            # E1,E2,G12,v12
            with open('./data/material_properties.csv') as csvfile:
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

    @classmethod
    def get_layers(cls, Q_0):
        # format of the layers file is:
        # angle1(pi radians),thickness1(m), F1t1, F1c1, F2t1, F2c1, F121
        # ...
        # angleN(pi radians),thicknessN(m), F1tN, F1cN, F2tN, F2cN, F12N
        layers = []
        height = 0
        with open('./data/layers.csv') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                height = height + float(row[1])
                layers.append(
                    Layer(float(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]),
                          float(row[6])))
        midplane = height / 2

        h_prev = 0
        for layer in layers:
            layer.set_Q_bar(Q_0)
            layer.set_heights(h_prev)
            h_prev = layer.height

        return (midplane, layers)

    def convert_to_local(self, global_ply_stress, global_ply_strain):
        local_ply_stress = np.zeros(shape=(3, 2 * len(self.layers)))
        local_ply_strain = np.zeros(shape=(3, 2 * len(self.layers)))

        for k in range(len(self.layers)):
            layer = self.layers[k]
            local_ply_strain[:, 2 * k] = np.dot(global_ply_strain[:, 2 * k], layer.T1)
            local_ply_strain[:, 2 * k + 1] = np.dot(global_ply_strain[:, 2 * k + 1], layer.T1)
            local_ply_stress[:, 2 * k] = np.dot(global_ply_stress[:, 2 * k], layer.T1)
            local_ply_stress[:, 2 * k + 1] = np.dot(global_ply_stress[:, 2 * k + 1], layer.T1)
        return (local_ply_stress, local_ply_strain)

    def failure_criterion(self, plystress):
        has_not_failed = True
        ply_failure_indices = []
        ply_failure_reciprocal_tsai_wu_values = []

        for i in range(len(self.layers)):
            tsai_wu_value = self.layers[i].tsai_wu(
                plystress[0][i],
                plystress[1][i],
                plystress[2][i]
            )

            if (tsai_wu_value >= 1):
                ply_failure_indices.append(i)
                has_not_failed = False
                ply_failure_reciprocal_tsai_wu_values.append(1 / tsai_wu_value)

        return (
            has_not_failed,
            ply_failure_indices,
            ply_failure_reciprocal_tsai_wu_values
        )


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


def check_flags():
    # handle a command line flag
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", action="store_true")
    all_flag_values = parser.parse_args()
    return all_flag_values.q


if __name__ == "__main__":
    stack = Stack()
    ABD = stack.get_ABD()
    global_ply_stress, global_ply_strain = stack.get_global(ABD)
    local_ply_stress, local_ply_strain = stack.convert_to_local(global_ply_stress, global_ply_strain)
    (has_not_failed, ply_failure_indices, ply_failure_reciprocal_tsai_wu_values) = stack.failure_criterion(local_ply_stress)
