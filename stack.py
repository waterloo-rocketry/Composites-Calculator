import math

import numpy as np


class Stack:
    def __init__(self, midplane, layers, force, moment):
        self.midplane = midplane
        self.layers = layers
        self.force = force
        self.moment = moment

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