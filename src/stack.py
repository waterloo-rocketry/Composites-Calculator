import math

import numpy as np


class Stack:
    def __init__(self, layers, force, moment):
        self.layers = layers
        self.force = force
        self.moment = moment

    def process_layers(self):
        total_height=0
        h_prev=0
        for layer in self.layers:
            layer.set_Q_bar()
            layer.set_heights(h_prev)
            h_prev = layer.height
            total_height=total_height+layer.thickness

        self.midplane=total_height/2


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
                A[i][j] = A_temp
                B[i][j] = 1 / 2 * B_temp
                D[i][j] = 1 / 3 * D_temp
        self.ABD = np.concatenate(
            (
                np.concatenate((A, B), axis=0),
                np.concatenate((np.transpose(B), D), axis=0)
            ),
            axis=1)
        return self.ABD

    def get_global(self):
        force_moment = np.append(self.force.astype(float), self.moment.astype(float))
        midstrain = np.dot(np.linalg.inv(self.ABD), force_moment)
        estrain_glob = midstrain[:3]
        print(print(midstrain))
        kstrain_glob = midstrain[3:6]


        for k in range(len(self.layers)):
            self.layers[k].get_global_values(estrain_glob, kstrain_glob);

    def convert_to_local(self):
        for k in range(len(self.layers)):
            self.layers[k].get_local_values()

    def failure_criterion(self):
        has_not_failed = True
        ply_failure_indices = []
        ply_failure_reciprocal_tsai_wu_values = []

        for i in range(len(self.layers)):
            tsai_wu_value = self.layers[i].tsai_wu()

            if (tsai_wu_value >= 1):
                ply_failure_indices.append(i)
                has_not_failed = False
                ply_failure_reciprocal_tsai_wu_values.append(1 / tsai_wu_value)

        return (
            has_not_failed,
            ply_failure_indices,
            ply_failure_reciprocal_tsai_wu_values
        )