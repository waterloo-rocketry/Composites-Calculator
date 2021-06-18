import math

import numpy as np


class Stack:
    def __init__(self, layers, force, moment):
        self.layers = layers
        self.force = force
        self.moment = moment

    def process_layers(self):
        total_height=0
        for layer in self.layers:
            total_height=total_height+layer.thickness
        self.midplane=total_height/2


        h_prev=-self.midplane
        for layer in self.layers:
            layer.set_Q_bar()
            layer.set_heights(h_prev)
            h_prev = layer.height


    def get_ABD(self):
        self.A=np.zeros(shape=(3,3))
        self.B=np.zeros(shape=(3,3))
        self.D=np.zeros(shape=(3,3))

        for layer in self.layers:
            self.A = self.A+layer.Q_bar*(layer.height-layer.h_prev)
            self.B = self.B+layer.Q_bar*(layer.height**2-layer.h_prev**2)/2
            self.D = self.D+layer.Q_bar*(layer.height**3-layer.h_prev**3)/3

        self.ABD = np.concatenate(
            (
                np.concatenate((self.A, self.B), axis=0),
                np.concatenate((np.transpose(self.B), self.D), axis=0)
            ),
            axis=1)
        return self.ABD

    def get_strains_and_stresses(self):
        force_moment = np.append(self.force, self.moment)
        # the rounding is here because of floating point strangeness 10^-10 precision was chosen with
        # an arbitrary process:
        # in testing, it was the value that made all of the values that should have been 0 by rounding
        # TODO: decide between
        #  A - continue to round to this amount, to ensure all rounding errors are removed(trig rounding, etc)
        #  B - decrease rounding, to ensure all floating point errors are handled, but not trig
        #  C - decrease rounding a lot, decide on another way of handling that sometimes 0*large_number != 0
        #  D - create special cases for 0s, by handling symmetric laminates, or evaluating trig expressions symbolically
        midstrain = np.round(np.matmul(np.linalg.inv(self.ABD), force_moment), decimals=6)

        self.estrain_glob = midstrain[:3]
        self.kstrain_glob = midstrain[3:]


        for layer in self.layers:
            layer.get_global_values(self.estrain_glob, self.kstrain_glob)
            layer.get_local_values()


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