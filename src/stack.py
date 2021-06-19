import numpy as np
from layer import Layer


class Stack:
    def __init__(self, stack_data):
        forces_data = stack_data['forces']
        layers_data = stack_data['layers']
        self.force = np.array([forces_data['Fx'],forces_data['Fy'],forces_data['Fz'],forces_data['Mx'],forces_data['My'], forces_data['Mxy']])
        self.layers = []

        total_height = 0
        for layer in layers_data:
            self.layers.append(Layer(float(layer['angle'])*np.pi, float(layer['thickness']), layer['material'], total_height))
            total_height = total_height + float(layer['thickness'])
        self.midplane = total_height / 2

        self.ABD = self.get_ABD()

        midstrain = self.get_strains_and_stresses()

        self.estrain_glob = midstrain[:3]
        self.kstrain_glob = midstrain[3:]
        self.failed, self.failed_layers_indices = self.failure_criterion()


    def get_ABD(self):
        A = np.zeros(shape=(3, 3))
        B = np.zeros(shape=(3, 3))
        D = np.zeros(shape=(3, 3))

        for layer in self.layers:
            h_0 = layer.height - self.midplane
            h_1 = layer.height + layer.thickness - self.midplane
            A = A + layer.Q_bar * (h_1 - h_0)
            B = B + layer.Q_bar * (h_1 ** 2 - h_0 ** 2) / 2
            D = D + layer.Q_bar * (h_1 ** 3 - h_0 ** 3) / 3

        ABD = np.concatenate(
            (
                np.concatenate((A, B), axis=0),
                np.concatenate((np.transpose(B), D), axis=0)
            ),
            axis=1)
        return ABD

    def get_strains_and_stresses(self):

        # the rounding is here because of floating point strangeness 10^-10 precision was chosen with
        # an arbitrary process:
        # in testing, it was the value that made all of the values that should have been 0 by rounding
        # TODO: decide between
        #  A - continue to round to this amount, to ensure all rounding errors are removed(trig rounding, etc)
        #  B - decrease rounding, to ensure all floating point errors are handled, but not trig
        #  C - decrease rounding a lot, decide on another way of handling that sometimes 0*large_number != 0
        #  D - create special cases for 0s, by handling symmetric laminates, or evaluating trig expressions symbolically
        midstrain = np.round(np.matmul(np.linalg.inv(self.ABD), self.force.astype(float)), decimals=6)

        for layer in self.layers:
            layer.get_ply_stress_strain(midstrain[:3], midstrain[3:], self.midplane)

        return midstrain

    def failure_criterion(self):
        failed = False
        ply_failure_indices = []

        for i in range(len(self.layers)):
            tsai_wu_value = self.layers[i].tsai_wu()

            if tsai_wu_value >= 1:
                ply_failure_indices.append(i)
                failed = True

        return (
            failed,
            ply_failure_indices,
        )
