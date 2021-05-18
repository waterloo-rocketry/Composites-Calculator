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
                    B_temp = B_temp + self.layers[k].Q_bar[i][j] * (self.layers[k].height ** 2 - self.layers[k].h_prev ** 2)
                    D_temp = D_temp + self.layers[k].Q_bar[i][j] * (self.layers[k].height ** 3 - self.layers[k].h_prev ** 3)
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
        # angle1(pi radians),thickness1(m)
        # ...
        # angleN(pi radians),thicknessN(m)
        layers = []
        height = 0
        with open('./data/layers.csv') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                height = height + float(row[1])
                layers.append(Layer(float(row[0]), float(row[1])))
        midplane = height / 2

        h_prev = 0
        for layer in layers:
            layer.set_Q_bar(Q_0)
            layer.set_heights(h_prev)
            h_prev = layer.height

        return (midplane, layers)


class Layer:
    def __init__(self, angle, thickness):
        self.angle = angle
        self.thickness = thickness

    def set_Q_bar(self, Q_0):
        m = np.cos(self.angle)
        n = np.sin(self.angle)
        T1 = np.array([[m ** 2, n ** 2, 2 * m * n], [n ** 2, m ** 2, -2 * m * n], [-m * n, m * n, m ** 2 - n ** 2]])
        T2 = np.array([[m ** 2, n ** 2, m * n], [n ** 2, m ** 2, -m * n], [-2 * m * n, 2 * m * n, m ** 2 - n ** 2]])
        Q_bar = np.matmul(np.matmul(np.linalg.inv(T1), Q_0), T2)
        self.Q_bar = Q_bar

    def set_heights(self, h_prev):
        self.h_prev=h_prev
        self.height = h_prev + self.thickness



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



