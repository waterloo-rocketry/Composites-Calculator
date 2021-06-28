import json

import numpy as np

from material import Material


class Layer:
    """Represents a layer of laminate material in a stack"""

    def __init__(self, angle, thickness, file, height):
        """
        The layer takes:
        an angle (pi radians,
        a thickness (meters),
        a path to the materials file,
        height in the stack (meters), with the height measured from the bottom of the stack to the bottom of the layer
        """
        self.angle = angle
        self.thickness = thickness
        with open(file) as f:
            material_data = json.load(f)
        self.material = Material(material_data)

        m = np.cos(self.angle)
        n = np.sin(self.angle)
        self.T1 = np.array(
            [[m ** 2, n ** 2, 2 * m * n], [n ** 2, m ** 2, -2 * m * n], [-m * n, m * n, m ** 2 - n ** 2]])
        self.T2 = np.array(
            [[m ** 2, n ** 2, m * n], [n ** 2, m ** 2, -m * n], [-2 * m * n, 2 * m * n, m ** 2 - n ** 2]])
        Q_bar = np.matmul(np.matmul(np.linalg.inv(
            self.T1), self.material.Q_0), self.T2)
        self.Q_bar = Q_bar

        self.height = height

        self.global_ply_strain = []
        self.global_ply_stress = []
        self.local_ply_strain = []
        self.local_ply_stress = []

    def get_ply_stress_strain(self, estrain, kstrain, midplane):
        """Evaluates the strains and stresses in a global and a local coordinate frame"""
        self.global_ply_strain = estrain + (self.height - midplane) * kstrain
        self.global_ply_stress = np.dot(self.Q_bar, self.global_ply_strain)
        self.local_ply_strain = np.dot(self.T1, self.global_ply_strain)
        self.local_ply_stress = np.dot(self.T2, self.global_ply_stress)

    def tsai_wu(self):
        """Applies the tsai wu failure criterion to the layer
          must be done after the stresses and strains have been evaluated"""

        f1 = 1 / self.material.F1t - 1 / self.material.F1c
        f11 = 1 / (self.material.F1t * self.material.F1c)
        f2 = 1 / self.material.F2t - 1 / self.material.F2c
        f22 = 1 / (self.material.F2t * self.material.F2c)
        f12 = -1 / 2 * (f11 * f22) ** 0.5
        f66 = 1 / (self.material.F12 ** 2)

        tsai_wu_criterion = (f1 * self.local_ply_stress[0] +
                             f2 * self.local_ply_stress[1] +
                             f11 * self.local_ply_stress[0] ** 2 +
                             f22 * self.local_ply_stress[1] ** 2 +
                             f66 * self.local_ply_stress[2] ** 2 +
                             2 * f12 * self.local_ply_stress[0] * self.local_ply_stress[1])

        return tsai_wu_criterion
