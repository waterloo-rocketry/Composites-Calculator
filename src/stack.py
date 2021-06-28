import numpy as np

from layer import Layer

"""Represents a stack consisting of any number of layers, and the forces applied to it"""


class Stack:
    """
    format of stack_data

    a dictionary with the following keys:
    name - the name of the stack
    forces - a dict with Fx,Fy,Fz,Mx,My,Mxy
    layers - an array of dicts with the keys thickness,angle,material_path

    material path is a relative path to the material file to be used for the layer
    example can be found at src/data/stacks/example1.json
    """

    def __init__(self, stack_data):
        self.name = stack_data['name']
        forces_data = stack_data['forces']
        layers_data = stack_data['layers']
        self.force = np.array([forces_data['Fx'], forces_data['Fy'], forces_data['Fz'],
                              forces_data['Mx'], forces_data['My'], forces_data['Mxy']])

        self.layers = []
        total_height = 0
        for layer in layers_data:
            self.layers.append(Layer(float(
                layer['angle'])*np.pi, float(layer['thickness']), layer['material'], total_height))
            total_height = total_height + float(layer['thickness'])
        self.midplane = total_height / 2

        self.ABD = self.get_ABD()

        midstrain = self.get_strains_and_stresses()

        self.estrain_glob = midstrain[:3]
        self.kstrain_glob = midstrain[3:]
        self.failed, self.failed_layers_indices = self.failure_criterion()

    """Evaluates the stiffness matrix (ABD matrix) for the laminate"""

    def get_ABD(self):
        A = np.zeros(shape=(3, 3))
        B = np.zeros(shape=(3, 3))
        D = np.zeros(shape=(3, 3))

        for layer in self.layers:
            h_0 = layer.height - self.midplane
            h_1 = layer.height + layer.thickness - self.midplane
            A += layer.Q_bar * (h_1 - h_0)
            B += layer.Q_bar * (h_1 ** 2 - h_0 ** 2) / 2
            D += layer.Q_bar * (h_1 ** 3 - h_0 ** 3) / 3

        ABD = np.concatenate(
            (
                np.concatenate((A, B), axis=0),
                np.concatenate((np.transpose(B), D), axis=0)
            ),
            axis=1)
        return ABD

    """
    evaluates the e and k strains based on the forces applied and the ABD matrix
    
    also iterates over the layers, having them calculate the per-ply strain and stresses,
    and converts them to local orientations
    """

    def get_strains_and_stresses(self):
        midstrain = np.round(np.matmul(np.linalg.inv(
            self.ABD), self.force.astype(float)), decimals=6)

        for layer in self.layers:
            layer.get_ply_stress_strain(
                midstrain[:3], midstrain[3:], self.midplane)

        return midstrain

    """
    iterates over the layers and evaluates each of their tsai wu failure criterion
    """

    def failure_criterion(self):
        failed = False
        ply_failure_indices = []

        for index, layer in enumerate(self.layers):
            tsai_wu_value = layer.tsai_wu()

            if tsai_wu_value >= 1:
                ply_failure_indices.append(index)
                failed = True

        return (
            failed,
            ply_failure_indices,
        )

    """"prints output to the console, should probably produce more output as new requirements are identified"""

    def produce_text_output(self):
        print(f"Displaying analysis for: {self.name}")

        print('ABD Matrix: ')
        np.set_printoptions(linewidth=np.inf)
        print(self.ABD)
        print('\n')

        if(not self.failed):
            print('Laminate did not fail!')
        else:
            print('Laminate failed\n')
            print('The failed laminates were: ')
            for i in range(len(self.failed_layers_indices)):
                print(
                    f'\t Layer {i} - Thickness: {self.layers[i].thickness} meters, Angle: {self.layers[i].angle} pi radians ')
