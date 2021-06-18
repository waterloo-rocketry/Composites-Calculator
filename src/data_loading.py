import csv
import numpy as np


class DataHelper:

    @classmethod
    def get_material(cls, materials_file):
        with open(materials_file) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            value = next(reader)
        return float(value[0]), float(value[1]), float(value[2]), float(value[3]), float(value[4]), value[5], value[6]

    @classmethod
    def get_layers(cls, layers_file):
        layers = []
        with open(layers_file) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                layers.append((float(row[0])*np.pi, float(row[1]), row[2]))
        return layers

    @classmethod
    def load_forces(cls, forces_file):
        with open(forces_file) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            force = np.array(reader.__next__())
        return force