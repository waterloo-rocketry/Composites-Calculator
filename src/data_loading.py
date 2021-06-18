import csv

import numpy as np

def get_material(materials_file):
    with open(materials_file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        value = next(reader)
    return value

def get_layers(layers_file):
    layers = []
    with open(layers_file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            layers.append(row)
    return layers

def load_forces(forces_file):
    with open(forces_file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        force = np.array(reader.__next__())
    return force