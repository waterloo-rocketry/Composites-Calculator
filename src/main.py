from data_loading import get_layers, load_forces
from stack import Stack

if __name__ == "__main__":

    layer_data = get_layers("../data/layers.csv")
    forces_data = load_forces("../data/forces.csv")
    # TODO: make the stack have a function to calculate all of this in a row
    # TODO: allow dynamic selection of files
    # TODO: consider ui over use of CLI, potentially allow saving and selecting of materials and stacks
    stack = Stack(layer_data, forces_data)
