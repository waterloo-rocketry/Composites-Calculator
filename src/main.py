from data_loading import get_layers, load_forces
from stack import Stack

if __name__ == "__main__":

    layer_data = get_layers("../data/layers.csv")
    forces_data = load_forces("../data/forces.csv")
    stack = Stack(layer_data, forces_data)
