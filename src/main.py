from data_loading import DataHelper
from stack import Stack

if __name__ == "__main__":

    layer_data = DataHelper.get_layers("../data/layers.csv")
    forces_data = DataHelper.load_forces("../data/forces.csv")
    stack = Stack(layer_data, forces_data)
