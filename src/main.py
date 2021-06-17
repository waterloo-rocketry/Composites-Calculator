import argparse

from data_loading import get_layers, load_forces, load_Q_0, calculate_Q_0
from stack import Stack

def check_flags():
    # handle a command line flag
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", action="store_true")
    all_flag_values = parser.parse_args()
    return all_flag_values.q

if __name__ == "__main__":
    should_load_q = check_flags()

    if should_load_q:
        Q_0 = load_Q_0()
    else:
        Q_0 = calculate_Q_0()

    layers = get_layers(Q_0)
    forces, moments = load_forces()

    stack = Stack(layers, forces, moments)
    stack.process_layers()
    stack.get_ABD()
    stack.get_global()
    stack.convert_to_local()

    stack.failure_criterion()
