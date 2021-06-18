from stack import Stack

if __name__ == "__main__":

    # TODO: make the stack have a function to calculate all of this in a row
    # TODO: allow dynamic selection of files
    # TODO: consider ui over use of CLI, potentially allow saving and selecting of materials and stacks
    stack = Stack("../data/layers.csv", "../data/forces.csv")
    stack.process_layers()
    stack.get_ABD()
    stack.get_strains_and_stresses()

    stack.failure_criterion()
