import json
from sys import argv

from stack import Stack

if __name__ == "__main__":
    if(len(argv) == 1):
        print('Please specify the file path for the stack')

    stack_path = argv[1]

    with open(stack_path) as f:
        stack_data = json.load(f)

    stack = Stack(stack_data)
    stack.produce_text_output()
