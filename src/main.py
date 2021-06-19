from sys import argv
from os import path
from stack import Stack
import json

if __name__ == "__main__":
    if(len(argv)==1):
        print('Please specify the file path for the stack')

    stack_path = argv[1]
    if(not path.exists(stack_path)):
        print('Specified file does not exist')

    with open(stack_path) as f:
        stack_data = json.load(f)

    stack = Stack(stack_data)
    stack.produce_text_output()
