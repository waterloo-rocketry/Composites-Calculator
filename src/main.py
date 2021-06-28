import json
import sys

from stack import Stack

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print('Please specify the file path for the stack')
        sys.exit()
    stack_path = sys.argv[1]

    with open(stack_path) as f:
        stack_data = json.load(f)

    stack = Stack(stack_data)
    stack.produce_text_output()
