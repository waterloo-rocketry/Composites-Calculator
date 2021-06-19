from material import Material
from stack import Stack
import json

if __name__ == "__main__":
    with open('data/stacks/example1.json') as f:
        stack_data = json.load(f)
    stack = Stack(stack_data)
