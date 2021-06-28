# Composites Calculator - WRT
To run:
```
python3 main.py [path_to_stack]
```
## Input
### Stack
The input is in the form of a json file in the data/stacks directory.

The stack file format can be found in the example1.json file in this directory. Any number of layers can be added to the layers array in this folder.

The material value should be a reference to a file in the data/materials directory.

### Materials
The materials file is in the form of json in the data/materials directory.

The maximum loads should be provided, and there are two options for evaluating the stiffness.

EITHER,

a two dimensional array, with key "Q" should be provided with the Q_0 matrix for the material, as in data/materials/carbon_epoxy.json

OR

the material properties from empirical testing can be provided, as in data/materials/example_of_calculated

## Work to Be Done
### Potential Directions
A better interface might be a consideration. Creating a UI or better CLI would allow us to use this a bit easier

An editor of some kind so that we do not have to edit the JSON files directly might also be a benefit

A way of producing comparisons, or a way for it to make recommendations as to the direction that more fibre is needed, etc could reduce even more of the analysis that has to be done manually

### Testing
More testing needs to be done to verify that the failure criterion are working. There are no examples of the failure criterion on drive, we should find examples with confirmed answers

More testing of strange configurations (those with multiple materials, uneven number of layers, etc)

### Cleanup
I expect that others looking at the code will have good suggestions for clarity. Fixes should be made to refine logic and improve readability.










