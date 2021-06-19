#Composites Calculator - WRT
##File format
###Note: This should probably be JSON instead
###Forces file
The file at ./data/forces.csv has one line with the format:

Fx(N),Fy(N),Fz(N),Mx(Nm),My(Nm),Mxy(Nm)

###Layers file
The file at ./data/layers.csv has a line for each layer in the stack with the format:

angle(pi radians),thickness(m),path_to_material

Where path_to_material is the path to the materials file that applies to that layer
###Materials file
The file at ./data/material.csv has the format:

F1t(N),F1c(N),F2t(N),F2c(N),F12(N),load/calc,path

ehere load/calc is a choice of whether to load Q or calculate, path is the path to either the Q file or the material props
###Q File
The file at ./data/forces.csv has the format:

Fx(N),Fy(N),Fz(N),Mx(Nm),My(Nm),Mxy(Nm)
###Properties file
The file at ./data/forces.csv has the format:

Fx(N),Fy(N),Fz(N),Mx(Nm),My(Nm),Mxy(Nm)
##Todo:
- Get the tests finished and working for all three, find one for failure criterion
- Create a pretty and useful output
- Convert to use JSON, ask for the single file, add explanation of format
