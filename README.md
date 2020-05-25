# Composites-Calcuator
Calculate ABD matrix and per-ply stress and strain in laminates based on Classical Laminate Plate Theory

The motivation behind implementing the calculator in Python and not Excel is that array functions in Excel are a nightmare. Python is (I think?) well-equipped to handle the matrix multiplication required in the calculation of the ABD matrices and per-ply stresses. By the same token, there's a blurred line between how much of this calculator would do better in Excel, since there's also a data entry aspect (for material properties and laminate parameters) that would (I think?) be easier if done in Excel.
