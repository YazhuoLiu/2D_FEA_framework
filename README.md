# 2D_FEA_framework
a FEA framework for two dimensional problems based on Matlab.
Code style: Object-oriented
Support for parallel computing: Yes

I have done:
Linear isotropic elastic material element creation.
T3 T6 Q4 Q8 element object creation.
Mesh object creation:
  generate computational space by code (generate simple geometry).
  Triangular meshes can be generated automatically.
  Quadrilateral meshs can only be import from txt file.
function to assemble stiffness matrices.
function to compute equivalent nodal force.
function to compute stress and strain by deformation.
function to achieve interpolation and data visualization.
function to write report.

I will develop:
more complex material model.
Complete the algorithm of drawing quadrilateral mesh.
A variety of triangular elements with different order.
function to generate complex geometry.
Implement this framework in multiple languages: Python/C++/Fortran etc.



