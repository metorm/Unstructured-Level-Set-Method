# Unstructured Level Set Method
Robust 3D Level Set Method for Evolving Fronts on Complex Unstructured Meshes

# Introduction

This repository contains demonstration code of our robust 3D level set method on complex unstructured meshes.

The main advantage of our implementation is its ability to handle discontinuous propagation speed, sharp geometric feature, and the geometries which start from the boundary of mesh.

# Usage

0. Download the code and `cd` to the downloaded directory in `MATLAB`.
1. In `src/config.m` edit the simulation configurations such as mesh scale.
2. Run `src/GenerateBasicData.m`, which generates a series of `*.mat` files which contains the generated unstructured mesh.
3. Run any of the `src/*Evolve.m`	 files, which perform the actual level set computation for corresponding problems. The results will be stored in `src/resultData` directory.
4. Run script `drawResult.m`, which draws the result of each example. Note that each example is corresponding to a boolean switch called `draw*` in this file. You may have to set the specific switch to `true` in order to see corresponding results.
5. `Reinitialize2D.m` contains an independent example, which is not described in our article. You can run the file directly (after running src/GenerateBasicData.m) to see the output figure.
6. Other files are auxiliary functions, of which the file names are quite self-explanatory.


# Acknowledgment

+ We have used [MESH2D](https://github.com/dengwirda/mesh2d) to generate the required mesh.
+ The [distance2curve](https://ww2.mathworks.cn/matlabcentral/fileexchange/34869-distance2curve) function is written by [John D'Errico](https://ww2.mathworks.cn/matlabcentral/profile/authors/869215-john-d-errico).
+ Functions `dunion`, `ddiff`, `dcircle`, `drectangle` are borrowed from [DistMesh](http://persson.berkeley.edu/distmesh/) toolkit.
+ Function [tricontour](https://ww2.mathworks.cn/matlabcentral/fileexchange/10408-contours-for-triangular-grids) is written by [Darren Engwirda](https://ww2.mathworks.cn/matlabcentral/profile/authors/870403-darren-engwirda).
