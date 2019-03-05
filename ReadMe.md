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

# Citation

If any of the provided code does help in your research, please cite the following paper:

[1] Wei, R., Bao, F., Liu, Y., Hui, W., 2018. Robust Three-Dimensional Level-Set Method for Evolving Fronts on Complex Unstructured Meshes. Mathematical Problems in Engineering 2018, 1–15. https://doi.org/10.1155/2018/2730829

BibTex item:

    @article{weiRobustThreeDimensionalLevelSet2018,
    title = {Robust {{Three}}-{{Dimensional Level}}-{{Set Method}} for {{Evolving Fronts}} on {{Complex Unstructured Meshes}}},
    volume = {2018},
    copyright = {All rights reserved},
    issn = {1024-123X, 1563-5147},
    doi = {10.1155/2018/2730829},
    abstract = {With a purpose to evolve the surfaces of complex geometries in their normal direction at arbitrarily defined velocities, we have developed a robust level-set approach which runs on three-dimensional unstructured meshes. The approach is built on the basis of an innovative spatial discretization and corresponding gradient-estimating approach. The numerical consistency of the estimating method is mathematically proven. A correction technology is utilized to improve accuracy near sharp geometric features. Validation tests show that the proposed approach is able to accurately handle geometries containing sharp features, computation regions having irregular shapes, discontinuous speed fields, and topological changes. Results of the test problems fit well with the reference results produced by analytical or other numerical methods and converge to reference results as the meshes refine. Compared to level-set method implementations on Cartesian meshes, the proposed approach makes it easier to describe jump boundary conditions and to perform coupling simulations.},
    language = {en},
    journal = {Mathematical Problems in Engineering},
    author = {Wei, Ran and Bao, Futing and Liu, Yang and Hui, Weihua},
    month = sep,
    year = {2018},
    pages = {1-15}
    }
