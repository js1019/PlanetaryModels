This repository provides MATLAB scripts for generating a planet model on a tetrahedal mesh,
as well as it reference gravity field. 


If you only need to build a tetrahedal mesh, please compile tetgen. 

For a linux machine, just go to ./tetgen1.5.0 and type make. 

If you need to compute the reference gravity field, please go to ./package/stfmmlib3d-1.2, and compile it. 

If you use a linux node, you may follow the README.md under ./package


This repository uses four different packages: 

1. TetGen by Hang Si

2. distmesh by Per-Olof Persson

3. fmmlib3d from Leslie Greengard and Zydrunas Gimbutas

4. a MATLAB script for vtk users by Shawn Walker

I put most of the original materials under different folders 
with some minor changes for this application.
