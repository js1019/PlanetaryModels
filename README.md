Planetary Model Builder 
================================================================

This repository provides MATLAB scripts to build a planet model on a tetrahedal mesh,
as well as its reference gravitational field. 

How to make it work for you? 
----------------------------------------------------------------
This repository uses several packages:  
1. TetGen by Hang Si  
2. distmesh by Per-Olof Persson  
3. several MATLAB scripts from Hesthaven & Warburton's book 'Nodal Discontinuous Galerkin Methods'  
4. fmmlib3d from Leslie Greengard and Zydrunas Gimbutas  
5. a MATLAB script for vtk users by Shawn Walker
I put most of the original materials under different folders 
with some minor changes for this application.

You may follow the README.md under ./packages to compile what you need. 
If you only need to build a tetrahedal mesh, please compile tetgen; if you are a linux user, just go to packages/tetgen1.5.0 and type make.   
If you need to compute the reference gravity, please go to packages/fmmlib3d-1.2, and compile it.  


Build your planetary models
-----------------------------------------------------------------

