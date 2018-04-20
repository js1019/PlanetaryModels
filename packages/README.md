Here are three packages that are used for model building.

1. distmesh by Per-Olof Persson

2. fmmlib3d from Leslie Greengard and Zydrunas Gimbutas

3. a MATLAB script for vtk users by Shawn Walker


If you need to compute the gravitational field, please compile 
fmmlib3d-1.2: 

For a linux machine, please do

module load matlab

make lib

make mwrap

make mex-matlab 

You can then call it in the MATLAB scripts. 
