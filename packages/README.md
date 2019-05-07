#  External packages
Here are four packages that are used for model building:  
### 1. [distmesh](http://persson.berkeley.edu/distmesh/) by Per-Olof Persson.   
~~~
@article{persson2004simple,
    title={A simple mesh generator in {MATLAB}},
    author={Persson, Per-Olof and Strang, Gilbert},
    journal={SIAM review},
    volume={46},
    number={2},
    pages={329--345},
    year={2004},
    publisher={SIAM}
}
~~~
### 2. [fmmlib3d](https://github.com/zgimbutas/fmmlib3d) from Leslie Greengard and Zydrunas Gimbutas.    
~~~
@misc{gimbutas2011fmmlib3d,
    title={{FMMLIB3D 1.2, FORTRAN libraries for fast multiple method in three dimensions}},
    author={Gimbutas, Z and Greengard, L},
    year={2011}
}
~~~
### 3. [a MATLAB script](https://www.mathworks.com/matlabcentral/fileexchange/58002-write-binary-vtk-file-for-tetrahedral-grid-with-scalar-and-vector-data?s_tid=prof_contriblnk) for vtk users by Shawn Walker.  
### 4. [TetGen](http://www.tetgen.org) by Hang Si.
~~~
@article{si2015tetgen,
    title={{TetGen, a Delaunay-based quality tetrahedral mesh generator}},
    author={Si, Hang},
    journal={ACM Transactions on Mathematical Software (TOMS)},
    volume={41},
    number={2},
    pages={11},
    year={2015},
    publisher={ACM}
}
~~~

We illustrate the tips for compiling these packages below. 

## How to install these packages? 
### 1. distmesh
It is a Matlab code with c++ mex and no compilation required.

### 2. fmmlib3d-1.2
This Fast Multiple Method is written in Fortran with OpenMP acceleration and other interfaces. You only need to compile this if you want to compute the reference gravity.

#### 2.1 Linux
For a Linux machine with Matlab installed, please do 
~~~
make lib; make mwrap; make mex-matlab 
~~~
check "make test-mex-matlab" and make sure that it works.
You can then call it in the MATLAB scripts.

#### 2.1 Mac
Same as above. 
If you get an error 'Actual argument contains too few elements for dummy argument', use the Fortran flag -std=legacy. E.g.
~~~
make FFLAGS='-std=legacy' lib
~~~
To compile the Matlab components (mwrap, test-mex-matlab) you need to be able to call Matlab from the command line:
~~~
export PATH=$PATH:/Applications/MATLAB_R2018b.app/bin/
~~~
You cannot 'make all' unless you have Octave installed.
If Matlab has a problem finding a *.dylib file, but you have the file somewhere on your system, a quick solution is to place a symbolic link to the file in the place where Matlab is looking, e.g.
~~~
ln -s /opt/local/lib/libgcc/libgfortran.3.dylib /usr/local/lib/libgfortran.3.dylib
~~~

#### 2.2 Other platforms
Not tested on other platforms. You may have to modify the makefiles.  

### 3. VTK scripts
It is a Matlab code and no compilation required.

### 4. TetGen
~~~
cd packages/tetgen1.5.0; make
~~~




