## Visualization of your normal modes 
Here we provide two scripts to visualize the normal modes with two simple demos. 
If the model does _not_ contain any fluid, you may use **visualCmain.m**; if the model contains fluid regions, you may use **visualEmain.m**.
The main difference is whether you obtain a file named name_vstat.dat. Please see demos for two output example. 
The name_vlist.dat contains the node ordering and name_vstat.dat contains node status, i.e., whether it is a fluid node, or a solid node, or a fluid-solid boundary node.  

You may directly run these two scripts to obtain the vtk files of your normal modes. 
If you want to plot other modes, please make sure that you set up the correct paths of different files. 
The setting provides information about the needed modes 
~~~
fmesh  = '../Demos/CONST/output/CONST3k/';
fout   = 'demos/CONST3k/';
fbase  = 'CONST_1L_3k.1';
fdtail = '0.00000000_1.00000000';
JOB = 1; pOrder = 1; nproc = 1; nth = 7; 
~~~
We note that JOB, pOrder are the inputs while you run this application; 
nproc denotes the number of processes; nth denotes the n-th modes of the mode collection. 

It is not difficult to visualize modes using ParaView's vector plots. 
We recommend readers to follow [this tutorial](http://www.bu.edu/tech/support/research/training-consulting/online-tutorials/paraview/#VECTOR). Please check "Vector Visualization Algorithms". 

Tips: use "Glyph": play "Scaling" options: "Scale Mode" -> "vector", 
Reset "Scale Factor" (to let the arrow sizes represent the vector magnitudes); 
play "Masking" options: set a different "Maximum Number Of Sample Points".
