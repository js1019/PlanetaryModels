## Demos for model building
The demos under CONST/ and PREM/ contain scripts to create radial ball models. You may directly run these demos. 

### Build similar models with various degrees of freedom 
Prior to build a planetary model, we utilize a radial model that is stored in ../../radialmodels/prem3L_noocean.mat. 
I will use CONST/CONST_mesh.m, you will need to set the path to you mesh and model with its name 
~~~
fmesh  = './output/CONST3k/CONST_1L_3k'; % note that CONST_1L_3k is the model name. 
~~~
You can then set up the finite element polynomial order 
~~~
pOrder  = 1; %(choose 1 or 2)
~~~
You then set up the value to contral the degrees of freedom 
~~~
a = 5e8; % you may play with this value to see the changes. 
~~~
Note that a here denotes the maximum volumn of any elements in km^3. The value _a_ controls the total number of the elements. 
You can then load a unit surface mesh 
~~~
load ../../unitspheres/data/Sph392.mat; % you have other options, see ../../unitspheres/data/. 
~~~
i.e., you may use 
~~~
load ../../unitspheres/data/Sph3k6.mat; or Sph6k, Sph15k, etc.  
~~~
The surface mesh provides basic resolution of the surface. You can then simply run 
~~~
run CONST_mesh
~~~
to obtain both meshes and models. 

### Computing the reference gravity 
Once we obtain the mesh and model, we compute the reference gravity, if needed, using Fast Multiple Method. 
Similarly, you can set up the path to the model and the finite element polynomial order 
~~~
fmesh  = './output/CONST3k/CONST_1L_3k'; % note that CONST_1L_3k is the model name. 
pOrder  = 1; %(choose 1 or 2)
~~~
The scaling factor is used for the vtk files and visualization only 
~~~
scaling = 6.371*10^3; % here we use the radial of the Earth. 
~~~
It will not affect the computation. 


