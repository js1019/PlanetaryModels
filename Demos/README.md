## Demos for model building
The demos under CONST/, PREM/, Moon/, RT_MarsDWAK/ contain scripts to create planetary models. 
You may directly run these demos. 

Four model generators are provided:
+ CONST: a constant perfectly spherically symmetric pure solid ball has a size of the Earth;
+ PREM: a standard radial Earth model from [Dziewonski & Anderson 1981](https://www.sciencedirect.com/science/article/pii/0031920181900467); 
+ Moon: a radial Moon model from [Weber et al., 2011](https://science.sciencemag.org/content/331/6015/309) with its 3D crust; 
+ RT_MarsDWAK: a radial Mars from [Khan et al., 2016](https://www.sciencedirect.com/science/article/pii/S0031920116300875) model with its 3D crust as well as its centrifugal acceleration due to its rotation; 

The discontinuities (surfaces and fluid-solid boundaries) are pre-computed. The data are stored in ../discontinuities/. You can use the scripts under folders in ../discontinuities/ to create your own discontinuities. 

### Building similar models with various degrees of freedom 
Prior to building a planetary model, we utilize a radial model that is stored in ../../radialmodels/prem3L_noocean.mat. 
Here, we use CONST/CONST_mesh.m as an example. You will need to set the path to your mesh and model with its name 
~~~
fmesh  = './output/CONST3k/CONST_1L_3k'; % note that CONST_1L_3k is the model name. 
~~~
You can then set up the finite element polynomial order 
~~~
pOrder  = 1; %(choose 1 or 2)
~~~
You then set up the value to control the degrees of freedom 
~~~
a = 5e8; % you may play with this value to see the changes. 
~~~
Note that a here denotes the maximum volumn of any elements in km^3. The value _a_ controls the total number of the elements via applying a maximum tetrhedron volume constraint, please see the manual in TetGen. 
You can then load a unit surface mesh 
~~~
load ../../discontinuities/data/Sph392.mat; % you have other options, see ../../discontinuities/data/. 
~~~
i.e., you may use 
~~~
load ../../discontinuities/data/Sph3k6.mat; or Sph6k, Sph15k, etc.  
~~~
The surface mesh provides basic resolution of the surface. You can then simply run 
~~~
run CONST_mesh
~~~
to obtain both meshes and models. 

### Computing the reference gravity 
Once we obtain the mesh and model, we compute the reference gravity, **if needed**, using the fast multipole method. 
Similarly, you can set up the path to the model and the finite element polynomial order 
~~~
fmesh  = './output/CONST3k/CONST_1L_3k'; % note that CONST_1L_3k is the model name. 
pOrder  = 1; %(choose 1 or 2)
~~~
The scaling factor is used for visualization only 
~~~
scaling = 6.371*10^3; % here we use the radial of the Earth. 
~~~
It will not affect the computation. You can then run 
~~~
run Gravity
~~~
If you want to compute the reference gravity for a huge model, the memory consumption may be an issue. 
Please see GravityETC/gravity.m for a simple implementation to reduce the memory costs. 
We split many bodies into several small groups and compute the interactions among them. 

Under GravityETC/, the effects of centrifugal force due to rotation are included as well. Please edit the 
period (in hours) for the target planet. 
