Planetary Model Builder 
================================================================
<img src="https://img.shields.io/github/languages/code-size/js1019/PlanetaryModels.svg?style=social"><img src="https://img.shields.io/github/issues/js1019/PlanetaryModels.svg?style=social"> <img src="https://img.shields.io/github/forks/js1019/PlanetaryModels.svg?style=social"> <img src="https://img.shields.io/github/stars/js1019/PlanetaryModels.svg?style=social"> <img src="https://img.shields.io/github/license/js1019/PlanetaryModels.svg?style=social">

This repository provides scripts to build a planetary model on a deformable tetrahedal mesh
as well as its reference gravity. It supports and provides input files for the repository ["normal mode computation at planetary scales"](https://github.com/js1019/NormalModes).  

<p align="center">
   <img src="figs/PREM_vp.gif" width="375"/> <img src="figs/PREM_vs.gif" width="375"/>
</p>

The above are two animations for [compressional wave speed of PREM](https://www.youtube.com/watch?v=4AeXhXGClcY) and [shear wave speed of PREM](https://www.youtube.com/watch?v=22yVo2G2e0k). 

Contents
----------------------------------------------------------------
+ **Demos**: scripts to build 3D planetary models; 
+ **bin**: a script to run the Matlab scripts in the Demos/; 
+ **figs**: figures and animation of different planetary models; 
+ **modelbuilder**: scripts to build finite-element matrices and other information;
+ **packages**: external packages to build the mesh, compute the reference gravity and visualize the modeling results; 
+ **radialmodels**: scripts and data that contain information about radial planetary models from the existing literature; 
+ **discontinuities**: scripts to build discontinuities using unstructured triangles as well as some precomputed meshes; 
+ **visual**: scripts to visualize the normal modes computed from the repository ["normal mode computation at planetary scales"](https://github.com/js1019/NormalModes).

How to make it work for you? 
----------------------------------------------------------------
This repository uses several packages, including 
+ [TetGen](http://www.tetgen.org) by Hang Si;   
+ [distmesh](http://persson.berkeley.edu/distmesh/) by Per-Olof Persson;  
+ [several Matlab scripts](https://github.com/tcew/nodal-dg) from Hesthaven & Warburton's book 'Nodal Discontinuous Galerkin Methods'; 
+ [fmmlib3d](https://github.com/zgimbutas/fmmlib3d) from Leslie Greengard and Zydrunas Gimbutas;  
+ [a Matlab script](https://www.mathworks.com/matlabcentral/fileexchange/58002-write-binary-vtk-file-for-tetrahedral-grid-with-scalar-and-vector-data?s_tid=prof_contriblnk) for vtk users by Shawn Walker; 

We've put most of the original codes under different folders 
with some minor changes for this application.


You may follow the README.md under packages to compile what you need. 
If you only need to build a tetrahedral mesh, please compile TetGen; if you are a Linux user, just go to packages/tetgen1.5.0 and type make. If you need to compute the reference gravity, please check the readme under packages/ and compile packages/fmmlib3d-1.2.  


Build your planetary models
-----------------------------------------------------------------
Please check the scripts under Demos/CONST for a constant ball model, Demos/PREM 
for a standard Earth model, Demos/Moon for a Moon model and Demos/RT_MarsDWAK for a Mars model. 
Under these folders, you may run 
~~~ 
run {CONST/PREM/M6Ltopo/MarsDWAK}_mesh; run Gravity;
~~~
to obtain the model and its reference gravity. 
Note that the 3D crusts used in this work are based on the results of several NASA missions, including 
the Lunar Orbiter Laser Altimeter ([LOLA](https://lola.gsfc.nasa.gov/)), Gravity Recovery and Interior Laboratory ([GRAIL](https://www.nasa.gov/mission_pages/grail/main/index.html)) and 
the Mars Orbiter Laser Altimeter ([MOLA](https://attic.gsfc.nasa.gov/mola/)). 
Please see discontinuities/ for more details. 

### Illustration of other terrestrial planets
<p align="center">
   <img src="figs/M6Ltopo1M_visual.png" width="250"/> <img src="figs/Mars_topo94k_axis.png" width="240"/><img src="figs/MDWAK200k_vs.png" width="240"/>
</p>   
Here, we show a Moon shear wave speed model in the top right figure, the Mars topography in the top middle figure and a Mars shear wave speed model in the top right figure. These models can be reproduced by the scripts provided in the Demos. 


The corresponding Moon and Mars models can be found via [this IEEE data port](https://ieee-dataport.org/documents/mars-and-moon-models-used-reproducibility-challenge-student-cluster-competition-sc19). 

### Include more heterogeneity
You can build more realistic models using similiar ideas and following the above tips! You may need to understand what distmesh can do. It is quite simple, please see [its Demos](http://persson.berkeley.edu/distmesh/).  
<p align="center">
   <img src="figs/CMI_94k-eps.png" width="350"/> <img src="figs/MIT4M_vp-eps.png" width="350"/>
</p>

The top right figure illustrates an Earth compressional wave speed model based on [MIT tomographic results](https://pubs.geoscienceworld.org/ssa/srl/article/79/3/384/367688/upper-mantle-heterogeneity-beneath-north-america) (Burdick et al. 2017) 
and [crust 1.0](https://igppweb.ucsd.edu/~gabi/crust1.html) (Laske et al. 2013). The top left one shows the topography of the Moho discontinuity under Tibetan Plateau.
If you would like to reproduce these models, you may need to download the model data and make some minor changes in the scripts to obtain and visualize your results. 


### A few remarks 
+ We use [ParaView](https://www.paraview.org/) to visualize the results (**vtk or vtu files**). 
+ To design your own models, you may change the settings, including **different discontinuities and model profiles**. 
+ To insert discontinuities, such as topography, interior boundaries, please check the folder discontinuities/ and utilize **surface distance functions** to build your own meshes.  

Postprocess for visualizing the modes
------------------------------------------------------------------
You can use scripts in **visual/** to visualize your computed normal modes from the repository ["normal mode computation at planetary scales"](https://github.com/js1019/NormalModes). 
Please check visual/README.md for more details. 



Reference
-------------------------------------------------------------------
+ _**Theory, discretization and validation:**_ Shi, Jia, Ruipeng Li, Yuanzhe Xi, Yousef Saad, and Maarten V. de Hoop. "A non-perturbative approach to computing seismic normal modes in rotating planets." arXiv preprint arXiv:1906.11082 (2021), [the preprint](https://arxiv.org/abs/1906.11082).
+ _**Reproducibility, summary and discussion of the student cluster competition results:**_ Shi, Jia, Ruipeng Li, Yuanzhe Xi, Yousef Saad, and Maarten V. de Hoop. "Planetary normal mode computation: Parallel algorithms, performance, and reproducibility." IEEE Transactions on Parallel and Distributed Systems, 32, no. 11 (2021): 2609-2622, [the IEEE TPDS paper link](https://ieeexplore.ieee.org/abstract/document/9319555). 
+ _**Parallel performace and algorithms:**_ Shi, Jia, Ruipeng Li, Yuanzhe Xi, Yousef Saad, and Maarten V. de Hoop. "**Computing planetary interior normal modes with a highly parallel polynomial filtering eigensolver.**" In SC18: International Conference for High Performance Computing, Networking, Storage and Analysis, pp. 894-906. IEEE, 2018, [the SC18 paper link](https://dl.acm.org/citation.cfm?id=3291751).


Other related references can be found in modelbuilder/, packages/, radialmodels/ and  discontinuities/. 


Report
-----------------------------------------------------------------
Please let us know any issues of this repository. Contributions are welcome. 
