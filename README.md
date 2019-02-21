Planetary Model Builder 
================================================================
This repository provides codes to build a planet model on a tetrahedal mesh
as well as its reference gravity. It supports [the planetary normal mode computation](https://github.com/js1019/NormalModes).  

<img src="figs/PREM_vp.gif" width="425"/> <img src="figs/PREM_vs.gif" width="425"/>

In the above figures, examples of compressional and shear wave speed models of our Earth are illustrated. 

How to make it work for you? 
----------------------------------------------------------------
This repository uses several packages, including 
+ [TetGen](http://www.tetgen.org) by Hang Si;   
+ [distmesh](http://persson.berkeley.edu/distmesh/) by Per-Olof Persson;  
+ [several MATLAB scripts](https://github.com/tcew/nodal-dg) from Hesthaven & Warburton's book 'Nodal Discontinuous Galerkin Methods'; 
+ [fmmlib3d](https://github.com/zgimbutas/fmmlib3d) from Leslie Greengard and Zydrunas Gimbutas;  
+ [a MATLAB script](https://www.mathworks.com/matlabcentral/fileexchange/58002-write-binary-vtk-file-for-tetrahedral-grid-with-scalar-and-vector-data?s_tid=prof_contriblnk) for vtk users by Shawn Walker. 

We've put most of the original materials under different folders 
with some minor changes for this application.


You may follow the README.md under packages to compile what you need. 
If you only need to build a tetrahedal mesh, please compile tetgen; if you are a linux user, just go to packages/tetgen1.5.0 and type make. If you need to compute the reference gravity, please go to packages/fmmlib3d-1.2, and compile it.  


Build your planetary models
-----------------------------------------------------------------
Please check the scripts under demos/CONST for a constant ball model or demos/PREM 
for a reference earth model. Under demos/CONST, you may run 
~~~ 
run CONST_mesh; run Gravity;
~~~
to obtain the model and its reference gravity. Under demos/PREM, you may also run 
~~~
run PREM_mesh; run Gravity;
~~~
to get the Preliminary Reference Earth Model (PREM) and its reference gravity.

### A few remarks
+ _If the reference gravity is not needed, 
you do not have to run Gravity._
+ Here are two animations for [compressional wave speed of PREM](https://www.youtube.com/watch?v=4AeXhXGClcY) and [shear wave speed of PREM](https://www.youtube.com/watch?v=22yVo2G2e0k). 
+ We use [Paraview](https://www.paraview.org/) to visualize the results. 
+ To design your own models, you may change the settings, including the different discontinuities and model profiles. 


Furthermore
------------------------------------------------------------------
You can build more realistic models using similiar ideas! 

<img src="figs/CMI_94k-eps.png" width="425"/> <img src="figs/MIT4M_vp-eps.png" width="425"/>

The top right figure illustrates an Earth compressional wave speed model based on [MIT tomographic results](https://pubs.geoscienceworld.org/ssa/srl/article/79/3/384/367688/upper-mantle-heterogeneity-beneath-north-america) (Burdick et al. 2017) 
and [crust 1.0](https://igppweb.ucsd.edu/~gabi/crust1.html) (Laske et al. 2013). The top left one shows the topography of the Moho discontinuity under Tibetan Plateau.
If you would like to reproduce these models, you may need to download the model data and make minor changes in the scripts to obtain and visualize your results. 



Reference
-------------------------------------------------------------------
The repository provides scripts to generate planetary models for [our SuperComputing (SC'18) paper](https://dl.acm.org/citation.cfm?id=3291751), see below for details. 

~~~
@inproceedings{shi2018computing,
  title={Computing planetary interior normal modes with a highly parallel polynomial filtering eigensolver},
  author={Shi, Jia and Li, Ruipeng and Xi, Yuanzhe and Saad, Yousef and de Hoop, Maarten V},
  booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage, and Analysis, {SC}'18, Dallas, TX, November 11-16, 2018},
  pages={71:1--71:13},
  year={2018},
  organization={ACM/IEEE}
}
~~~


Contact 
-----------------------------------------------------------------
Please report issues under this repository. 
