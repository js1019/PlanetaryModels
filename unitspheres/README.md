## Build discontinuities for planets

### 1. Perfectly spherically symmetric surfaces 
Please check data/unit_sphere.m as a demo to build a surface via applying a surface distant function. 
In this test case, we simply use a unit sphere. 

### 2. Ellipsoidal surfaces
Please check EarthEllp/ and MarsEllp/ as demos to build ellipsoidal surfaces via using a surface distant function. 

### 3. 3D crusts 
Please check mooncrust/ or Marscrust/ as demos to build a 3D crust, including topography and thickness. 
The lunar crust model is based on the results of two NASA missions: the Lunar Orbiter Laser Altimeter [(LOLA)](https://lola.gsfc.nasa.gov/) and Gravity Recovery and Interior Laboratory [(GRAIL)](https://www.nasa.gov/mission_pages/grail/main/index.html). 
The Mars crust model is based on the results of the NASA mission, the Mars Orbiter Laser Altimeter [(MOLA)](https://attic.gsfc.nasa.gov/mola/). 

### Reference 
~~~
@article{smith2010initial,
  title={{Initial observations from the lunar orbiter laser altimeter (LOLA)}},
  author={Smith, David E and Zuber, Maria T and Neumann, Gregory A and Lemoine, Frank G and Mazarico, Erwan and Torrence, Mark H and McGarry, Jan F and Rowlands, David D and Head, James W and Duxbury, Thomas H and others},
  journal={Geophysical Research Letters},
  volume={37},
  number={18},
  year={2010},
  publisher={Wiley Online Library}
}
~~~
for [the lunar topography](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010GL043751). 
~~~
@article{wieczorek2013crust,
  title={{The crust of the Moon as seen by GRAIL}},
  author={Wieczorek, Mark A and Neumann, Gregory A and Nimmo, Francis and Kiefer, Walter S and Taylor, G Jeffrey and Melosh, H Jay and Phillips, Roger J and Solomon, Sean C and Andrews-Hanna, Jeffrey C and Asmar, Sami W and others},
  journal={Science},
  volume={339},
  number={6120},
  pages={671--675},
  year={2013},
  publisher={American Association for the Advancement of Science}
}
~~~
for [the lunar crust](https://science.sciencemag.org/content/339/6120/671). 
~~~
@article{goossens2017evidence,
  title={{Evidence for a low bulk crustal density for Mars from gravity and topography}},
  author={Goossens, Sander and Sabaka, Terence J and Genova, Antonio and Mazarico, Erwan and Nicholas, Joseph B and Neumann, Gregory A},
  journal={Geophysical research letters},
  volume={44},
  number={15},
  pages={7686--7694},
  year={2017},
  publisher={Wiley Online Library}
}
~~~
for [the Mars crust](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL074172). 
