## Build discontinuities of planets
In this folder, we provide scripts to create triangle meshes for discontinuities of the target planet.


### Perfectly spherically symmetric surfaces 
Please check data/unit_sphere.m as a demo to build a surface via applying a surface distant function. 
In this test case, we simply use a unit sphere. 

### Ellipsoidal surfaces
Please check EarthEllp/ and MarsEllp/ as demos to build ellipsoidal exterior and interior surfaces via using surface distant functions. 

### Topography and 3D crusts 
Please check mooncrust/ or Marscrust/ as demos to build a 3D crust, including topography and thickness. 
The lunar crust model is based on the results of two NASA missions: the Lunar Orbiter Laser Altimeter ([LOLA](https://lola.gsfc.nasa.gov/)) and Gravity Recovery and Interior Laboratory ([GRAIL](https://www.nasa.gov/mission_pages/grail/main/index.html)). 
The Martian crust model is based on the results of the NASA mission, the Mars Orbiter Laser Altimeter ([MOLA](https://attic.gsfc.nasa.gov/mola/)). 

## Reference 
Here are some references for modeling the 3D crusts. 
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
@article{smith1999global,
  title={The global topography of Mars and implications for surface evolution},
  author={Smith, David E and Zuber, Maria T and Solomon, Sean C and Phillips, Roger J and Head, James W and Garvin, James B and Banerdt, W Bruce and Muhleman, Duane O and Pettengill, Gordon H and Neumann, Gregory A and others},
  journal={Science},
  volume={284},
  number={5419},
  pages={1495--1503},
  year={1999},
  publisher={American Association for the Advancement of Science}
}
~~~
for [the Martian topography](https://science.sciencemag.org/content/284/5419/1495). 
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
for [the Martian crust](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL074172). 
