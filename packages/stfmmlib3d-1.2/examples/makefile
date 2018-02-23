
stfmm3dpart:
	make -f stfmm3dpart.make -j4

stfmm3dtria:
	make -f stfmm3dtria.make -j4

lfmm3dpartquad:
	make -f lfmm3dpartquad.make -j4

st3green:
	make -f st3green.make -j4



stfmm3dpart-openmp:
	make -f stfmm3dpart.make -j4 HOST=linux-gfortran-openmp

stfmm3dtria-openmp:
	make -f stfmm3dtria.make -j4 HOST=linux-gfortran-openmp

lfmm3dpartquad-openmp:
	make -f lfmm3dpartquad.make -j4 HOST=linux-gfortran-openmp


clean: 
	rm -f *.o
 
distclean: 
	rm -f int2 fort.* *.o

