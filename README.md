# StokesMFS3D
Solves the Stokes resistance and mobility problems using the method of fundamental solutions. See demo.m. This is a minimal working example. I recommend putting a sufficiently large separation, delta, between the particles to guarantee accuracy. The use of images is not implemented for the mobility solve, and hence, only the basic algorithm with a proxy surface of source points is used. 

# Dependencies

* FMM3D from flatiron (if solving problems with many particles). I think I normally would switch at, say, 40 particles. 
* SE_unified from Joar Bagge. Needed for evaluating direct sums of the fundamental solutions. If no images are in use, this would be only the Stokeslets. 
