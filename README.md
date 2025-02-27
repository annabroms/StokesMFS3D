
# StokesMFS3D

Solves the Stokes resistance and mobility problems using the method of fundamental solutions. See demo.m. This is a minimal working example. I recommend putting a sufficiently large separation, delta, between the particles to guarantee accuracy. The use of images is not implemented for the mobility solve, and hence, only the basic algorithm with a proxy surface of source points is used both for resistance and mobility.

See publications:
* Accurate close interactions of Stokes spheres using lubrication-adapted image systems, A. Broms, A.H. Barnett, A.-K. Tornberg, JCP, (2024) https://doi.org/10.1016/j.jcp.2024.113636
* A Method of Fundamental Solutions for Large-Scale 3D Elastance and Mobility Problems, A. Broms, A.H. Barnett, A.-K. Tornberg, under review in ACOM, https://www.arxiv.org/abs/2409.04215

# Dependencies
* FMM3D from Flatiron (if solving problems with many particles). I think I normally would switch at, say, 40 particles.
* SE_unified from Joar Bagge. Needed for evaluating direct sums of the fundamental solutions. If no images are in use, which is currently the case, this would only be needed for Stokeslets. You will then basically only need to compile SE0P_Stokeslet_direct.c.
* Spherical design nodes. See README in the geometry folder.

Disclaimer! This is a private repo under construction. The code is not very organised.
