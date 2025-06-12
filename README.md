
# StokesMFS3D

Solves the Stokes resistance and mobility problems using the method of fundamental solutions. See demo.m. This is a minimal working example. It's recommended to put a sufficiently large separation, delta, between the particles to guarantee accuracy. In this demo, image enhancement, as presented in paper 1 below, is not yet in use, and we demonstrated the basic algorithms with a proxy surface of source points only, for both resistance and mobility.

See publications:
1. Accurate close interactions of Stokes spheres using lubrication-adapted image systems, A. Broms, A.H. Barnett, A.-K. Tornberg, JCP, (2024) https://doi.org/10.1016/j.jcp.2024.113636
2. A Method of Fundamental Solutions for Large-Scale 3D Elastance and Mobility Problems, A. Broms, A.H. Barnett, A.-K. Tornberg, to appear in ACOM (2025), https://www.arxiv.org/abs/2409.04215

# Dependencies
* FMM3D from Flatiron (if solving problems with many particles). (Recommended switch at ~40 particles with the default settings).
* Stokes_Direct. Needed for evaluating direct sums of the fundamental solutions. If no images are in use, which is currently the case, this would only be needed for Stokeslets. Only SE0P_Stokeslet_direct.c has to be compiled. Precompiled version exists in Stokes_Direct.
* Spherical design nodes. See README in the geometry folder.

Image enhancement and code for ellipsoidal particles will be added.

![My Image](cluster.png)
