clear;
close all
%%%%%% Demo for Stokes MFS in 3D for a random configuration of P spherical
% particles:
%
% We solve 1. A resistance problem with no-slip bc set from random RBM
% (rigid body motion). Returns vector of net forces/torques on every
% particle. 2. A mobility problem with the computed net forces/torques from
% 1. as input. Returns computed RBM.
%
% The "2-way" error is determined by comparing the given RBM (the input to
% 1. to the computed RBM from 2.
%
% This demo code can be modified to generate 
%
% * Example 1 in paper [1] (Accurate close interactions of Stokes spheres
% using lubrication-adapted image systems, JCP 2024). To do so, generate
% configurations with q = set_position(P,L,delta); below and solve only the
% resistance problem with given randomly sampled translational and angular
% velocities.
%
% * Example 2 in paper [1]. To do so, generate
% configurations with [q,B] = grow_cluster(P,delta) below and set boundary
% conditions from a background flow in the resistance problem.
% 
% * Examples 1 and 2 in paper [2] (A Method of Fundamental Solutions for Large-Scale 3D Elastance and Mobility
% Problems, ACOM, 2025). To do so, generate configurations with [q,B] =
% grow_cluster(P,delta) below and loop over delta or P. Solve a resistance
% problem, followed by a mobility problem. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate center coordinates for the particles
P = 20; %number of bodies
delta = 1; %smallest particle particle distance 
%q = [0 0 0; 2+delta 0 0]; %center coordiante matrix for P particles, x,y,z: size P x 3

%random configurations
%L = 10; %set size of domain
%q = set_position(P,L,delta); %Random in a qube or in a layer, with minimum
%distance
[q,B] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
  
fmm = 1; %only activate if many particles (say, more than 40)

%% Solve resistance problem first (given velocities)
disp('Start with resistance: ')
Uref = rand(6*P,1); 
Fref = rand(6*P,1);

%Note, for resistance, the number of GMRES iters will grow with P. For both
%resistance and mobility, GMRES iters will increase with decreasing delta.
[rvec_in,rvec_out,opt] = init_spheres(q);
opt.fmm = fmm;
[Fvec,it_res,lambda_norm_res,err_res] = solve_resistance(q,rvec_in,rvec_out,Uref, opt); 
Rp = 1-1.05*(1-opt.Rp);
[rvec_in,rvec_out,opt] = init_spheres(q,Rp);
[U,it_mob,lambda_norm_mob,err_mob]  = solve_mobility(q,rvec_in,rvec_out,Fvec, opt); 

%might want to change proxy radius a little to get fair 2-way error, when solving resistance followed by mobility. 
% Use Rp as extra argument to solve_mobility (see commented code below)

disp('Two way error')
norm(U-Uref,inf)/norm(Uref,inf) 
disp('Residual in mobility problem')
err_mob
disp('Residual in resistance problem')
err_res





