clear;
close all
%% Demo for Stokes MFS in 3D for a random configuration of P spherical particles: 
% We solve 
% 1. A resistance problem with no-slip bc set from random RBM (rigid body
% motion). Returns vector of net forces/torques on every particle. 
% 2. A mobility problem with the computed net forces/torques from 1. as
% input. Returns computed RBM.
%
% The "2-way" error is determined by comparing the given RBM (the input to
% 1. to the computed RBM from 2. 
%
% This demo code can be used to generate Examples 1 and 2 in paper [2] 
% (A Method of Fundamental Solutions for Large-Scale 3D Elastance and
% Mobility Problems). To do so, generate configurations with [q,B] = grow_cluster(P,delta);
% below and loop over delta or P. 
%%  Generate center coordinates for the particles
P = 40; %number of bodies
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

%Note, for resistance, the number of GMRES iters will grow with P. For both
%resistance and mobility, GMRES iters will increase with decreasing delta.

[Fvec,it_res,lambda_norm_res,err_res] = solve_resistance(q,Uref,fmm); 

[U,it_mob,lambda_norm_mob,err_mob]  = solve_mobility(q,Fvec,fmm); 

%might want to change proxy radius a little to get fair 2-way error, when solving resistance followed by mobility. 
% Use Rp as extra argument to solve_mobility (see commented code below)

disp('Two way error')
norm(U-Uref,inf)/norm(Uref,inf) 
disp('Residual in mobility problem')
err_mob
disp('Residual in resistance problem')
err_res


% %% Solve mobility problem first 
% disp('Start with mobility: ')
% Fref = rand(6*P,1); 
% Rp = 0.63;
% U = solve_mobility(q,Fref,fmm,1.01*Rp);
% 
% [Fvec,it,lambda_norm2] = solve_resistance(q,U,fmm,Rp); 
% 
% disp('Two way error')
% norm(Fvec-Fref,inf)/norm(Fref,inf)



