% Test of long-range preconditioning for spheres usng deflation
clear;
close all;

%% Generate geometry
P = 200; %number of bodies
delta = 1; %smallest particle particle distance 
%q = [0 0 0; 2+delta 0 0; 0 2+delta 0]+[1 1 2]; %center coordiante matrix for P particles, x,y,z: size P x 3
%q = [0 0 0];
%random configurations
%L = 10; %set size of domain
%q = set_position(P,L,delta); %Random in a qube or in a layer, with minimum
%distance
[q,B] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
  
fmm = 0; %only activate if many particles (say, more than 40)

%% Solve resistance problem first (given velocities)
rng(6);
disp('Start with resistance: ')
Uref = rand(6*P,1); 
Fref = rand(6*P,1);

%Note, for resistance, the number of GMRES iters will grow with P without long-range precond. For both
%resistance and mobility, GMRES iters will increase with decreasing delta.
Rp = 0.2;
N = 100;

[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N);
opt.fmm = fmm;
opt.lr = 0; % no long range precond
disp('...solve without deflation...')

[Fvec0,it_res0,lambda_norm_res0,err_res0] = solve_resistance(q,rvec_in,rvec_out,Uref, opt); 
% Now, with long range precond
opt.lr = 1;  
disp('...solve with deflation 1...')
[Fvec1,it_res1,lambda_norm_res1,err_res1] = solve_resistance(q,rvec_in,rvec_out,Uref, opt);
opt.lr = 2; 
disp('...solve with deflation 2...')
[Fvec,it_res2,lambda_norm_res2,err_res2] = solve_resistance(q,rvec_in,rvec_out,Uref, opt); 
Rp = 1-1.05*(1-opt.Rp);

%Compare to the number of iterations needed for mobility
disp('...solve mobility problem...')
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N); %use a very coarse resolution for now
[U,it_mob,lambda_norm_mob,err_mob]  = solve_mobility(q,rvec_in,rvec_out,Fvec, opt); 

%might want to change proxy radius a little to get fair 2-way error, when solving resistance followed by mobility. 
% Use Rp as extra argument to solve_mobility (see commented code below)

disp('Two way error')
norm(U-Uref,inf)/norm(Uref,inf) 
disp('Residual in mobility problem')
err_mob
disp('Residual in resistance problem')
err_res0
err_res1
err_res2