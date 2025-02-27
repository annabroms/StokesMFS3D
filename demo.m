clear;
close all

%%  Generate center coordinates for the particles
P = 10; %number of bodies
delta = 1; %smallest particle particle distance 
q = [0 0 0; 2+delta 0 0]; %Size P x 3

%random configurations
%L = 10; 
%q = set_position(P,L,delta); %Random in a qube or in a layer, with minimum
%distance
[q,B] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
  
fmm = 0; %only activate if many particles (say, more than 40)

%% Solve resistance problem first (given velocities)
disp('Start with resistance: ')
Uref = rand(6*P,1); 

%Note, for resistance, the number of GMRES iters will grow with P. For both
%resistance and mobility, GMRES iters will increase with decreasing delta.

[Fvec,it_res,lambda_norm_res,err_res] = solve_resistance(q,Uref,fmm); 

[U,it_mob,lambda_norm_mob,err_mob]  = solve_mobility(q,Fvec,fmm); 

%might want to change proxy radius a little to get fair 2-way error, when solving resistance followed by mobility. 
% Use Rg as extra argument to solve_mobility (see commented code below)

disp('Two way error')
norm(U-Uref,inf)/norm(Uref,inf) 
disp('Residual in mobility problem')
err_mob
disp('Residual in resistance problem')
err_res


% %% Solve mobility problem first 
% disp('Start with mobility: ')
% Fref = rand(6*P,1); 
% Rg = 0.63;
% U = solve_mobility(q,Fref,fmm,1.01*Rg);
% 
% [Fvec,it,lambda_norm2] = solve_resistance(q,U,fmm,Rg); 
% 
% disp('Two way error')
% norm(Fvec-Fref,inf)/norm(Fref,inf)



