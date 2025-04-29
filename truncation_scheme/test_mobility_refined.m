clear;
close all
rng(6);
%%  Generate center coordinates for the particles
P = 8; %number of bodies
delta = 1; %smallest particle particle distance 
q = [0 0 0; 2+delta 0 0]; %Size P x 3
Lcut = 2+2*delta; 

%random configurations
%L = 10; 
%q = set_position(P,L,delta); %Random in a qube or in a layer, with minimum
%distance
[q,B] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
%q = [0 0 0];
  
fmm = 0; %only activate if many particles (say, more than 40)
fprintf('Using fmm %u\n',fmm)

%% Solve resistance problem first (given velocities)
% P = 1; 
% q = [0 0 0];
P = 3; 
q = [0 0 0; 2+delta 0 0; 0 100 0];


disp('Start with resistance: ')
Fvec = rand(6*P,1); 

%Note, for resistance, the number of GMRES iters will grow with P. For both
%resistance and mobility, GMRES iters will increase with decreasing delta.

[U,iters,lambda_norm,uerr] = solve_mobility_refined_RPY(q,Fvec,fmm,Lcut);
%[U1,it_mob1,lambda_norm_mob1,err_mob1]  = solve_mobility_refined(q,Fvec,fmm,Lcut);  

[U2,it_mob2,lambda_norm_mob2,err_mob2]  = solve_mobility(q,Fvec,fmm); 

%might want to change proxy radius a little to get fair 2-way error, when solving resistance followed by mobility. 
% Use Rg as extra argument to solve_mobility (see commented code below)

norm(U-U2)/norm(U2)

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



