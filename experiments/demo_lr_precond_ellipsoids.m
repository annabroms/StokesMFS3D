% Long-range preconditioning for ellipsoids
clear;
close all;
 
P = 40; %number of particles
Nv = 20; %sets resolution. This will give pour accuracy, but much faster run with lower resolution
%Nv = 40; %resolution used in paper [2]
delta = 0.5;
%delta = 1; % minimum particle separation, configuration will be generated at random
visualise = 1;
solve_res = 1; %solve both resistance and mobility
read_name = [];
%% Solve without deflation (long-range preconditioning) 
save_name = "test_without_deflation";
lr = 0; %first solve without long range preconditioning
rng(6);
ellipsoid_mobility_run(P,delta,Nv,visualise,solve_res,read_name,save_name,lr)
%% Solve with deflation 
read_name = save_name; %use same configuration and 
save_name = 'test_with_deflation';
lr = 2;
%rng(6); %want to create the same geometry, and solve with long range precond
ellipsoid_mobility_run(P,delta,Nv,visualise,solve_res,read_name,save_name,lr)

%Compare iter_r from these two runs