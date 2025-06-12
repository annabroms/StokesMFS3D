function [Fvec,iters,lambda_norm,err_res] = solve_resistance(q,U,fmm,Rp,N)
%SOLVE_RESISTANCE Solve a Stokes resistance problem for spherical particles.
%
%   [Fvec, iters, lambda_norm, err_res] = SOLVE_RESISTANCE(q, U, fmm, Rp,N)
%
%   Computes the hydrodynamic forces and torques on spherical particles
%   immersed in Stokes flow, given prescribed translational and rotational
%   velocities. The function is based on the Method of Fundamental Solutions (MFS)
%   with proxy source points placed on inner spherical surfaces and boundary
%   velocity conditions enforced at a slightly larger number of collocation
%   points at the true particle surfaces.
%
%   INPUTS:
%       q       - Px3 matrix of particle centers, with P the number of particles.
%       U       - 6P x 1 vector of prescribed translational and angular velocities.
%                 Format: [u1; omega1; u2; omega2; ...], each ui and omegai is 3x1.
%       fmm     - Logical flag: if true, use FMM3D acceleration.
%       Rp      - (Optional) Radius of the proxy surface around each particle.
%                 Default is 0.68.
%       N       - (Optional) Number of proxy source points per particle.
%                 Default is 700.
%
%   OUTPUTS:
%       Fvec        - 6P x 1 vector of computed forces and torques on each particle.
%                     Format matches input velocity vector.
%       iters       - Number of GMRES iterations used to reach tolerance.
%       lambda_norm - Infinity norm of the source density vector; a large value
%                     may indicate loss of accuracy or ill-conditioning.
%       err_res     - Maximum relative error in velocity on the particle surface.
%
%   METHOD OVERVIEW:
%       - Builds inner (proxy) and outer (collocation) grids for each particle.
%       - Constructs velocity boundary conditions from given velocities.
%       - Applies a one-body preconditioned GMRES solver to solve for source strengths.
%       - Computes forces and torques 
%       - Validates accuracy by evaluating flow at checkpoint nodes and
%       comparing to known Dirichlet boundary data.
%
%   NOTES:
%       - Assumes all particles are identical and spherical.
%       - The function can visualize the configuration if 'visualise' is enabled.
%       - Uses spherical design nodes
%       - Matrix-free GMRES 
%
%   DEPENDENCIES:
%       - init_MFS, getDesignGrid, getKmat, oneBodyPrecondRes,
%         helsing_gmres, getFlow
%
%   EXAMPLE USAGE:
%       q = [0 0 0; 3 0 0]; % Two particles separated by delta = 1.
%       U = [1; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0]; 
%       [F, iters, lambda_norm, err] = solve_resistance(q, U, false);
%
%   Anna Broms, June 12, 2025

%% Parameters 
% Handle optional inputs Rg (proxy radius) and N (number of proxy points).
% If not provided, use default values.
if nargin < 4
    Rp = 0.68; 
    N = 700;
elseif nargin < 5
    N = 700;
end

% Initialize MFS options 
opt = init_MFS(N);
opt.Rp = Rp;
opt.fmm = fmm;          % Enable or disable FMM acceleration
opt.maxit = 600;        % GMRES max iterations
gmres_tol = 1e-7;       % GMRES tolerance

P = size(q,1); %number of spheres


%% Particle geometry
% Build inner (proxy) and outer (collocation) surfaces for a unit sphere.
% These will later be translated to each particle center.
[rin, rout] = getDesignGrid(opt.Rp, opt);

% Get number of proxy and collocation points
N = size(rin,1);
M = size(rout,1);

rvec_in = [];   % proxy source points on everybody
rvec_out = [];  % collocation points on surfaces of everybody

for k = 1:P
    rvec_in = [rvec_in; rin + q(k,:)];     % translate inner points
    rvec_out = [rvec_out; rout + q(k,:)];  % translate outer points
end


%% Visualize geometry
% Optional block for displaying the particle configuration 
visualise = 1;
if visualise
    [XX,YY,ZZ] = sphere(12);
    r = 1; % assumed unit sphere
    figure()
    for k = 1:size(q,1)
        hSurface = surf(r*XX+q(k,1), r*YY+q(k,2), r*ZZ+q(k,3));
        set(hSurface, 'FaceColor', [1 0 0], ...
                      'FaceAlpha', 0.9, ...
                      'FaceLighting', 'gouraud', ...
                      'EdgeColor', 'none');
        hold on
        axis equal
    end
    camlight
end


%% Assign RHS in resistance problem
Kout = getKmat(rout,[0,0,0]);
%For each particle, get data at surface, given rigid body motion
for k = 1:P
    u_bndry((k-1)*3*M+1:3*k*M) = Kout*U((k-1)*6+1:k*6);
end

%% Compute preconditioning. It's enough to do this for a single particle 
%as it's assumed that everyone has the same shape.

[Yk,UUk] = oneBodyPrecondRes(rin,rout);

%The format is used to prepare for the case when different shapes are is
%use
Y{1} = Yk;
UU{1} = UUk; 
disp('Done preconditioning')

if opt.profile
    memorygraph('label','start matvec resistance');
end

%% Solve problem
[mu_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvecStokesMFS(x,rvec_in,rvec_out,q,UU,Y,opt,1,[]),u_bndry',3*size(rvec_out,1),opt.maxit,gmres_tol,0);

if opt.profile
    memorygraph('label','done matvec resistance, remap and determine force')
end

%% Determine source strengths on proxy sources from the solution at the boundary,
% lambda <- mu. Then, determine forces and torques on particles, given lambda

Fvec = zeros(6*size(q,1),1);
Kin = getKmat(rin,[0,0,0]);
for i = 1:P
    lambda_i = Y{1}*(UU{1}*(mu_gmres((i-1)*3*M+1:i*3*M)));
    lambda_gmres(3*(i-1)*N+1:i*3*N) = lambda_i;
    Fvec(6*(i-1)+1:6*i) = Kin'*lambda_i; 
end

% lambda_norm gives a sanity check on the source distribution. If large,
% the representation for the problem is not optimal.
lambda_norm = norm(lambda_gmres, inf);

%% Estimate surface residual
% Evaluate how well the flow generated by the sources matches the prescribed boundary velocity (from bc).
if opt.profile
    memorygraph('label','done solving resistance, compute velocities')
end

% Get new nodes for evaluating velocity residuals
opt.des_n = 2000; %fine grid
[~, rout_check] = getDesignGrid(opt.Rp, opt);
n_check = size(rout_check,1); 

% Assemble all check points
rcheck = [];
for k = 1:P
    x = q(k,:) + rout_check;
    rcheck = [rcheck; x];
end

% Assign expected velocity at check points using prescribed U
ucheck = zeros(n_check*3*P,1);
for k = 1:P
    Kcheck = getKmat(rcheck(n_check*(k-1)+1:k*n_check,:), q(k,:));
    ucheck((k-1)*3*n_check+1:3*k*n_check) = Kcheck * U((k-1)*6+1:k*6);
end

% Evaluate flow from solution to resistance problem
ubdry = getFlow(lambda_gmres, rvec_in, rcheck, opt);

% Compute relative residual
uerr_vec = vecnorm(reshape(ucheck - ubdry, 3, []), 2, 1) ./ ...
           max(vecnorm(reshape(ucheck, 3, []), 2, 1));
err_res = max(uerr_vec);



end