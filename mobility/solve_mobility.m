function [U, iters, lambda_norm, uerr] = solve_mobility(q,rvec_in,rvec_out,Fvec, opt,R,E0)
%SOLVE_MOBILITY Solve a Stokes mobility problem for a configuration of ellipsoidal particles using MFS.
%
%   [U, iters, lambda_norm, uerr] = SOLVE_MOBILITY(q, rvec_in, rvec_out, Fvec, opt, R, E0)
%
%   Given the forces and torques applied to each particle, the function computes the resulting translational and rotational
%   velocities.
%
%   INPUTS:
%       q         - P × 3 matrix of particle center positions.
%       rvec_in   - (N*P) × 3 array of proxy source points (stacked by particle).
%       rvec_out  - (M*P) × 3 array of collocation points on particle surfaces (stacked by particle).
%       Fvec      - 6P × 1 vector of applied forces and torques, format: [F1; T1; F2; T2; ...].
%       opt       - Struct containing solver options (e.g., gmres tolerance, fmm flag).
%       R         - P x 1 cell array of rotation matrices for the P
%                   particles 
%       E0        - 1 × 3 vector of semiaxes [a, b, c] of the ellipsoidal particles.
%
%   OUTPUTS:
%       U           - 6P × 1 vector of resulting rigid body velocities: [u1; omega1; u2; omega2; ...].
%       iters       - Number of GMRES iterations until convergence.
%       lambda_norm - Infinity norm of the final density vector (for diagnostic use).
%       uerr        - Maximum relative residual of the velocity field on the surface.
%
%   METHOD OVERVIEW:
%       - Builds MFS representation from proxy and collocation surfaces.
%       - Uses a "completion source" to represent force and torque.
%       - Computes rigid body motion from the computed source density.
%       - Determines the surface residuals in new points.
%
%   DEPENDENCIES:
%       init_MFS, getDesignGrid, getCompletionSource, matvecStokesMFS,
%       oneBodyPrecondMob, helsing_gmres, getKmat, getStokesletFlow,
%       rotate_vector, ellipsoid_param, setupsurfquad, getTractionFast
%
%   See also: SOLVE_RESISTANCE
%
%   NOTES:
%       - If opt.plot is true, a visualization of the geometry is produced.
%       - Surface velocity and traction visualization is computed only when opt.plot is true.
%
%   Anna Broms, June 13, 2025

if nargin==0, test_solve; 
    return; end


if nargin < 6
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 7
    E0 = [1 1 1];
end

P = size(q,1);

if ~isfield(opt,'gmres_verbose')
    opt.gmres_verbose = 1;
end


%% One-body preconditioning
N = size(rvec_in,1)/P; %number of sources per particle
M = size(rvec_out,1)/P; %numer of collocation points per particle

opt.N = N;
opt.M = M;

if opt.ellipsoid
    [rin_block, rout_block] = get_body_frame_block_data( ...
        rvec_in(1:N,:), rvec_out(1:M,:), q(1,:), R{1});
    opt.Sblock = stokes_SLP_mat(rin_block, rout_block);
else
    opt.Sblock = stokes_SLP_mat(rvec_in(1:N,:), rvec_out(1:M,:));
end

%Create pseudoinverse of self-interaction matrix,
if opt.ellipsoid
    [Y,UU,LL,Kin,~] = oneBodyPrecondMob((R{1}'*rvec_in(1:N,:)')',...
        (R{1}'*rvec_out(1:M,:)')',q(1,:));
else
    [Y,UU,LL,Kin,~] = oneBodyPrecondMob(rvec_in(1:N,:),...
        rvec_out(1:M,:),q(1,:));
end

%The format is used to prepare for the case when different shapes are is
%use
Yii{1} = Y;
UUii{1} = UU; 

%% Assemble completion source, given force and torque

lambda_vec = []; 

for k = 1:P

    %Create right hand side, given forces and torques on the particles
    F = Fvec(6*(k-1)+1:6*(k-1)+3);
    T = Fvec(6*(k-1)+4:6*k);
    
    if opt.ellipsoid
        Rk = R{k};
        lambda_k = getCompletionSource(Rk'*F,Rk'*T,Kin);
        lambda_vec = [lambda_vec; rotate_vector(lambda_k,Rk)];
    else
        lambda_k = getCompletionSource(F,T,Kin);
        lambda_vec = [lambda_vec; lambda_k];
    end

    
 
end

%% Get flow field due to completion source.
uvec = getStokesletFlow(lambda_vec,rvec_in,rvec_out,opt); 
uvec = -uvec;


%% Solve for source strengths
[x_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvecStokesMFS(x,rvec_in,rvec_out,q,UUii,Yii,opt,0,R,LL),uvec,3*size(rvec_out,1),opt.maxit,opt.gmres_tol,opt.gmres_verbose,0);

%check residual
%abs_res = norm(matvecStokesMFS(x_gmres,rvec_in,rvec_out,q,UUii,Yii,opt,0,R,LL)-uvec);


%% Map back to the sought density in source points, determine rigid body velocities 
U = zeros(6*P,1);
lambda_gmres = zeros(N*P*3,1);
for i = 1:P
    if opt.ellipsoid
        temp_i = Y*(UU*(rotate_vector(x_gmres((i-1)*3*M+1:i*3*M),R{i}')));
        lambda_i = rotate_vector(temp_i,R{i});
        Kin_i = getKmat(rvec_in(N*(i-1)+1:N*i,:),q(i,:));
        U(6*(i-1)+1:6*i) = -Kin_i'*lambda_i; 
    else
        lambda_i = Y*(UU*x_gmres((i-1)*3*M+1:i*3*M));
        U(6*(i-1)+1:6*i) = -Kin'*lambda_i; 
    end
    
    lambda_gmres((i-1)*3*N+1:i*3*N) = lambda_i;    
     
end
lambda_norm = norm(lambda_gmres,inf);


%% Check residual
% Get new nodes for evaluating velocity residuals
% Set up a baseline ellipsoid at the origin, axis-aligne
b = ellipsoid_param(E0(1),E0(2),E0(3));   % semi-axes a=1, b=1, c=1
% Discretize the ellipsoid surface with specified resolution
b = setupsurfquad(b, [46, 55]);  % [# nodes in t-direction, s-direction]


% Initialize array for all check points
rcheck = [];

% Assemble check points for all P particles
for k = 1:P
    if opt.ellipsoid
        x = q(k,:)+(R{k}*b.x)';
    else
        x = q(k,:)+b.x';
    end
    % Apply translation to the surface points of particle k
    
    rcheck = [rcheck; x];  % add transposed surface nodes
 
end

% Number of surface nodes on one ellipsoid
n_check = size(b.x, 2);


%Assign velocities at checkpoints
ucheck = zeros(n_check*3*P,1); 
for k = 1:P
    Kcheck = getKmat(rcheck(n_check*(k-1)+1:k*n_check,:),q(k,:));
    ucheck((k-1)*3*n_check+1:3*k*n_check) = Kcheck*U((k-1)*6+1:k*6);   
end

for i =1:P  
    if opt.ellipsoid
        densityK_particle = lambda_gmres(3*(i-1)*N+1:i*3*N)-...  %better
            rotate_vector(LL*rotate_vector(lambda_gmres(3*(i-1)*N+1:i*3*N),R{i}'),R{i})+lambda_vec(3*(i-1)*N+1:i*3*N);
    else
        densityK_particle = (eye(3*N)-LL)*lambda_gmres(3*(i-1)*N+1:i*3*N)+lambda_vec(3*(i-1)*N+1:i*3*N);
    end
    densityK(3*(i-1)*N+1:i*3*N) = densityK_particle;
end

%get flow and compare RHS and LHS of representation
ubdry_check = getStokesletFlow(densityK,rvec_in,rcheck,opt);
uerr_vec = vecnorm(reshape(ucheck-ubdry_check,3,[]),2,1)/max(vecnorm(reshape(ucheck,3,[]),2,1));
uerr = max(uerr_vec);

%% Optional visualization: color by surface velocity and traction magnitudes
if opt.plot
    ubdry_surf = getStokesletFlow(densityK, rvec_in, rvec_out, opt);
    umag = vecnorm(reshape(ubdry_surf,3,[]),2,1).';
    nvec = rvec_out - kron(q, ones(M,1));
    nvec = nvec ./ vecnorm(nvec,2,2);
    traction = getTractionFast(densityK, rvec_in, rvec_out, nvec, opt);
    tmag = vecnorm(reshape(traction,3,[]),2,1).';

    plot_surface_scalar(rvec_out, M, P, umag, ...
        'Mobility: surface velocity magnitude', 'parula');
    plot_surface_scalar(rvec_out, M, P, tmag, ...
        'Mobility: traction magnitude on surface', 'hot');
end



end

function test_solve
% Self-test: resistance -> mobility -> compare rigid-body motion.
rng(5); %reproducable
close all

P = 8; %number of bodies
delta = 1; %smallest particle particle distance
[q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta

fmm = 0; %only activate if many particles (say, more than 40)

Uref = rand(6*P,1);

[rvec_in,rvec_out,opt] = init_spheres(q);
opt.fmm = fmm;
opt.lr = 0;
opt.gmres_tol = 1e-10;
opt.plot = 1; 

[Fvec,it_res,lambda_norm_res,err_res] = solve_resistance(q,rvec_in,rvec_out,Uref, opt);
[U,it_mob,lambda_norm_mob,err_mob]  = solve_mobility(q,rvec_in,rvec_out,Fvec, opt);

rel_err = norm(U-Uref,inf)/norm(Uref,inf);

fprintf('\nSelf-test summary (solve_mobility)\n');
fprintf('  P = %d, delta = %.3g, N = %d, M = %d\n', P, delta, size(rvec_in,1)/P, size(rvec_out,1)/P);
fprintf('  gmres_tol = %.1e\n', opt.gmres_tol);
fprintf('  Resistance: iters = %d, lambda_norm = %.3e, rel surf residual = %.3e\n', it_res, lambda_norm_res, err_res);
fprintf('  Mobility:   iters = %d, lambda_norm = %.3e, rel surf residual = %.3e\n', it_mob, lambda_norm_mob, err_mob);
fprintf('  Relative rigid-body motion error = %.3e\n', rel_err);

alignfigs; 

end
