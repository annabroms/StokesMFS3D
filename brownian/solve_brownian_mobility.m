function [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec,opt,fluct_vel,R,E0)
%SOLVE_BROWNIAN_MOBILITY Solve a Stokes mobility problem with fluctuating velocity field for a configuration of ellipsoidal particles using MFS.
%This uses the non-standard DLP/SLP formulation (TG) for the solve.
%
%   [U, iters, lambda_norm] = SOLVE_BROWNIAN_MOBILITY(q, rvec_in, rvec_out, Fvec, opt, fluct_vel)
%
%   Given the forces and torques applied to each particle, the function computes the resulting translational and rotational
%   velocities.
%
%   INPUTS:
%       q         - P × 3 matrix of particle center positions.
%       rvec_in   - NP × 3 vector of collocation points on particle surfaces (stacked).
%       rvec_out  - MP × 3 vector of proxy source points (stacked).
%       Y         - Factor in one-body pseudoinverse
%       UU        - Factor in one-body pseudoinverse, such that Y*U' is
%                   TG_L^{(ii}).
%       LL        - Projection matrix for sources not to contribute to net
%                   force or torque.
%       Tblock    - DLP matrix for single body
%       Fvec      - 6P × 1 vector of applied forces and torques, format: [F1; T1; F2; T2; ...].
%       opt       - Struct containing solver options (e.g., gmres tolerance, fmm flag).
%       fluct_vel - Vector of size 3NP x 1 containing the fluctuating velocity
%                   field on the proxy surface (e.g. sqrt(C)*dW). If this
%                   is already scaled by sqrt(2/dt), set opt.fluct_vel_scaled = true.
%       R         - P x 1 cell array of rotation matrices for the P
%                   ellipsoids [NOT YET SUPPORTED]
%       E0        - 1 × 3 vector of semiaxes [a, b, c] of the ellipsoidal particles.
%
%   OUTPUTS:
%       U           - 6P × 1 vector of resulting rigid body velocities: [u1; omega1; u2; omega2; ...].
%       iters       - Number of GMRES iterations until convergence.
%       lambda_norm - Infinity norm of the final density vector (for diagnostic use).
%       uerr        - Maximum relative residual of the velocity field on the surface.
%
%   METHOD OVERVIEW:
%       - Builds MFS representation based on compositiong TG (stresslet
%       flow with Stokeslet flow data from proxy forces).
%       - Uses a "completion source" to represent net force and torque.
%       - Computes rigid body motion from the computed source density.
%       - Determines the surface residuals in new points.
%
%   DEPENDENCIES:
%       solve_mobility_with_DLP, oneBodyPrecondMobDLP, stokes_DLP_mat, getKmat
%
%   See also: LARGE_ELLIPSOID_EX, SOLVE_RESISTANCE
%
%   Anna Broms, Feb 7, 2026

if nargin==0, test_solve; 
    return; end

P = size(q,1);
N = size(rvec_in,1)/P; %number of sources per particle
M = size(rvec_out,1)/P; %numer of collocation points per particle

opt.M = M; 
opt.N = N;

if nargin < 11
    R = eye(3);
    E0 = [1 1 1];
    fluct_vel = [];
elseif nargin < 12
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 13
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 14
    E0 = [1 1 1];
end



%% Use solve_mobility_with_DLP with precomputed factors and optional noise
opt.precond.Y = Y;
opt.precond.UU = UU;
opt.precond.LL = LL;
opt.precond.Kin = Kin;
opt.precond.Tblock = Tblock;
opt.compute_residual = false;

% Optional fluctuating velocity contribution (added to TG*lambda_0)
if isempty(fluct_vel)
    opt.extra_uvec_T = [];
else
    opt.extra_uvec_T = -sqrt(2/opt.dt)*fluct_vel;
end

[U, iters, lambda_norm] = solve_mobility_with_DLP(q,rvec_in,rvec_out,Fvec,opt,R,E0);

end

function test_solve
% Self-test: compare solve_brownian_mobility (no fluct_vel) to solve_mobility_with_DLP.
rng(5); %reproducable

P = 2;
delta = 2;
q = [0 0 0; 2+delta 0 0];

Rp = 0.15;
N = 50;
a = 2;

[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a);
opt.fmm = 0;
opt.lr = 0;
opt.gmres_tol = 1e-10;

N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P;
Tblock = stokes_DLP_mat(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:)-q(1,:));
[Y,UU,LL,Kin] = oneBodyPrecondMobDLP(rvec_in(1:N,:), rvec_out(1:M,:), q(1,:));

Fref = rand(6*P,1);
fluct_vel = [];

[U1,it1,ln1] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fref,opt,fluct_vel);
[U2,it2,ln2,err2] = solve_mobility_with_DLP(q,rvec_in,rvec_out,Fref,opt);

rel_err = norm(U1-U2,inf)/norm(U2,inf);

fprintf('\nSelf-test summary (solve_brownian_mobility)\n');
fprintf('  P = %d, delta = %.3g, N = %d, M = %d\n', P, delta, N, M);
fprintf('  gmres_tol = %.1e\n', opt.gmres_tol);
fprintf('  Brownian (no fluct_vel): iters = %d, lambda_norm = %.3e\n', it1, ln1);
fprintf('  Mobility DLP/SLP:     iters = %d, lambda_norm = %.3e, uerr = %.3e\n', it2, ln2, err2);
fprintf('  Relative rigid-body motion error = %.3e\n', rel_err);

end
