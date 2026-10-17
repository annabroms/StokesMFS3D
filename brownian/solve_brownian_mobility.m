function [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,nout,wout,Y,UU,LL,Kin,Tblock,Fvec,opt,fluct_vel,R,E0)
%SOLVE_BROWNIAN_MOBILITY Solve a Stokes mobility problem with fluctuating velocity field for a configuration of ellipsoidal particles using MFS.
%This uses the non-standard DLP/SLP formulation (TG) for the solve.
%
%   [U, iters, lambda_norm] = SOLVE_BROWNIAN_MOBILITY(q, rvec_in, rvec_out, nout, wout, Fvec, opt, fluct_vel)
%
%   Given the forces and torques applied to each particle, the function computes the resulting translational and rotational
%   velocities.
%
%   INPUTS:
%       q         - P × 3 matrix of particle center positions.
%       rvec_in   - NP × 3 vector of collocation points on particle surfaces (stacked).
%       rvec_out  - MP × 3 vector of proxy source points (stacked).
%       nout      - MP × 3 DLP directions/normals on rvec_out.
%       wout      - MP × 1 quadrature weights on rvec_out.
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

if nargin < 13
    R = eye(3);
    E0 = [1 1 1];
    fluct_vel = [];
elseif nargin < 14
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 15
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 16
    E0 = [1 1 1];
end



%% Use solve_mobility_with_DLP with precomputed factors and optional noise
opt.precond.Y = Y;
opt.precond.UU = UU;
opt.precond.LL = LL;
opt.precond.Kin = Kin;
opt.precond.Tblock = Tblock;
if isfield(opt,'outer_force') && logical(opt.outer_force)
    if isfield(opt,'ellipsoid') && logical(opt.ellipsoid)
        if ~iscell(R)
            if isequal(size(R),[3 3])
                R = repmat({R},P,1);
            else
                error('solve_brownian_mobility:badRotations', ...
                    'R must be a cell array of 3-by-3 rotations for ellipsoids.');
            end
        end
        rin_pre = (R{1}' * (rvec_in(1:N,:) - q(1,:))')';
        rout_pre = (R{1}' * (rvec_out(1:M,:) - q(1,:))')';
        nout_pre = (R{1}' * nout(1:M,:)')';
        q_pre = [0 0 0];
    else
        rin_pre = rvec_in(1:N,:);
        rout_pre = rvec_out(1:M,:);
        nout_pre = nout(1:M,:);
        q_pre = q(1,:);
    end
    opt.precond.Fmap = getDLPForceMap( ...
        rin_pre,rout_pre,nout_pre,wout(1:M),q_pre,opt);
end
if isfield(opt,'Sblock')
    opt.precond.Sblock = opt.Sblock;
end
if isfield(opt,'TSblock')
    opt.precond.TSblock = opt.TSblock;
end
opt.compute_residual = false;
opt.transform_slip = false;

% Optional fluctuating boundary-velocity contribution, already expressed
% on the proxy grid and added to the mobility right-hand side.
if isempty(fluct_vel)
    u_slip = [];
else
    u_slip = sqrt(2/opt.dt)*fluct_vel;
end

[U, iters, lambda_norm] = solve_mobility_with_DLP( ...
    q,rvec_in,rvec_out,nout,wout,Fvec,u_slip,opt,R,E0);

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
opt.add_rank1 = false;
opt.outer_force = false;

N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P;
nout = rvec_out - kron(q,ones(M,1));
nout = nout ./ vecnorm(nout,2,2);
wout = (4*pi/M) * ones(P*M,1);
Tblock = stokes_DLP_mat(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:)-q(1,:));
opt_precond = opt;
opt_precond.inner_only = true;
[Y,UU,LL,Kin] = oneBodyPrecondMobDLP( ...
    rvec_in(1:N,:), rvec_out(1:M,:), opt_precond, q(1,:), wout(1:M), nout(1:M,:));

Fref = rand(6*P,1);
fluct_vel = [];

[U1,it1,ln1] = solve_brownian_mobility( ...
    q,rvec_in,rvec_out,nout,wout,Y,UU,LL,Kin,Tblock,Fref,opt,fluct_vel);
[U2,it2,ln2,err2] = solve_mobility_with_DLP( ...
    q,rvec_in,rvec_out,nout,wout,Fref,[],opt);

rel_err = norm(U1-U2,inf)/norm(U2,inf);

fprintf('\nSelf-test summary (solve_brownian_mobility)\n');
fprintf('  P = %d, delta = %.3g, N = %d, M = %d\n', P, delta, N, M);
fprintf('  gmres_tol = %.1e\n', opt.gmres_tol);
fprintf('  Brownian (no fluct_vel): iters = %d, lambda_norm = %.3e\n', it1, ln1);
fprintf('  Mobility DLP/SLP:     iters = %d, lambda_norm = %.3e, uerr = %.3e\n', it2, ln2, err2);
fprintf('  Relative rigid-body motion error = %.3e\n', rel_err);

end
