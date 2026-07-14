function [U, iters, lambda_norm, uerr] = solve_mobility_with_DLP(q,rvec_in,rvec_out,nout,wout,Fvec,opt,R,E0)
%SOLVE_MOBILITY_WITH_DLP Solve a Stokes mobility problem for a
%configuration of ellipsoidal particles using "non-standard" MFS where the
%ansatz is a u = TWG \lambda, with G a Stokeslet (SLP) mapping from proxy
%surface to true surface, W is a diagonal matrix of quadrature weights 
%and T a Stresslet (DLP) mapping from true surface back to proxy surface. 
%Everyhting is analogous to the standard solve_mobility function. This is
%written in preparation for simulations of Brownian motion, where it is
%beneficial to work with a symmetric saddle point system. 
%
%   [U, iters, lambda_norm, uerr] = SOLVE_MOBILITY_WITH_DLP(q, rvec_in, rvec_out, nout, wout, Fvec, opt, R, E0)
%
%   Given the forces and torques applied to each particle, the function computes the resulting translational and rotational
%   velocities.
%
%   INPUTS:
%       q         - P × 3 matrix of particle center positions.
%       rvec_in   - (N*P) × 3 array of proxy source points (stacked by particle).
%       rvec_out  - (M*P) × 3 array of collocation points on particle surfaces (stacked by particle).
%       nout      - (M*P) × 3 DLP directions/normals on rvec_out.
%       wout      - (M*P) × 1 quadrature weights on rvec_out.
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
%       - Solves for a source density in the TWG representation.
%       - Applies block-diagonal right preconditioning. 
%       - Computes rigid body motion from the computed source density.
%       - Determines the surface residuals at check points.
%
%   DEPENDENCIES:
%       init_MFS, getDesignGrid, getCompletionSourceFromMap, matvecStokesMFS_DLP,
%       oneBodyPrecondMob_DLP, helsing_gmres, getKmat, getStokesletFlow,
%       getStressletFlow, stokes_DLP_mat, rotate_vector, ellipsoid_param,
%       setupsurfquad
%
%   See also: SOLVE_RESISTANCE_WITH_DLP, SOLVE_MOBILITY
%
%   Anna Broms, Nov 28, 2025

if nargin==0, test_solve; 
    return; end

if nargin < 7
    error('solve_mobility_with_DLP requires q, rvec_in, rvec_out, nout, wout, Fvec, and opt.');
end

if nargin < 8
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 9
    E0 = [1 1 1];
end

P = size(q,1);

if ~isfield(opt,'compute_residual')
    opt.compute_residual = true;
end
if ~isfield(opt,'gmres_verbose')
    opt.gmres_verbose = 1;
end
if ~isfield(opt,'ellipsoid')
    opt.ellipsoid = false;
end
if ~isfield(opt,'inner_only')
    opt.inner_only = false;
end
opt.inner_only = logical(opt.inner_only);
if ~isfield(opt,'symmetrize_weighted')
    opt.symmetrize_weighted = false;
end
opt.symmetrize_weighted = logical(opt.symmetrize_weighted);
if ~isfield(opt,'add_rank1')
    opt.add_rank1 = false;
end
opt.add_rank1 = logical(opt.add_rank1);
if ~isfield(opt,'rank1_scale')
    opt.rank1_scale = 1;
end
if ~isfield(opt,'outer_force')
    opt.outer_force = false;
end
opt.outer_force = logical(opt.outer_force);
if ~isfield(opt,'debug')
    opt.debug = false;
end
if ~isfield(opt,'profile')
    opt.profile = false;
end
if opt.ellipsoid && ~iscell(R)
    if isequal(size(R),[3 3])
        R = repmat({R},P,1);
    else
        error('solve_mobility_with_DLP:badRotations', ...
            'R must be a cell array of 3-by-3 rotations for ellipsoids.');
    end
end


%% One-body preconditioning
N = size(rvec_in,1)/P; %number of sources per particle
M = size(rvec_out,1)/P; %numer of collocation points per particle

if abs(N-round(N)) > eps || abs(M-round(M)) > eps
    error('rvec_in and rvec_out must contain the same number of nodes for each particle.');
end
N = round(N);
M = round(M);

if size(nout,2) ~= 3 || size(nout,1) ~= M*P
    error('solve_mobility_with_DLP:badNormals', ...
        'nout must be a (M*P)-by-3 array.');
end
wout = full(wout(:));
if numel(wout) ~= M*P
    error('solve_mobility_with_DLP:badWeights', ...
        'wout must contain one quadrature weight for each row of rvec_out.');
end

opt.N = N;
opt.M = M;

nin = [];
if opt.add_rank1
    nin = get_inner_normals(rvec_in,q,opt,R);
    opt.nin = nin;
end

if opt.ellipsoid
    rin_block = (R{1}' * (rvec_in(1:N,:) - q(1,:))')';
    rout_block = (R{1}' * (rvec_out(1:M,:) - q(1,:))')';
else
    rin_block = rvec_in(1:N,:);
    rout_block = rvec_out(1:M,:);
end

nin_body = [];
if opt.add_rank1
    nin_body = nin(1:N,:);
    if opt.ellipsoid
        nin_body = (R{1}' * nin_body')';
    end
end

% Prepare self-correction block needed in block-diagonal right precond
opt.Sblock = stokes_SLP_mat(rin_block, rout_block);
n_self = nout(1:M,:);
if opt.ellipsoid
    n_self = (R{1}' * n_self')';
end
Tself = stokes_DLP_mat(rout_block,rin_block,n_self);

Wself = diag(repelem(wout(1:M),3));
opt.TWSblock = Tself * Wself * opt.Sblock;
if opt.symmetrize_weighted
    opt.TWSblock = 0.5 * (opt.TWSblock + opt.TWSblock');
end
if opt.add_rank1
    opt.TWSblock = opt.TWSblock + normal_rank1_block(nin_body,opt.rank1_scale);
end

%Create pseudoinverse of self-interaction matrix (or use precomputed, which 
% saves computational time in any dynamic setting)
if isfield(opt,'precond') && ~isempty(opt.precond)
    Y = opt.precond.Y;
    UU = opt.precond.UU;
    LL = opt.precond.LL;
    Kin = opt.precond.Kin;
    if opt.outer_force
        if isfield(opt.precond,'Fmap')
            Fmap_body = opt.precond.Fmap;
        elseif isfield(opt.precond,'Kforce')
            Fmap_body = opt.precond.Kforce;
        else
            error('solve_mobility_with_DLP:missingPrecondForceMap', ...
                'opt.outer_force=true with opt.precond requires precond.Fmap.');
        end
    else
        Fmap_body = Kin;
    end
    if isfield(opt.precond,'Sblock')
        opt.Sblock = opt.precond.Sblock;
    end
    if isfield(opt.precond,'TSblock')
        opt.TWSblock = opt.precond.TSblock;
    end
    if isfield(opt.precond,'TWSblock')
        opt.TWSblock = opt.precond.TWSblock;
    end
else
    opt_precond = opt;
    if opt.ellipsoid
        n_precond = (R{1}' * nout(1:M,:)')';
        [Y,UU,LL,Kin,~,Fmap_body] = oneBodyPrecondMobDLP( ...
            rin_block,rout_block,opt_precond,[0 0 0],wout(1:M),n_precond,nin_body);
    else
        [Y,UU,LL,Kin,~,Fmap_body] = oneBodyPrecondMobDLP(rvec_in(1:N,:),...
            rvec_out(1:M,:),opt_precond,q(1,:),wout(1:M),nout(1:M,:),nin_body);
    end
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
        lambda_k = getCompletionSourceFromMap(Rk'*F,Rk'*T,Fmap_body);
        lambda_vec = [lambda_vec; rotate_vector(lambda_k,Rk)];
    else
        lambda_k = getCompletionSourceFromMap(F,T,Fmap_body);
        lambda_vec = [lambda_vec; lambda_k];
    end

    
end

%% Get flow field due to completion source: TWG*lambda_0
uvec_T = -apply_weighted_B( ...
    lambda_vec,rvec_in,rvec_out,nout,wout,opt,opt.symmetrize_weighted,nin);
if isfield(opt,'extra_uvec_T') && ~isempty(opt.extra_uvec_T)
    uvec_T = uvec_T + opt.extra_uvec_T;
    warning('check extra flow field')
end

% Debug: look at system matrix
if opt.debug
    s = length(rvec_in)*3;
    e = zeros(s,1);
    syst_mat = zeros(s);
    for i = 1:s
        i
        e(:) = 0;
        e(i) = 1;
        syst_mat(:,i) = matvecStokesMFS_DLP(e,rvec_in,rvec_out,q,wout,nout,UUii,Yii,opt.TWSblock,opt,0,R,LL);
    end
    [V,D] = eig(syst_mat);
    d = diag(D);
    figure()
    plot(real(d),imag(d),'+');
end


%% Solve for source strengths
[x_gmres,iters,~,~] = helsing_gmres(@(x) matvecStokesMFS_DLP(x,rvec_in,rvec_out,q,wout,nout,UUii,Yii,opt.TWSblock,opt,0,R,LL),uvec_T,3*size(rvec_in,1),opt.maxit,opt.gmres_tol,opt.gmres_verbose,0);


%% Map back to the sought density in source points, determine rigid body velocities 
U = zeros(6*P,1);
lambda_gmres = zeros(N*P*3,1);
for i = 1:P
    if opt.ellipsoid
        temp_i = Y*(UU*(rotate_vector(x_gmres((i-1)*3*N+1:i*3*N),R{i}')));
        lambda_i = rotate_vector(temp_i,R{i});
        if opt.outer_force
            lambda_body_i = rotate_vector(lambda_i,R{i}');
            U_body = -Fmap_body' * lambda_body_i;
            U(6*(i-1)+1:6*i) = [R{i} * U_body(1:3); R{i} * U_body(4:6)];
        else
            Kin_i = getKmat(rvec_in(N*(i-1)+1:N*i,:),q(i,:));
            U(6*(i-1)+1:6*i) = -Kin_i' * lambda_i;
        end
    else
        lambda_i = Y*(UU*x_gmres((i-1)*3*N+1:i*3*N));
        if opt.outer_force
            inds_n = N*(i-1)+1:N*i;
            inds_m = M*(i-1)+1:M*i;
            Fmap_i = getDLPForceMap( ...
                rvec_in(inds_n,:),rvec_out(inds_m,:),nout(inds_m,:),wout(inds_m),q(i,:),opt);
            U(6*(i-1)+1:6*i) = -Fmap_i' * lambda_i;
        else
            U(6*(i-1)+1:6*i) = -Kin' * lambda_i;
        end
    end
    
    lambda_gmres((i-1)*3*N+1:i*3*N) = lambda_i;    
     
end
lambda_norm = norm(lambda_gmres,inf);

if opt.compute_residual
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
ubdry = getStokesletFlow(densityK,rvec_in,rcheck,opt);
uerr_vec = vecnorm(reshape(ucheck-ubdry,3,[]),2,1)/max(vecnorm(reshape(ucheck,3,[]),2,1));
uerr = max(uerr_vec);
else
    uerr = NaN;
end

end



function y = apply_node_weights(x,w)
y = reshape(x,3,[]);
y = y .* reshape(w,1,[]);
y = y(:);
end

function y = apply_weighted_B(x,rvec_in,rvec_out,nout,wout,opt,symmetrize_weighted,nin)
y = apply_weighted_B_core(x,rvec_in,rvec_out,nout,wout,opt);
if symmetrize_weighted
    yt = apply_weighted_BT_core(x,rvec_in,rvec_out,nout,wout,opt);
    y = 0.5 * (y + yt);
end
if isfield(opt,'add_rank1') && logical(opt.add_rank1)
    if nargin < 8 || isempty(nin)
        error('solve_mobility_with_DLP:missingInnerNormals', ...
            'opt.add_rank1=true requires inner normals for the weighted operator.');
    end
    y = y + apply_normal_rank1(x,nin,opt.rank1_scale);
end
end

function y = apply_weighted_B_core(x,rvec_in,rvec_out,nout,wout,opt)
u = getStokesletFlow(x,rvec_in,rvec_out,opt);
u = apply_node_weights(u,wout);
y = getStressletFlow(rvec_out,rvec_in,nout,u,numel(wout),opt);
end

function y = apply_weighted_BT_core(x,rvec_in,rvec_out,nout,wout,opt)
traction_opt = opt;
traction_opt.fmm = 0;
u = getTractionFast(x,rvec_in,rvec_out,nout,traction_opt);
u = apply_node_weights(u,wout);
y = getStokesletFlow(u,rvec_out,rvec_in,opt);
end

function test_solve
rng(5); %reproducable

P = 10; %number of bodies
delta = 1; %smallest particle particle distance 
%q = [0 0 0; 2+delta 0 0]; %center coordiante matrix for P particles, x,y,z: size P x 3
%q = [0 0 0]; 

%random configurations
[q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
  
fmm = 0; %only activate if many particles (say, more than 40)

%% Solve mobility problem (given forces/torques)
Fref = rand(6*P,1); 

%Test first with very low resolution
Rp = 0.30;
N = 100; 
a = 1.2; 
%a = 2; 


% Rp = 0.15; %Proxy sphere radius -- very coarse resolution
% N = 50;  % Number of proxy sources
% a = 2; 
%a = 2; %or play with SVD truncation level

%Rp = 0.68; %proxy radius
%N = 700; % approximate number of proxy sources on every particle


[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a);
%[rvec_in,rvec_out,opt] = init_spheres(q);
M = size(rvec_out,1)/P;
nout = rvec_out - kron(q,ones(M,1));
nout = nout ./ vecnorm(nout,2,2);
wout = (4*pi/M) * ones(P*M,1);
opt.fmm = fmm;
opt.lr = 0; % no long range precond for the standard method, to make it more comparable to the DLP version.
opt.gmres_tol = 1e-12; %does this matter for the accuracy? 
opt.inner_only = 0; % use the weighted inner-grid one-body preconditioner
opt.debug = 1; 
opt.add_rank1 = 0; 
opt.outer_force = false;
[Uvec,it_mob,lambda_norm_mob,err_mob] = solve_mobility_with_DLP( ...
    q,rvec_in,rvec_out,nout,wout,Fref,opt);
opt.gmres_tol = 1e-10; 
[Uvec2,it_mob2,lambda_norm_mob2,err_mob2] = solve_mobility(q,rvec_in,rvec_out,Fref, opt);
rel_err = norm(Uvec-Uvec2,inf)/norm(Uvec2,inf);

fprintf('\nTest summary (solve_mobility_with_DLP)\n');
fprintf('  P = %d, delta = %.3g, N = %d, M = %d\n', P, delta, size(rvec_in,1)/P, size(rvec_out,1)/P);
fprintf('  Rp = %.3g, a = %.3g, gmres_tol = %.1e\n', Rp, a, opt.gmres_tol);
fprintf('  DLP/SLP combo: iters = %d, lambda_norm = %.3e, uerr = %.3e\n', it_mob, lambda_norm_mob, err_mob);
fprintf('  Standard solve: iters = %d, lambda_norm = %.3e, uerr = %.3e\n', it_mob2, lambda_norm_mob2, err_mob2);
fprintf('  Relative rigid-body motion error vs standard = %.3e\n', rel_err);

end
