function [U, iters, lambda_norm, uerr] = solve_mobility_sym_left(q,rvec_in,rvec_out,nout,wout,Fvec,opt,R,E0)
%SOLVE_MOBILITY_SYM_LEFT Solve a 3D Stokes mobility problem written with the symmetrized 
% formulation using left preconditioning.
%
%   [U,iters,lambda_norm,uerr] = SOLVE_MOBILITY_SYM_LEFT( ...
%       q,rvec_in,rvec_out,nout,wout,Fvec,opt,R,E0)
%
%   This mobility solver uses the weighted DLP/SLP operator
%   B = T*W*S, with the physical source density lambda as the GMRES
%   unknown. The one-body block preconditioner is applied from the left.
%   The recompleted operator is B*(I-L) plus a rigid-body completion.
%
%   Set opt.inner_only = true to use the inner-grid completion Kin*Kin'.
%   Set opt.inner_only = false to use the quadrature approximation 
%   T*W*Kout*Kin' for the completion. The default is true, as it results in higher accuracy.
%   Set opt.symmetrize_weighted = true to replace B by 0.5*(B+B'), using
%   the matrix-free adjoint action from the Brownian TWS code.

if nargin == 0
    test_solve;
    return;
end

if nargin < 7
    error('solve_mobility_sym_left requires q, rvec_in, rvec_out, nout, wout, Fvec, and opt.');
end

P = size(q,1);
N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P;

if abs(N-round(N)) > eps || abs(M-round(M)) > eps
    error('rvec_in and rvec_out must contain the same number of nodes for each particle.');
end
N = round(N);
M = round(M);

if nargin < 8 || isempty(R)
    R = repmat({eye(3)},P,1);
elseif ~iscell(R)
    if isequal(size(R),[3 3]) && P == 1
        R = {R};
    else
        R = repmat({eye(3)},P,1);
    end
end

if nargin < 9 || isempty(E0)
    E0 = [1 1 1];
end

if isempty(opt)
    opt = struct();
end
opt = set_default(opt,'gmres_verbose',1);
opt = set_default(opt,'compute_residual',true);
opt = set_default(opt,'profile',false);
opt = set_default(opt,'ellipsoid',false);
opt = set_default(opt,'inner_only',true);
opt = set_default(opt,'symmetrize_weighted',false);
opt = set_default(opt,'add_rank1',false);
opt = set_default(opt,'rank1_scale',1);
opt = set_default(opt,'outer_force',false);
opt = set_default(opt,'debug',false);
opt.inner_only = logical(opt.inner_only);
opt.symmetrize_weighted = logical(opt.symmetrize_weighted);
opt.add_rank1 = logical(opt.add_rank1);
opt.outer_force = logical(opt.outer_force);
opt.N = N;
opt.M = M;
opt.P = P;

nout = validate_normals(nout,M,P);
wout = validate_weights(wout,M,P);
nin = [];
if opt.add_rank1
    nin = get_inner_normals(rvec_in,q,opt,R);
    opt.nin = nin;
end

symmetrize_weighted = logical(opt.symmetrize_weighted);

%% One-body weighted left preconditioner
[rin_body,rout_body,nout_body,nin_body] = get_first_body_frame_data( ...
    q,rvec_in,rvec_out,nout,nin,N,M,opt,R);


[Y,UU,LL,Kin,~,Fmap] = oneBodyPrecondMobDLP( ...
    rin_body,rout_body,opt,[0 0 0],wout(1:M),nout_body,nin_body);
ImLL = eye(3*N) - LL;

%% Assemble completion source for the prescribed force/torque
lambda_vec = zeros(3*N*P,1);
for k = 1:P
    F = Fvec(6*(k-1)+1:6*(k-1)+3);
    T = Fvec(6*(k-1)+4:6*k);

    if opt.ellipsoid
        lambda_k = getCompletionSourceFromMap(R{k}'*F,R{k}'*T,Fmap);
        lambda_k = rotate_vector(lambda_k,R{k});
    else
        lambda_k = getCompletionSourceFromMap(F,T,Fmap);
    end

    lambda_vec(3*N*(k-1)+1:3*N*k) = lambda_k;
end

rhs = -apply_weighted_B(lambda_vec,rvec_in,rvec_out,nout,wout,opt,symmetrize_weighted,nin);
rhs = apply_left_precond(rhs,Y,UU,N,P,opt,R);
%%
matvec = @(lambda) apply_left_precond( ...
    apply_weighted_B(project_mobility_source(lambda,ImLL,N,P,opt,R), ...
        rvec_in,rvec_out,nout,wout,opt,symmetrize_weighted,nin) + ...
    build_rigid_completion(lambda,Kin,Fmap,q,rvec_in,rvec_out,nout,wout,N,M,P,opt,R), ...
    Y,UU,N,P,opt,R);

% Debug: look at system matrix
if opt.debug
    s = length(rvec_in)*3;
    e = zeros(s,1);
    syst_mat = zeros(s);
    for i = 1:s
        i;
        e(:) = 0;
        e(i) = 1;
        syst_mat(:,i) = matvec(e);
    end
end


%% Solve the left-preconditioned mobility system.

[lambda_gmres,iters,~,~] = helsing_gmres( ...
    matvec,rhs,3*size(rvec_in,1),opt.maxit,opt.gmres_tol,opt.gmres_verbose,0);

%% Recover rigid-body velocities.
U = zeros(6*P,1);
for k = 1:P
    rows_n = 3*N*(k-1)+1:3*N*k;
    if opt.outer_force
        if opt.ellipsoid
            lambda_body = rotate_vector(lambda_gmres(rows_n),R{k}');
            U_body = -Fmap' * lambda_body;
            U(6*(k-1)+1:6*k) = [R{k} * U_body(1:3); R{k} * U_body(4:6)];
        else
            U(6*(k-1)+1:6*k) = -Fmap' * lambda_gmres(rows_n);
        end
    else
        inds_n = N*(k-1)+1:N*k;
        Kin_i = getKmat(rvec_in(inds_n,:),q(k,:));
        U(6*(k-1)+1:6*k) = -Kin_i' * lambda_gmres(rows_n);
    end
end

lambda_norm = norm(lambda_gmres,inf);

if opt.compute_residual
    uerr = compute_surface_residual(q,rvec_in,lambda_gmres,lambda_vec, ...
        LL,N,P,opt,R,E0,U);
else
    uerr = NaN;
end

end

function nout = validate_normals(nout,M,P)
nout = full(nout);
if size(nout,2) ~= 3 || size(nout,1) ~= M*P
    error('nout must be a (M*P)-by-3 array.');
end
end

function wout = validate_weights(wout,M,P)
wout = full(wout(:));
if numel(wout) ~= M*P
    error('wout must contain one quadrature weight for each row of rvec_out.');
end
end

function [rin_body,rout_body,nout_body,nin_body] = get_first_body_frame_data(q,rvec_in,rvec_out,nout,nin,N,M,opt,R)
rin_block = rvec_in(1:N,:);
rout_block = rvec_out(1:M,:);
nout_block = nout(1:M,:);
nin_body = [];
if ~isempty(nin)
    nin_body = nin(1:N,:);
end

if opt.ellipsoid
    rin_body = (R{1}' * (rin_block - q(1,:))')';
    rout_body = (R{1}' * (rout_block - q(1,:))')';
    nout_body = (R{1}' * nout_block')';
    if ~isempty(nin_body)
        nin_body = (R{1}' * nin_body')';
    end
else
    rin_body = rin_block - q(1,:);
    rout_body = rout_block - q(1,:);
    nout_body = nout_block;
end
end

function y = apply_weighted_B(x,rvec_in,rvec_out,nout,wout,opt,symmetrize_weighted,nin)
if symmetrize_weighted
    y = apply_sym_Amat(x,rvec_in,rvec_out,nout,wout,opt);
else
    y = apply_Amat(x,rvec_in,rvec_out,nout,wout,opt);
end
if isfield(opt,'add_rank1') && logical(opt.add_rank1)
    if nargin < 8 || isempty(nin)
        error('solve_mobility_sym_left:missingInnerNormals', ...
            'opt.add_rank1=true requires inner normals for the weighted operator.');
    end
    y = y + apply_normal_rank1(x,nin,opt.rank1_scale);
end
end

function y = apply_left_precond(x,Y,UU,N,P,opt,R)
y = zeros(size(x));
if opt.ellipsoid
    for k = 1:P
        rows_n = 3*N*(k-1)+1:3*N*k;
        x_body = rotate_vector(x(rows_n),R{k}');
        y(rows_n) = rotate_vector(Y * (UU * x_body),R{k});
    end
else
    xmat = reshape(x,3*N,P);
    ymat = Y * (UU * xmat);
    y = ymat(:);
end
end

function source = project_mobility_source(lambda,ImLL,N,P,opt,R)
source = zeros(size(lambda));
if opt.ellipsoid
    for k = 1:P
        rows_n = 3*N*(k-1)+1:3*N*k;
        lambda_body = rotate_vector(lambda(rows_n),R{k}');
        source_body = ImLL * lambda_body;
        source(rows_n) = rotate_vector(source_body,R{k});
    end
else
    lambda_mat = reshape(lambda,3*N,P);
    source_mat = ImLL * lambda_mat;
    source = source_mat(:);
end
end

function source = build_rigid_completion(lambda,Kin,Fmap,q,rvec_in,rvec_out,nout,wout,N,M,P,opt,R)
if ~opt.inner_only
    source = build_outer_rigid_completion(lambda,Fmap,q,rvec_in,rvec_out,nout,wout,N,M,P,opt,R);
    return;
end

source = zeros(size(lambda));
if opt.ellipsoid
    for k = 1:P
        rows_n = 3*N*(k-1)+1:3*N*k;
        lambda_body = rotate_vector(lambda(rows_n),R{k}');
        source_body = Kin * (Fmap' * lambda_body);
        source(rows_n) = rotate_vector(source_body,R{k});
    end
else
    lambda_mat = reshape(lambda,3*N,P);
    source_mat = Kin * (Fmap' * lambda_mat);
    source = source_mat(:);
end
end

function source = build_outer_rigid_completion(lambda,Fmap,q,rvec_in,rvec_out,nout,wout,N,M,P,opt,R)
uout = zeros(3*M*P,1);

for k = 1:P
    rows_n = 3*N*(k-1)+1:3*N*k;
    rows_m = M*(k-1)+1:M*k;

    if opt.ellipsoid
        lambda_body = rotate_vector(lambda(rows_n),R{k}');
        U_body = Fmap' * lambda_body;
        rout_body = (R{k}' * (rvec_out(rows_m,:) - q(k,:))')';
        Kout_body = getKmat(rout_body,[0 0 0]);
        uout(3*M*(k-1)+1:3*M*k) = rotate_vector(Kout_body * U_body,R{k});
    else
        Kout_i = getKmat(rvec_out(rows_m,:),q(k,:));
        uout(3*M*(k-1)+1:3*M*k) = Kout_i * (Fmap' * lambda(rows_n));
    end
end

uout = apply_quad_weights(uout,wout);
source = getStressletFlow(rvec_out,rvec_in,nout,uout,numel(wout),opt);
end

function uerr = compute_surface_residual(q,rvec_in,lambda_gmres,lambda_vec,LL,N,P,opt,R,E0,U)
b = ellipsoid_param(E0(1),E0(2),E0(3));
b = setupsurfquad(b,[46,55]);

rcheck = [];
for k = 1:P
    if opt.ellipsoid
        x = q(k,:) + (R{k} * b.x)';
    else
        x = q(k,:) + b.x';
    end
    rcheck = [rcheck; x]; %#ok<AGROW>
end

n_check = size(b.x,2);
ucheck = zeros(n_check*3*P,1);
for k = 1:P
    rows_check = n_check*(k-1)+1:k*n_check;
    Kcheck = getKmat(rcheck(rows_check,:),q(k,:));
    ucheck(3*n_check*(k-1)+1:3*n_check*k) = Kcheck * U(6*(k-1)+1:6*k);
end

densityK = zeros(size(lambda_gmres));
for k = 1:P
    rows_n = 3*N*(k-1)+1:3*N*k;
    lambda_i = lambda_gmres(rows_n);
    if opt.ellipsoid
        lambda_body = rotate_vector(lambda_i,R{k}');
        density_i = lambda_i - rotate_vector(LL * lambda_body,R{k}) + lambda_vec(rows_n);
    else
        density_i = (eye(3*N) - LL) * lambda_i + lambda_vec(rows_n);
    end
    densityK(rows_n) = density_i;
end

ubdry = getStokesletFlow(densityK,rvec_in,rcheck,opt);
uerr_vec = vecnorm(reshape(ucheck-ubdry,3,[]),2,1) ./ ...
    max(vecnorm(reshape(ucheck,3,[]),2,1));
uerr = max(uerr_vec);
end

function opt = set_default(opt,name,value)
if ~isfield(opt,name)
    opt.(name) = value;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_solve
rng(5);

P = 1;
delta = 2;
if P > 1
    [q,~] = grow_cluster(P,delta);
else
    q = [0 0 0];
end
Fref = rand(6*P,1);
Rp = 0.30;
Ndes = 100;
a = 1.2; 

[rvec_in,rvec_out,opt] = init_spheres(q,Rp,Ndes,a);
M = size(rvec_out,1)/P;
nout = rvec_out - kron(q,ones(M,1));
nout = nout ./ vecnorm(nout,2,2);
wout = (4*pi/M) * ones(P*M,1);

opt.fmm = 0;
opt.lr = 0;
opt.gmres_tol = 1e-10;
opt.gmres_verbose = 0;
opt.maxit = 400;
opt.compute_residual = true;
opt.symmetrize_weighted = true;
opt.debug = 1; 
opt.inner_only = false;
opt.add_rank1 = false;
opt.outer_force = false;

[U_left,it_left,lambda_norm_left,err_left] = ...
    solve_mobility_sym_left(q,rvec_in,rvec_out,nout,wout,Fref,opt);

opt_dlp = opt;
opt_dlp.inner_only = true;
opt_dlp.symmetrize_weighted = false;
opt_dlp.debug = 0; 
[U_dlp,it_dlp,lambda_norm_dlp,err_dlp] = ...
    solve_mobility_with_DLP(q,rvec_in,rvec_out,nout,wout,Fref,opt_dlp);
[U_std,it_std,lambda_norm_std,err_std] = ...
    solve_mobility(q,rvec_in,rvec_out,Fref,opt);

rel_left_dlp = norm(U_left-U_dlp,inf) / max(norm(U_dlp,inf),eps);
rel_left_std = norm(U_left-U_std,inf) / max(norm(U_std,inf),eps);

fprintf('\nSelf-test summary (solve_mobility_sym_left)\n');
fprintf('  P = %d, delta = %.3g, N = %d, M = %d\n', ...
    P,delta,size(rvec_in,1)/P,size(rvec_out,1)/P);
fprintf('  Rp = %.3g, a = %.3g, gmres_tol = %.1e\n',Rp,a,opt.gmres_tol);
fprintf('  Left weighted: iters = %d, lambda_norm = %.3e, uerr = %.3e\n', ...
    it_left,lambda_norm_left,err_left);
fprintf('  Right DLP:     iters = %d, lambda_norm = %.3e, uerr = %.3e\n', ...
    it_dlp,lambda_norm_dlp,err_dlp);
fprintf('  Standard:      iters = %d, lambda_norm = %.3e, uerr = %.3e\n', ...
    it_std,lambda_norm_std,err_std);
fprintf('  rel U error vs right DLP = %.3e\n',rel_left_dlp);
fprintf('  rel U error vs standard  = %.3e\n',rel_left_std);

assert(isfinite(err_left), ...
    'Left-preconditioned mobility solve produced a non-finite surface residual.');
assert(isfinite(rel_left_dlp), ...
    'Left-preconditioned and right-preconditioned DLP mobility comparison is non-finite.');

if err_left >= 1e-2
    warning('solve_mobility_sym_left:selfTestSurfaceResidual', ...
        'Left-preconditioned surface residual is %.3e for this test case.',err_left);
end
if rel_left_dlp >= 5e-2
    warning('solve_mobility_sym_left:selfTestDlpMismatch', ...
        'Left and right DLP mobility outputs differ by %.3e for this test case.',rel_left_dlp);
end
end
