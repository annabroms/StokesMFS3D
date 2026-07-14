function [Y, UU, LL, Kin, Kout, Fmap, Kouter] = oneBodyPrecondMobDLP(rin, rout, opt, q, wout, nout, nin)
%oneBodyPrecondMobDLP Determines pseudoinverse factors for the single body MFS
% Stokes mobility system matrix constructed as T*S, with T a
% double layer (stresslet) matrix and S a target-from-source matrix of Stokeslets
%
%   [Y, UU, LL, Kin] = oneBodyPrecondMobDLP(rin, rout, opt, q, wout, nout, nin)
%  
%   The function builds and factorizes a one-body matrix that does not
%   contribute to net force and torque and that utilizes an unused subspace (the range of Kin) to
%   express rigid body motions in terms of the unknown source vector,
%
%   INPUTS:
%       rin  - N x 3 matrix of proxy source points.
%       rout - M x 3 matrix of collocation points.
%       opt  - (Optional) struct. If opt.inner_only is true, use the
%              weighted inner-grid completion block
%              B*(I-LL) + Kin*Kin', where B = T*W*S.
%       q    - (Optional) 1 x 3 vector for the particle center. Default: [0 0 0].
%       wout - M x 1 quadrature weights for rout.
%       nout - M x 3 DLP directions/normals.
%       nin  - (Optional) N x 3 inner proxy normals. Required for
%              ellipsoids when opt.add_rank1 is true.
%
%   OUTPUTS:
%       U - Matrix of left singular vectors of the 1-body system matrix 
%           corresponding to retained singular values
%       Y - Product VS⁺, where:
%           - S⁺ is a diagonal matrix with entries 1/σ for retained
%           singular values of the 1-body system matrix
%           - V contains the corresponding right singular vectors
%       LL     - Force/torque projection matrix built from the selected
%                force map Fmap.
%       Kin    - Inner source force/torque map.
%       Kout   - Rigid body velocity map on rout.
%       Fmap   - Selected source force/torque map; equals Kin unless
%                opt.outer_force is true, then equals Kouter.
%       Kouter - Outer DLP force/torque map T*W*Kout.
%
%   DEPENDENCIES:
%       - stokes_SLP_mat, getKmat,
%       getPseudoFactors,oneBodyPrecondResDLP
%
% Anna Broms, Nov 28, 2025

if nargin<1
    self_test();
    return;
end
if nargin < 3
    opt = [];
end
if nargin < 4
    q = [];
end
if nargin < 5
    wout = [];
end
if nargin < 6
    nout = [];
end
if nargin < 7
    nin = [];
end

if nargin < 3 || isempty(opt)
    opt = struct();
    q = [0 0 0];
elseif isnumeric(opt)
    q_old = opt;
    if nargin >= 4 && isstruct(q)
        opt_old = q;
    else
        opt_old = struct();
    end
    opt = opt_old;
    q = q_old;
elseif isstruct(opt)
    if nargin < 4 || isempty(q)
        q = [0 0 0];
    end
else
    error('oneBodyPrecondMobDLP:badOpt', ...
        'Third argument must be an options struct.');
end

if ~isfield(opt,'inner_only')
    opt.inner_only = false;
end
if ~isfield(opt,'add_rank1')
    opt.add_rank1 = false;
end
if ~isfield(opt,'rank1_scale')
    opt.rank1_scale = 1;
end
if ~isfield(opt,'outer_force')
    opt.outer_force = false;
end
inner_only = logical(opt.inner_only);
symmetrize_weighted = isfield(opt,'symmetrize_weighted') && logical(opt.symmetrize_weighted);
add_rank1 = logical(opt.add_rank1);

% Build MFS and force/torque maps.
S    = stokes_SLP_mat(rin, rout);
[Fmap, Kin, Kout, Kouter, T, W] = getDLPForceMap(rin,rout,nout,wout,q,opt);

% Construct the force/torque projection matrix.
LL = Fmap * ((Fmap' * Fmap) \ Fmap');

% Build one-body operator that does not contribute to force/torque. The
% last term takes care of velocity constraint.
B = T * W * S;
if symmetrize_weighted
    B = 0.5 * (B + B');
end
if add_rank1
    if isempty(nin)
        nin = get_inner_normals(rin,q,opt,[]);
    end
    B = B + normal_rank1_block(nin,opt.rank1_scale);
end

if inner_only %RBM on inner surfaces
    Cmap = Kin;
else
    Cmap = Kouter;
end
A = B * (eye(size(LL)) - LL) + Cmap * Fmap';

% Compute pseudoinverse factors 
if isfield(opt,'tol')
    tol = opt.tol;
else
    tol = 1e-14;
end
visualise = 0; 
[Y, UU] = getPseudoFactors(A,tol,visualise);
UU = UU';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function self_test()
% Self-test for oneBodyPrecondMob and getCompletionSource.
% Verifies projection behavior.

% Setup proxy and boundary points, normals and quadrature weights
Rp = 0.7;
opt.des_n = 100;
opt.a_glob = 1.2;
opt.inner_only = 0; 

[rbase_in, rbase_out] = getDesignGrid(Rp, opt);  % Unit test geometry
M = size(rbase_out,1);
nout = rbase_out;
nout = nout ./ vecnorm(nout,2,2);
wout = (4*pi/M) * ones(M,1);

% -Compute pseudoinverse and projection matrix 
[~, ~, LL, Kin, ~, Fmap] = oneBodyPrecondMobDLP(rbase_in, rbase_out,opt,[0 0 0],wout,nout);

% Define prescribed net force and torque 
F = [1; 0; 0];   % Unit force in x
T = [0; 0; 1];   % Unit torque about z

% Construct completion source 
lambda0 = getCompletionSourceFromMap(F, T, Fmap);

% Check that (I - L) * lambda0 = 0 
I = eye(size(LL));
residual = norm((I - LL) * lambda0);
fprintf('Relative residual ||(I - L) * lambda0|| = %.2e\n', residual);

if residual > 1e-12
    warning('(I - L) * lambda0 is not small — test failed!');
end

%Verify that Fmap' * lambda0 = [F; T] 
checkFT = Fmap' * lambda0;
fprintf('Error in F: %g, T: %g\n', norm(checkFT(1:3)-F), norm(checkFT(4:6)-T));

%Check LL is a projection (i.e., L^2 = L) 
Lerr = norm(LL*LL - LL);
fprintf('Projection error ||L^2 - L|| = %.2e\n', Lerr);

M = size(rbase_out,1);
opt_inner = opt;
opt_inner.inner_only = true;
wout = (4*pi/M) * ones(M,1);
nout = rbase_out;
nout = nout ./ vecnorm(nout,2,2);
[Yin, UUin, LLin, Kinin] = oneBodyPrecondMobDLP( ...
    rbase_in, rbase_out, opt_inner, [0 0 0], wout, nout);
fprintf('Inner-only weighted preconditioner sizes: Y %dx%d, UU %dx%d\n', ...
    size(Yin,1),size(Yin,2),size(UUin,1),size(UUin,2));
assert(isequal(size(LLin),size(LL)), ...
    'Inner-only LL size changed unexpectedly.');
assert(isequal(size(Kinin),size(Kin)), ...
    'Inner-only Kin size changed unexpectedly.');

opt_outer = opt_inner;
opt_outer.outer_force = true;
opt_outer.add_rank1 = false;
[~, ~, LLout, ~, ~, Fmap_out, Kouter] = oneBodyPrecondMobDLP( ...
    rbase_in, rbase_out, opt_outer, [0 0 0], wout, nout);
lambda0_out = getCompletionSourceFromMap(F, T, Fmap_out);
lambda_test = randn(size(Fmap_out,1),1);
[~, ~, Kout_test, Kouter_test, T_test, W_test] = getDLPForceMap( ...
    rbase_in, rbase_out, nout, wout, [0 0 0], opt_outer);
proj_force_err = norm(Fmap_out' * (eye(size(LLout)) - LLout),'fro');
completion_err = norm(Fmap_out' * lambda0_out - [F; T]);
dense_force_err = norm(Kouter_test' * lambda_test - Kout_test' * W_test * T_test' * lambda_test);
selection_err = norm(Fmap_out - Kouter,'fro');
fprintf('Outer-force projector error = %.2e\n', proj_force_err);
fprintf('Outer-force completion error = %.2e\n', completion_err);
fprintf('Outer-force map selection error = %.2e\n', dense_force_err);
fprintf('Outer-force selected-map error = %.2e\n', selection_err);
if proj_force_err > 1e-11 || completion_err > 1e-11 || ...
        dense_force_err > 1e-12 || selection_err > 1e-12
    error('oneBodyPrecondMobDLP:outer_force_self_test', ...
        'outer_force projector or completion test failed.');
end

opt_rank = opt_inner;
opt_rank.add_rank1 = true;
nin = rbase_in ./ vecnorm(rbase_in,2,2);
[Yrank, UUrank] = oneBodyPrecondMobDLP( ...
    rbase_in, rbase_out, opt_rank, [0 0 0], wout, nout, nin);
fprintf('Rank-one weighted preconditioner sizes: Y %dx%d, UU %dx%d\n', ...
    size(Yrank,1),size(Yrank,2),size(UUrank,1),size(UUrank,2));

end
