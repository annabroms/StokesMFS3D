function [U,iters,lambda_norm,uerr,lambda,details] = ...
    solve_mobility_wall_images(q,rin,rout,Fvec,opt)
%SOLVE_MOBILITY_WALL_IMAGES Mobility with wall ST_tD and close-pair SDR.
%   [U,iters,lambda_norm,uerr,lambda,details] = ...
%       solve_mobility_wall_images(q,rin,rout,Fvec,opt)
%   solve_mobility_wall_images() runs its own plotted particle self-test.
%
% q is P-by-3 physical centres; rin/rout are ordinary proxy/surface grids
% stacked by particle. N counts initial proxy points; M counts initial
% targets. Image lines and contact caps are appended internally; final M(k)
% and N_total(k) can differ by particle. N_image is points per image line.
% Fvec is 6P-by-1 [F;T] per particle ON THE FLUID; hydrodynamic particle
% drag is -Fvec. U returns [translation;angular velocity] per particle.
%
% lambda contains the total PHYSICAL coefficients, including completion.
% Default blocks are
% [proxy S;wall S;wall T_t;wall D;pair 1 S;pair 1 D;pair 1 R;...] per body.
% Each block interleaves Cartesian components by point. S has force units,
% T_t/R force*length and D force*length^2. details.lambda is dimensionless,
% with physical conversion details.sourceUnits. details.rin/rout are final
% physical grids and details.problem retains normalized geometry/factors.
% lambda_norm is norm(details.lambda,inf). uerr is the independent maximum
% surface velocity discrepancy/max rigid speed, absolute for zero motion;
% NaN when compute_residual=false. Exceeding 1e-6 is reported, never asserted.
%
% Kin' maps all coefficients, including pair rotlets, to six force/torque
% components. completion=Kin/(Kin'*Kin) maps six moments to coefficients.
% lambda_vec is the source carrying the prescribed moments. lambda_gmres
% is the auxiliary source reconstructed from collocation GMRES data. Its
% moment-carrying part (the main solver's LL action) is removed before
% adding lambda_vec; rigid velocity is -Kin'*lambda_i BEFORE projection.
%
% The rectangular self matrix is
%   A = Sblock+(Kout-Sblock*completion)*Kin'.
% Sblock maps local mixed sources to local complete-wall velocity; Kout maps
% rigid motion to surface velocity. Y{k}*(UU{k}*x) applies A's truncated
% pseudoinverse. UU holds retained LEFT singular vectors in rows; Y holds
% RIGHT singular vectors divided by retained singular values. U denotes
% physical rigid velocities throughout, never singular vectors.
%
% Options retain solve_mobility naming: tol, maxit, gmres_tol, gmres_verbose,
% fmm, fmm_tol, plot, compute_residual. See prepare_wall_images_problem for
% cap and unit controls.
% opt.viscosity=1 by default (mu is accepted as a legacy alias).
% opt.SVDsolve=false uses one-body preconditioning and helsing_gmres; true
% factors the global rectangular recompleted matrix using direct kernels,
% with no local SVDs, identity replacement or FMM dependency. In this mode
% iters=0, real_res/converged are NaN, and linearSystemResidual reports the
% global rectangular fitting residual rather than a GMRES residual.
% Select opt.wallImageFamilies="ST_tD" and
% opt.pairImageFamilies="SDR" independently (defaults). Accepted tokens:
% S,R,D,T_t,T_n,T_n2; T includes all three stresslet orientations. t is
% along the relevant wall/pair line; n/n2 are transverse. "none" disables
% one line type. Actual coefficient blocks are recorded in body{k}.basis.
% fmm=true is hybrid: ONE batched complete-wall S
% FMM plus blocked DIRECT T_t/D/R. Optional families retain quadratic body
% scaling. StokesHalfspaceFMM must be on the path; FMM3D is required for S
% acceleration. No slip/background flow, unequal spheres or overlaps.
% Use init_spheres for initial grids; explicitly select opt.fmm=1 if its
% third output is used (the legacy initializer defaults to direct).

if nargin == 0
    [U,iters,lambda_norm,uerr,lambda,details] = self_test();
    return
end
if nargin < 5, opt = struct(); end
validateattributes(Fvec,{'double'},{'column','real','finite','numel',6*size(q,1)});
problem = prepare_wall_images_problem(q,rin,rout,opt,'mobility');
opt = problem.options;
forceUnits = repmat([opt.mu*opt.radius*ones(3,1); ...
    opt.mu*opt.radius^2*ones(3,1)],problem.P,1);
Fvec = Fvec./forceUnits;

%% Completion source: carry all prescribed force and torque exactly.
lambda_vec = zeros(problem.sourceOffsets(end),1);
for k = 1:problem.P
    cols = problem.sourceOffsets(k)+(1:3*problem.N_total(k));
    lambda_vec(cols) = problem.body{k}.completion*Fvec(6*(k-1)+(1:6));
end
uvec = -apply_wall_images_velocity(lambda_vec,problem);

%% Solve for auxiliary sources by global SVD or one-body Helsing GMRES.
timer = tic;
if opt.SVDsolve
    [lambda_gmres,info] = solve_wall_images_svd(problem,uvec);
    iters = 0; resvec = []; real_res = NaN;
    problem.ranks = info.rank; problem.singularValues = {info.singularValues};
    problem.timing.globalAssemblySeconds = info.assemblySeconds;
    problem.timing.globalFactorizationSeconds = info.factorizationSeconds;
elseif ~any(uvec)
    x_gmres = zeros(size(uvec)); iters = 0; resvec = 0; real_res = 0;
else
    % The callback projects force/torque before the physical flow evaluation.
    matvec = @(x) wall_images_matvec(x,problem);
    [x_gmres,iters,resvec,real_res] = helsing_gmres(matvec,uvec, ...
        numel(uvec),min(opt.maxit,numel(uvec)),opt.gmres_tol, ...
        0,opt.gmres_verbose,logical(opt.gmres_verbose),true);
end
solveSeconds = toc(timer);

%% Reconstruct rigid motion from unprojected sources, then enforce moments.
U = zeros(6*problem.P,1);
if ~opt.SVDsolve, lambda_gmres = zeros(problem.sourceOffsets(end),1); end
for k = 1:problem.P
    cols = problem.sourceOffsets(k)+(1:3*problem.N_total(k));
    body = problem.body{k};
    if ~opt.SVDsolve
        rows = problem.targetOffsets(k)+(1:3*problem.M(k));
        lambda_gmres(cols) = problem.Y{k}*(problem.UU{k}*x_gmres(rows));
    end
    lambda_i = lambda_gmres(cols);
    moments = body.Kin.'*lambda_i;
    U(6*(k-1)+(1:6)) = -moments;
end
lambda = lambda_vec;
for k = 1:problem.P
    cols = problem.sourceOffsets(k)+(1:3*problem.N_total(k));
    body = problem.body{k};
    lambda_i = lambda_gmres(cols);
    % LL*lambda_i is the moment-carrying source; never form the LL matrix.
    lambda(cols) = lambda_i-body.completion*(body.Kin.'*lambda_i)+lambda_vec(cols);
end
lambda_norm = norm(lambda,inf);
details = validate_wall_images_solution(problem,lambda,U);
uerr = details.surfaceError;
details.resvec = resvec; details.real_res = real_res;
details.converged = real_res <= opt.gmres_tol;
details.linearSystemResidual = real_res;
if opt.SVDsolve
    details.converged = NaN;
    details.linearSystemResidual = info.linearSystemResidual;
    details.systemSize = info.systemSize;
end
details.momentError = norm(details.realizedForceTorque./forceUnits-Fvec,inf);
details.timing.solveSeconds = solveSeconds;
if opt.SVDsolve
    details.timing.matvecSeconds = NaN;
elseif any(uvec)
    timer = tic;
    wall_images_matvec(uvec/norm(uvec),problem);
    details.timing.matvecSeconds = toc(timer);
else
    details.timing.matvecSeconds = 0;
end
if opt.plot, details.figures = plot_wall_images_solution(problem,lambda,details); end
% Return dimensional angular velocity and total physical coefficients.
U = U./repmat([1;1;1;opt.radius;opt.radius;opt.radius],problem.P,1);
lambda = lambda.*problem.sourceUnits;
if opt.gmres_verbose
    if opt.SVDsolve
        fprintf('Mobility SVD: rank %d, rectangular residual %.3e, independent %.3e.\n', ...
            info.rank,info.linearSystemResidual,uerr);
    else
        fprintf('Mobility: %d iterations, GMRES %.3e, collocation %.3e, independent %.3e.\n', ...
            iters,real_res,details.collocationError,uerr);
    end
    if ~opt.SVDsolve && ~details.converged
        fprintf('  Recomputed GMRES residual did not reach the requested %.3e tolerance.\n', ...
            opt.gmres_tol);
    end
end
end

function [U,iters,lambda_norm,uerr,lambda,details] = self_test()
% Drive a configurable near-wall particle problem with independent loads.
%q = [0,0,1.001;2.001,0,1.001;1,4,1.002];
q = [0,0,1.001;1,4,1.002];
[rin,rout] = init_spheres(q,0.63,700,2);
opt = struct('N_image',40,'wallImageFamilies',"ST_tD", ...
    'pairImageFamilies',"SDR",'viscosity',1,'SVDsolve',false, ...
    'fmm',1,'fmm_tol',1e-12,'plot',1, ...
    'tol',1e-10,'gmres_tol',1e-9,'maxit',400);
Fvec = [1;0;-1;0;0.2;0; -1;0;0.5;0.1;0;0];% 0;1;0;0;0;0.5];
[U,iters,lambda_norm,uerr,lambda,details] = ...
    solve_mobility_wall_images(q,rin,rout,Fvec,opt);
if ~all(isfinite([U;lambda])) || details.momentError > 1e-8*max(1,norm(Fvec,inf))
    error('wall:ImagesSelfTest','Mobility coefficients are nonfinite or moments are incorrect.');
end
fprintf(['Mobility self-test: wall gap %.3g, pair gap %.3g, %d iterations.\n', ...
    '  Independent error %.3e (goal 1e-6; reported without an accuracy assertion).\n', ...
    '  Moment error %.3e.\n'], ...
    details.problem.minWallGap,details.problem.minBodyGap,iters,uerr,details.momentError);
if opt.SVDsolve
    fprintf('  Global SVD rank %d; rectangular residual %.3e; solve %.3gs.\n', ...
        details.retainedRanks,details.linearSystemResidual,details.timing.solveSeconds);
else
    fprintf('  Warm hybrid matvec %.3gs; setup %.3gs; solve %.3gs.\n', ...
        details.timing.matvecSeconds,details.timing.preparationSeconds,details.timing.solveSeconds);
end
end
