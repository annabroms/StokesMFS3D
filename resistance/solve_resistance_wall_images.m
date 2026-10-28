function [Fvec,iters,lambda_norm,err_res,lambda,details] = ...
    solve_resistance_wall_images(q,rin,rout,U,opt)
%SOLVE_RESISTANCE_WALL_IMAGES Resistance with wall ST_tD and close-pair SDR.
%   [Fvec,iters,lambda_norm,err_res,lambda,details] = ...
%       solve_resistance_wall_images(q,rin,rout,U,opt)
%   solve_resistance_wall_images() runs the plotted three-sphere self-test.
%
% q: P-by-3 physical centres of equal spheres above z=opt.wall_height.
% rin: N*P-by-3 ordinary interior proxy points, grouped by particle.
% rout: M*P-by-3 ordinary surface collocation points. Contact caps and image
%       lines are appended internally; details.rin/rout contain final grids.
% U: 6P-by-1 prescribed [vx;vy;vz;omega_x;omega_y;omega_z] per particle.
% Fvec: 6P-by-1 [Fx;Fy;Fz;Tx;Ty;Tz] ON THE FLUID; particle drag is -Fvec.
% iters: GMRES iterations; lambda_norm: dimensionless coefficient infinity
% norm; err_res: independent maximum pointwise velocity error divided by
% maximum rigid surface speed (absolute error for zero motion). NaN if the
% independent check is disabled. Missing a 1e-6 goal is reported, not fatal.
%
% lambda: total physical mixed-source coefficients, with interleaved xyz at
% each node. Default body/family blocks: [proxy S;wall S;wall T_t;wall D;
% pair 1 S;pair 1 D;pair 1 R;...]. S has force units, T_t/R force*length,
% and D force*length^2. details.lambda stores dimensionless coefficients;
% details.sourceUnits converts them to lambda. details.problem contains
% normalized geometry/factors for flow evaluation; do not pass physical
% lambda to apply_wall_images_velocity without dividing by sourceUnits.
%
% One-body SVD factors are Y{k} and UU{k}: Y{k}*(UU{k}*x) maps collocation
% data to coefficients. UU holds retained left singular vectors transposed;
% Y holds right singular vectors divided by retained singular values. Kin'
% maps ALL source families to particle force/torque; Kout maps rigid motion
% to surface velocity. Sources on pair lines belong to their owning body.
% Every source includes its complete no-slip wall correction. No slip or
% background flow is imposed; particles may approach but must not overlap.
%
% Main options retain solve_resistance names: tol (relative SVD cutoff),
% maxit, gmres_tol, gmres_verbose, fmm, fmm_tol, plot, compute_residual.
% opt.viscosity=1 by default (mu is accepted as a legacy alias).
% opt.SVDsolve=false uses one-body preconditioning and helsing_gmres.
% true instead solves the full rectangular physical matrix by truncated
% global SVD, with direct kernels and no one-body factors or FMM dependency.
% In this mode iters=0, real_res/converged are NaN (GMRES not applicable),
% and details.linearSystemResidual is the relative rectangular fit residual.
% N_image counts points per line (default 30), not proxy count N. See
% prepare_wall_images_problem for cap and physical-unit option defaults.
% Select line families independently with opt.wallImageFamilies="ST_tD"
% and opt.pairImageFamilies="SDR" (defaults). Tokens are S,R,D,T_t,T_n,T_n2;
% T includes all three stresslet orientations. t is along each wall/pair
% line; n/n2 are transverse. "none" disables one line type. Actual selected
% coefficient blocks and their ordering are in details.problem.body{k}.basis.
% fmm=true uses a HYBRID action: batched S FMM plus direct T_t/D/R. The
% latter retains quadratic cost in body count. StokesHalfspaceFMM must be
% on the path; fmm=true also requires its compiled FMM3D dependency.
%
% Example grids: [rin,rout]=init_spheres(q/a,0.63,700,2); then multiply
% both by a, or use init_spheres(q,0.63,700,2,a). Explicitly set opt.fmm=1
% to select the hybrid path when using options returned by init_spheres.

if nargin == 0
    [Fvec,iters,lambda_norm,err_res,lambda,details] = self_test();
    return
end
if nargin < 5, opt = struct(); end
validateattributes(U,{'double'},{'column','real','finite','numel',6*size(q,1)});
problem = prepare_wall_images_problem(q,rin,rout,opt,'resistance');
opt = problem.options;
% Rotation scales with length: normalized rigid data are [V;a*Omega].
U_scaled = U.*repmat([1;1;1;opt.radius;opt.radius;opt.radius],problem.P,1);
u_bndry = zeros(problem.targetOffsets(end),1);
for k = 1:problem.P
    rows = problem.targetOffsets(k)+(1:3*problem.M(k));
    u_bndry(rows) = problem.body{k}.Kout*U_scaled(6*(k-1)+(1:6));
end

%% Solve by global SVD or by the one-body-preconditioned Helsing iteration.
timer = tic;
if opt.SVDsolve
    [lambda,info] = solve_wall_images_svd(problem,u_bndry);
    iters = 0; resvec = []; real_res = NaN;
    problem.ranks = info.rank; problem.singularValues = {info.singularValues};
    problem.timing.globalAssemblySeconds = info.assemblySeconds;
    problem.timing.globalFactorizationSeconds = info.factorizationSeconds;
elseif ~any(u_bndry)
    mu_gmres = zeros(size(u_bndry)); iters = 0; resvec = 0; real_res = 0;
else
    % The callback replaces each factored self block with identity.
    matvec = @(x) wall_images_matvec(x,problem);
    [mu_gmres,iters,resvec,real_res] = helsing_gmres(matvec,u_bndry, ...
        numel(u_bndry),min(opt.maxit,numel(u_bndry)),opt.gmres_tol, ...
        0,opt.gmres_verbose,logical(opt.gmres_verbose),true);
end
solveSeconds = toc(timer);

%% Reconstruct all sources, then sum force/torque with the mixed-source map.
if ~opt.SVDsolve, lambda = zeros(problem.sourceOffsets(end),1); end
Fvec = zeros(6*problem.P,1);
for k = 1:problem.P
    cols = problem.sourceOffsets(k)+(1:3*problem.N_total(k));
    if ~opt.SVDsolve
        rows = problem.targetOffsets(k)+(1:3*problem.M(k));
        lambda(cols) = problem.Y{k}*(problem.UU{k}*mu_gmres(rows));
    end
    lambda_i = lambda(cols);
    Fvec(6*(k-1)+(1:6)) = problem.body{k}.Kin.'*lambda_i;
end
lambda_norm = norm(lambda,inf);
details = validate_wall_images_solution(problem,lambda,U_scaled);
err_res = details.surfaceError;
details.resvec = resvec; details.real_res = real_res;
details.converged = real_res <= opt.gmres_tol;
details.linearSystemResidual = real_res;
if opt.SVDsolve
    details.converged = NaN;
    details.linearSystemResidual = info.linearSystemResidual;
    details.systemSize = info.systemSize;
end
details.timing.solveSeconds = solveSeconds;
if opt.SVDsolve
    details.timing.matvecSeconds = NaN;
elseif any(u_bndry)
    timer = tic;
    wall_images_matvec(u_bndry/norm(u_bndry),problem);
    details.timing.matvecSeconds = toc(timer);
else
    details.timing.matvecSeconds = 0;
end
if opt.plot, details.figures = plot_wall_images_solution(problem,lambda,details); end
% Physical conversion happens exactly once, after normalized diagnostics.
Fvec = Fvec.*repmat([opt.mu*opt.radius*ones(3,1); ...
    opt.mu*opt.radius^2*ones(3,1)],problem.P,1);
lambda = lambda.*problem.sourceUnits;
if opt.gmres_verbose
    if opt.SVDsolve
        fprintf('Resistance SVD: rank %d, rectangular residual %.3e, independent %.3e.\n', ...
            info.rank,info.linearSystemResidual,err_res);
    else
        fprintf('Resistance: %d iterations, GMRES %.3e, collocation %.3e, independent %.3e.\n', ...
            iters,real_res,details.collocationError,err_res);
    end
    if ~opt.SVDsolve && ~details.converged
        fprintf('  Recomputed GMRES residual did not reach the requested %.3e tolerance.\n', ...
            opt.gmres_tol);
    end
end
end

function [Fvec,iters,lambda_norm,err_res,lambda,details] = self_test()
% Show simultaneous wall and pair lubrication with one-body preconditioning.
q = [0,0,1.001;2.001,0,1.001;1,4,1.002];
[rin,rout] = init_spheres(q,0.63,700,2);
opt = struct('N_image',30,'wallImageFamilies',"ST_tD", ...
    'pairImageFamilies',"SDR",'viscosity',1,'SVDsolve',true, ...
    'fmm',1,'fmm_tol',1e-12,'plot',1, ...
    'tol',1e-12,'gmres_tol',1e-9,'maxit',400);
U = [1;0;-1;0;0.2;0; -1;0;0.5;0.1;0;0; 0;0.5;0;0;0;0.2];
[Fvec,iters,lambda_norm,err_res,lambda,details] = ...
    solve_resistance_wall_images(q,rin,rout,U,opt);
if ~all(isfinite([Fvec;lambda]))
    error('wall:ImagesSelfTest','Resistance produced nonfinite coefficients or forces.');
end
fprintf(['Resistance self-test: wall gap %.3g, pair gap %.3g, %d iterations.\n', ...
    '  Independent error %.3e (goal 1e-6; reported without an accuracy assertion).\n'], ...
    details.problem.minWallGap,details.problem.minBodyGap,iters,err_res);
if opt.SVDsolve
    fprintf('  Global SVD rank %d; rectangular residual %.3e; solve %.3gs.\n', ...
        details.retainedRanks,details.linearSystemResidual,details.timing.solveSeconds);
else
    fprintf('  Warm hybrid matvec %.3gs; setup %.3gs; solve %.3gs.\n', ...
        details.timing.matvecSeconds,details.timing.preparationSeconds,details.timing.solveSeconds);
end
end
