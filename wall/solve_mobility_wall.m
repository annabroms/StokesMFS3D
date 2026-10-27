function [U,it,lnorm,berr,lambda,details] = solve_mobility_wall(q,rin,rout,FT,opt)
%SOLVE_MOBILITY_WALL Solve zero-slip mobility above a no-slip wall.
% FT is force/torque on the fluid; U=[V;Omega] for each particle.
if nargin == 0
    self_test(); U=[]; it=[]; lnorm=[]; berr=[]; lambda=[]; details=[];
    return
end

opt.kernel = "blake";
opt = wall_solver_options(opt);
validateattributes(FT,{'double'},{'real','finite','column','numel',6*size(q,1)});

setupTimer = tic;
cacheHit = ~isempty(opt.preparedProblem);
if cacheHit
    problem = opt.preparedProblem;
    opt.preparedProblem = [];
    assert(isequaln(opt,problem.options) && isequal(q,problem.physical.q) && ...
        isequal(rin,problem.physical.rin) && isequal(rout,problem.physical.rout), ...
        'wall:PreparedMismatch','Prepared geometry or options do not match.');
else
    problem = prepare_wall_problem(q,rin,rout,opt);
end
preparationCallSeconds = toc(setupTimer);
P = problem.P;
N = problem.N;
M = problem.M;

timer = tic;
force = FT./problem.forceUnits;
lambda0 = zeros(3*N*P,1);
for p = 1:P
    ns = (p-1)*3*N+(1:3*N);
    rb = (p-1)*6+(1:6);
    lambda0(ns) = problem.body{p}.C*force(rb);
end
rhs = -apply_stokes_kernel(lambda0,problem.kernel);
rhsSeconds = toc(timer);

timer = tic;
if opt.solverMode == "denseSVD"
    factor = problem.dense.mobility;
    z = factor.Y*(factor.U.'*rhs);
    it = 0;
    residualHistory = [];
    trueResidual = NaN;
    gmresFlag = 0;
    ranks = factor.rank;
    singularValues = {factor.singularValues};
else
    systemOperator = @(x) stokes_mfs_kernel_matvec(x,problem,'mobility');
    % helsing_gmres applies f(x)-x internally and restores the identity in H.
    gmresOperator = systemOperator;
    if norm(rhs) == 0
        x = zeros(size(rhs));
        it = 0;
        residualHistory = 0;
        trueResidual = 0;
    else
        limit = min(opt.maxIterations,numel(rhs));
        [x,it,residualHistory] = helsing_gmres(gmresOperator,rhs, ...
            numel(rhs),limit,opt.gmresTolerance,0,0,true);
        trueResidual = norm(systemOperator(x)-rhs)/norm(rhs);
    end
    gmresFlag = double(trueResidual > opt.gmresTolerance);
    z = zeros(3*N*P,1);
    ranks = zeros(P,1);
    singularValues = cell(P,1);
    for p = 1:P
        ns = (p-1)*3*N+(1:3*N);
        ms = (p-1)*3*M+(1:3*M);
        factor = problem.body{p}.mobility;
        z(ns) = factor.Y*(factor.U.'*x(ms));
        ranks(p) = factor.rank;
        singularValues{p} = factor.singularValues;
    end
end
solveSeconds = toc(timer);

timer = tic;
velocityInternal = zeros(6*P,1);
physicalSource = lambda0;
for p = 1:P
    ns = (p-1)*3*N+(1:3*N);
    rb = (p-1)*6+(1:6);
    body = problem.body{p};
    velocityInternal(rb) = -body.H*z(ns);
    physicalSource(ns) = physicalSource(ns)+z(ns)-body.C*(body.H*z(ns));
end
U = velocityInternal.*problem.velocityUnits;
lambda = physicalSource*problem.strengthUnit;
lnorm = norm(lambda,inf);
reconstructionSeconds = toc(timer);

timer = tic;
physicalVelocity = apply_stokes_kernel(physicalSource,problem.kernel);
desiredVelocity = zeros(3*M*P,1);
FTrealized = zeros(6*P,1);
for p = 1:P
    ms = (p-1)*3*M+(1:3*M);
    rb = (p-1)*6+(1:6);
    body = problem.body{p};
    desiredVelocity(ms) = body.K*velocityInternal(rb);
    ns = (p-1)*3*N+(1:3*N);
    FTrealized(rb) = body.H*physicalSource(ns);
end
collocationAbsolute = norm(physicalVelocity-desiredVelocity);
collocationResidual = collocationAbsolute/max(norm(desiredVelocity),realmin);
validation = struct();
berr = NaN;
if opt.computeResidual
    validation = wall_surface_validation(problem,lambda,U);
    berr = validation.error;
end
validationSeconds = toc(timer);

timing = problem.timing;
timing.cachedPreparation = cacheHit;
timing.preparationCallSeconds = preparationCallSeconds;
timing.rhsSeconds = rhsSeconds;
timing.solveSeconds = solveSeconds;
timing.reconstructionSeconds = reconstructionSeconds;
timing.validationSeconds = validationSeconds;
details = struct('options',opt,'mode','mobility','N',N,'M',M,'P',P, ...
    'converged',gmresFlag == 0,'gmresFlag',gmresFlag, ...
    'residualHistory',residualHistory,'trueGmresResidual',trueResidual, ...
    'linearSystemResidual',collocationResidual, ...
    'collocationAbsoluteResidual',collocationAbsolute,'retainedRanks',ranks, ...
    'singularValues',{singularValues},'realizedForceTorque',FTrealized, ...
    'rigidVelocity',U,'validation',validation,'timing',timing, ...
    'minBodyGap',problem.minBodyGap,'minWallGap',problem.minWallGap, ...
    'strengthUnit',problem.strengthUnit,'velocityUnits',problem.velocityUnits, ...
    'forceUnits',problem.forceUnits);
end

function self_test()
% Verify prescribed moments of the reconstructed physical source vector.
a = 1;
q = [0,0,2];
[s,t,o] = init_spheres_wall(q,a, ...
    struct('N',300,'backend',"blake",'computeResidual',false));
f = (1:6).'/6;
[u,it,~,~,lambda,d] = solve_mobility_wall(q,s,t,f,o);
assert(all(isfinite(u)) && norm(getKmat(s,q).'*lambda-f)<1e-11);
fprintf('solve_mobility_wall: PASS; %d iterations, force/torque constraint error %.3e, collocation residual %.3e.\n', ...
    it,norm(getKmat(s,q).'*lambda-f),d.linearSystemResidual);
fprintf('  Translation: [% .4g % .4g % .4g]; angular velocity: [% .4g % .4g % .4g].\n',u);
end
