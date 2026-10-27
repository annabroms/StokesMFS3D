function [FT,it,lnorm,berr,lambda,details] = solve_resistance_wall(q,rin,rout,U,opt)
%SOLVE_RESISTANCE_WALL Solve resistance above a no-slip wall with matrix-free GMRES.
% FT is force/torque on the fluid; lambda contains physical point-force strengths.
if nargin == 0
    self_test(); FT=[]; it=[]; lnorm=[]; berr=[]; lambda=[]; details=[];
    return
end

opt.kernel = "blake";
opt = wall_solver_options(opt);
validateattributes(U,{'double'},{'real','finite','column','numel',6*size(q,1)});

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
prescribed = U./problem.velocityUnits;
rhs = zeros(3*M*P,1);
for p = 1:P
    ms = (p-1)*3*M+(1:3*M);
    rb = (p-1)*6+(1:6);
    rhs(ms) = problem.body{p}.K*prescribed(rb);
end
rhsSeconds = toc(timer);

timer = tic;
if opt.solverMode == "denseSVD"
    factor = problem.dense.resistance;
    z = factor.Y*(factor.U.'*rhs);
    it = 0;
    residualHistory = [];
    trueResidual = NaN;
    gmresFlag = 0;
    ranks = factor.rank;
    singularValues = {factor.singularValues};
else
    systemOperator = @(x) stokes_mfs_kernel_matvec(x,problem,'resistance');
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
        factor = problem.body{p}.resistance;
        z(ns) = factor.Y*(factor.U.'*x(ms));
        ranks(p) = factor.rank;
        singularValues{p} = factor.singularValues;
    end
end
solveSeconds = toc(timer);

timer = tic;
lambda = z*problem.strengthUnit;
FT = zeros(6*P,1);
for p = 1:P
    ns = (p-1)*3*N+(1:3*N);
    rb = (p-1)*6+(1:6);
    FT(rb) = problem.body{p}.H*z(ns);
end
FT = FT.*problem.forceUnits;
lnorm = norm(lambda,inf);
reconstructionSeconds = toc(timer);

timer = tic;
physicalVelocity = apply_stokes_kernel(z,problem.kernel);
collocationAbsolute = norm(physicalVelocity-rhs);
collocationResidual = collocationAbsolute/max(norm(rhs),realmin);
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
details = struct('options',opt,'mode','resistance','N',N,'M',M,'P',P, ...
    'converged',gmresFlag == 0,'gmresFlag',gmresFlag, ...
    'residualHistory',residualHistory,'trueGmresResidual',trueResidual, ...
    'linearSystemResidual',collocationResidual, ...
    'collocationAbsoluteResidual',collocationAbsolute,'retainedRanks',ranks, ...
    'singularValues',{singularValues},'realizedForceTorque',FT, ...
    'rigidVelocity',U,'validation',validation,'timing',timing, ...
    'minBodyGap',problem.minBodyGap,'minWallGap',problem.minWallGap, ...
    'strengthUnit',problem.strengthUnit,'velocityUnits',problem.velocityUnits, ...
    'forceUnits',problem.forceUnits);
end

function self_test()
% Exercise one resistance and one mobility solve, then display their geometry.
q = [0,0,1.5;3,0,1.5];
[s,t,o] = init_spheres_wall(q,1, ...
    struct('N',400,'backend',"blake",'validationResolution',[24,16]));
problem = prepare_wall_problem(q,s,t,o);
o.preparedProblem = problem;
U = [(1:6).'/6;-(1:6).'/6];
FT = [(1:6).'/6;(1:6).'/12];
[~,itR,~,~,lambda,dR] = solve_resistance_wall(q,s,t,U,o);
[u,itM,~,~,lambdaM,dM] = solve_mobility_wall(q,s,t,FT,o);
resistanceValidation = wall_surface_validation(problem,lambda,U);
mobilityValidation = wall_surface_validation(problem,lambdaM,u);
surfaceValidation = resistanceValidation;
surfaceValidation.pointErrors = max([resistanceValidation.pointErrors, ...
    mobilityValidation.pointErrors],[],2);
surfaceValidation.error = max(surfaceValidation.pointErrors);
fprintf('solve_resistance_wall self-test: resistance %d iterations, mobility %d iterations.\n', ...
    itR,itM);
fprintf('  Collocation residuals: resistance %.3e, mobility %.3e.\n', ...
    dR.linearSystemResidual,dM.linearSystemResidual);
plot_wall_geometry(q,1,s,struct('wall_height',0,'delta_body',1, ...
    'delta_wall',0.5,'visible',get(groot,'DefaultFigureVisible')),surfaceValidation);
end
