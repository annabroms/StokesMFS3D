function solution = bubble_explicit_solve(opt,Fg,geometry,settings)
%BUBBLE_EXPLICIT_SOLVE Scaled least squares for independent physical lambda,U.
% Force rows are weighted observations: force balance is measured, not repaired.
if nargin==0
    self_test();
    solution = [];
    return
end
opt = bubble_options(opt);
if nargin<3 || isempty(geometry), geometry = bubble_geometry(opt); end
if nargin<4, settings = struct(); end
assert(isstruct(settings) && isscalar(settings) && ...
    all(ismember(fieldnames(settings),{'forceWeight'})), ...
    'bubble_explicit_solve:BadSettings','Only forceWeight is supported.');
if ~isfield(settings,'forceWeight'), settings.forceWeight = 1; end
validateattributes(settings.forceWeight,{'numeric'}, ...
    {'scalar','real','finite','positive'});
assert(~opt.useAreaWeights && ~opt.useColumnScaling && ...
    ~opt.useConditionCountScaling && isequal(opt.conditionWeights,[1 1 1]), ...
    'bubble_explicit_solve:UnsupportedWeighting', ...
    'Explicit solve currently supports isotropic unit condition weights without area/column/count scaling.');
Fg = Fg(:);
assert(numel(Fg)==3 && all(isfinite(Fg)) && norm(Fg)>0);
rows = 3*geometry.collocationCount+3;
columns = geometry.coefficientCount+3;
memory = bubble_svd_memory(rows,columns);
assert(memory.estimatedGB<=opt.maxSvdMemoryGB, ...
    'bubble_explicit_solve:MemoryBudget', ...
    'Estimated dense solve %.3g GB exceeds budget %.3g GB.', ...
    memory.estimatedGB,opt.maxSvdMemoryGB);
fprintf(['Explicit bubble preflight: %d base + %d cone Stokeslets; ', ...
    '%d cap + %d mobile targets; matrix %d x %d; estimated %.3g GB.\n'], ...
    geometry.baseSourceCount,geometry.coneSourceCount, ...
    geometry.capCount,geometry.mobileCount,rows,columns,memory.estimatedGB);
system = bubble_explicit_system(geometry,opt,Fg);
Uref = norm(Fg)/stagnant_cap_exact_drag(opt.phi,0,opt.muOut,opt.radius);
tractionScale = opt.muOut*Uref/opt.radius;
forceScale = norm(Fg);
rowScale = [ones(3*geometry.capCount+geometry.mobileCount,1)/Uref; ...
    ones(2*geometry.mobileCount,1)/tractionScale; ...
    ones(3,1)*settings.forceWeight/forceScale];
physicalRowScale = rowScale;
physicalRowScale(system.forceRows) = 1/forceScale;
A = rowScale.*system.matrix;
b = rowScale.*system.rhs;
fprintf(['Explicit SVD matrix %d x %d; rows cap/mobile-normal/', ...
    'mobile-traction/force = %d/%d/%d/%d; cutoff %.3g; force weight %.3g.\n'], ...
    rows,columns,system.blockCounts,opt.svdTolerance,settings.forceWeight);
[rightOverSingular,left,singularValues] = getPseudoFactors(A,opt.svdTolerance,0);
assert(~isempty(singularValues),'SVD discarded every singular value.');
x = rightOverSingular*(left.'*b);
lambda = x(1:end-3);
U = x(end-2:end);
FH = -system.forceMap*lambda;
forceDefect = FH+Fg;
fit = system.matrix*x-system.rhs;
solution.lambda = lambda;
solution.U = U;
solution.FH = FH;
solution.Fg = Fg;
solution.torque = -system.momentMap(4:6,:)*lambda;
solution.forceDefect = forceDefect;
solution.forceResidual = norm(forceDefect)/forceScale;
solution.transverseDrift = norm(U(1:2))/max(abs(U(3)),eps);
solution.basis = system.basis;
solution.grid = geometry;
solution.Uref = Uref;
solution.tractionScale = tractionScale;
solution.forceScale = forceScale;
solution.matrixSize = [rows,columns];
solution.retainedRank = numel(singularValues);
solution.svdTolerance = opt.svdTolerance;
solution.sigmaMax = singularValues(1);
solution.svdAbsoluteThreshold = opt.svdTolerance*singularValues(1);
solution.smallestRetainedSingularValue = singularValues(end);
solution.svdInputMB = 1000*memory.inputGB;
solution.svdKnownDenseMB = 1000*memory.knownGB;
solution.svdMemoryEstimateGB = memory.estimatedGB;
solution.conditionWeights = opt.conditionWeights;
solution.useConditionCountScaling = false;
solution.useColumnScaling = false;
solution.useAreaWeights = false;
solution.rowScale = rowScale;
solution.columnScale = ones(columns,1);
solution.settings = settings;
solution.fitRMS = norm(physicalRowScale.*fit)/sqrt(rows);
solution.fitMaximum = norm(physicalRowScale.*fit,inf);
solution.weightedFitRMS = norm(rowScale.*fit)/sqrt(rows);
solution.weightedFitMaximum = norm(rowScale.*fit,inf);
solution.blockCounts = system.blockCounts;
fprintf(['Explicit sigmaMax %.6g; absolute cutoff %.6g; retained rank %d/%d. ', ...
    'Relative force residual %.6g; FH+Fg = [% .4g % .4g % .4g].\n'], ...
    solution.sigmaMax,solution.svdAbsoluteThreshold,solution.retainedRank, ...
    min(rows,columns),solution.forceResidual,forceDefect);
end

function self_test()
% Compare a small direct pseudoinverse, physical force defect and drag limits.
o = bubble_options(struct('sourceRequest',12,'phi',pi/3, ...
    'n_level',1,'n_az_source',6,'n_az_colloc',8,'svdTolerance',1e-10));
g = bubble_geometry(o);
s = bubble_explicit_solve(o,[0;0;-1],g,struct('forceWeight',2));
a = bubble_explicit_system(g,o,s.Fg);
A = s.rowScale.*a.matrix;
fprintf('Self-test reference pseudoinverse: matrix %d x %d; cutoff %.3g.\n', ...
    size(A,1),size(A,2),s.svdAbsoluteThreshold);
x = pinv(A,s.svdAbsoluteThreshold)*(s.rowScale.*a.rhs);
assert(norm([s.lambda;s.U]-x)<1e-9*max(1,norm(x)));
assert(abs(s.forceResidual-norm(a.forceMap*s.lambda-s.Fg))<1e-13);
trial = zeros(size(s.lambda));
assert(norm(a.forceMap*trial-s.Fg)/norm(s.Fg)==1);
for phi = [0 pi]
    endpoint = bubble_options(struct('sourceRequest',50,'phi',phi, ...
        'svdTolerance',1e-11,'muOut',1.7,'radius',1.2));
    e = bubble_explicit_solve(endpoint,[0;0;-1]);
    theory = -1/stagnant_cap_exact_drag(phi,0,endpoint.muOut,endpoint.radius);
    assert(abs(e.U(3)/theory-1)<5e-3 && e.forceResidual<5e-3);
end
o.useConditionCountScaling = true;
try
    bubble_explicit_solve(o,[0;0;-1],g);
    error('bubble_explicit_solve:SelfTest','Unsupported weighting accepted.');
catch ME
    assert(strcmp(ME.identifier,'bubble_explicit_solve:UnsupportedWeighting'));
end
fprintf('bubble_explicit_solve self-test passed.\n');
end
