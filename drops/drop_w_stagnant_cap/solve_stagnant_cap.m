function solution = solve_stagnant_cap(opt,drive,geometry)
%SOLVE_STAGNANT_CAP Solve prescribed-force or prescribed-speed drop flow.
if nargin == 0
    self_test();
    return
end
opt = stagnant_cap_options(opt);
if nargin<3 || isempty(geometry)
    geometry = stagnant_cap_geometry(opt);
end
system = stagnant_cap_system(geometry,opt,drive);
M = system.M;
rowScale = [ones(4*M,1)/system.velocityScale; ...
    ones(2*M,1)/system.tractionScale];
A = rowScale.*system.matrix;
b = rowScale.*system.rhs;
% Real-double storage for this economy SVD. Both system.matrix and A are
% live here; LAPACK may allocate additional workspace beyond this bound.
[svdRows,svdCols] = size(A);
svdThin = min(svdRows,svdCols);
svdInputMB = 8*svdRows*svdCols/1e6;
svdFactorsMB = 8*(svdRows*svdThin + ...
    svdCols*svdThin + svdThin^2)/1e6;
svdKnownDenseMB = 2*svdInputMB + svdFactorsMB;
if isfield(geometry,'exteriorBasis')
    exteriorFamilies = unique([geometry.exteriorBasis.family]);
    interiorFamilies = unique([geometry.interiorBasis.family]);
    fprintf(['Source points inside/outside %d/%d (%d distinct); ', ...
        'families %d/%d; basis blocks %d/%d; scalar unknowns %d.\n'], ...
        size(geometry.innerSources,1),size(geometry.outerSources,1), ...
        size(geometry.innerSources,1)+size(geometry.outerSources,1), ...
        numel(exteriorFamilies),numel(interiorFamilies), ...
        numel(geometry.exteriorBasis),numel(geometry.interiorBasis), ...
        svdCols);
end
fprintf(['Boundary blocks: continuity %d, normal %d, cap tangent %d, ', ...
    'mobile traction %d rows; exterior/interior coefficient columns ', ...
    '%d/%d.\n'],3*M,M,2*system.Mc,2*system.Mm, ...
    system.outerColumns,system.innerColumns);
fprintf(['SVD %s: %d rows x %d cols; input matrix %.1f MB; ', ...
    'economy U/S/V %.1f MB; known live dense storage >= %.1f MB ', ...
    '(real double, 1 MB=10^6 bytes; excludes SVD workspace).\n'], ...
    system.kind,svdRows,svdCols,svdInputMB,svdFactorsMB, ...
    svdKnownDenseMB);
[rightOverSingular,left,singularValues] = ...
    getPseudoFactors(A,opt.svdTolerance,0);
assert(~isempty(singularValues),'SVD discarded every singular value.');
coefficients = rightOverSingular*(left.'*b);
fitResidual = A*coefficients-b;
outerColumns = system.outerColumns;
if system.kind=="force"
    alpha = coefficients(1:outerColumns);
    lambda = alpha-system.forceCompletion*(system.forceMap*alpha) ...
        +system.lambda0;
    translation = -opt.completionScale*system.forceMap*alpha;
else
    lambda = coefficients(1:outerColumns);
    translation = [0;0;drive.Uz];
end
gamma = coefficients(outerColumns+1:end);
force = -system.forceMap*lambda;
torque = -system.momentMap(4:6,:)*lambda;
solution.kind = system.kind;
solution.lambda = lambda;
solution.gamma = gamma;
solution.exteriorBasis = system.exteriorBasis;
solution.interiorBasis = system.interiorBasis;
solution.U = translation;
solution.FH = force;
solution.torque = torque;
solution.sourceForce = -force;
solution.forceResidual = [];
if system.kind=="force"
    solution.Fg = system.force;
    solution.forceResidual = norm(force+system.force);
else
    solution.Fg = [];
end
solution.velocityScale = system.velocityScale;
solution.tractionScale = system.tractionScale;
solution.svdTolerance = opt.svdTolerance;
solution.retainedRank = numel(singularValues);
solution.matrixSize = size(A);
solution.svdInputMB = svdInputMB;
solution.svdFactorsMB = svdFactorsMB;
solution.svdKnownDenseMB = svdKnownDenseMB;
solution.fitRMS = norm(fitResidual)/sqrt(numel(fitResidual));
solution.fitMaximum = norm(fitResidual,inf);
solution.grid = geometry;
end

function self_test()
%SELF_TEST Check force balance, resistance reciprocity, and zero transverse drift.
o = stagnant_cap_options(struct('sourceRequest',50,'ringLevels',2, ...
    'svdTolerance',1e-11));
g = stagnant_cap_geometry(o);
res = solve_stagnant_cap(o,struct('kind','speed','Uz',-1),g);
mob = solve_stagnant_cap(o,struct('kind','force', ...
    'Fg',[0 0 res.sourceForce(3)]),g);
assert(res.sourceForce(3)<0);
assert(mob.forceResidual<1e-10*norm(mob.Fg));
assert(norm(mob.U-res.U)<0.2*norm(res.U));
assert(norm(mob.U(1:2))<0.1*norm(mob.U));
assert(isequal(res.matrixSize,[6*g.collocationCount,6*g.N]));
assert(abs(res.svdInputMB-8*prod(res.matrixSize)/1e6)<eps);
assert(res.svdKnownDenseMB>2*res.svdInputMB);
o.n_level = 1;
o.n_az_source = 6;
o.n_az_colloc = 8;
refined = stagnant_cap_refined_geometry(o);
driven = solve_stagnant_cap(o,struct('kind','force', ...
    'Fg',[0 0 -1]),refined);
assert(isequal(driven.matrixSize,[6*refined.collocationCount, ...
    refined.exteriorCoefficientCount+refined.interiorCoefficientCount]));
assert(driven.forceResidual<1e-10 && all(isfinite(driven.U)));
assert(isfinite(norm(driven.torque)));
% The three explicit stresslet orientations must pass through a force solve
% without the legacy T block or a coefficient-order mismatch.
o.refinementFamilies = "ST_tT_nT_n2";
o.n_az_source = 4;
o.n_az_colloc = 8;
oriented = stagnant_cap_refined_geometry(o);
orientedSolve = solve_stagnant_cap(o,struct('kind','force', ...
    'Fg',[0 0 -1]),oriented);
assert(numel(orientedSolve.exteriorBasis)==5 && ...
    all(isfinite(orientedSolve.U)) && ...
    orientedSolve.forceResidual<1e-10);
o.refinementFamilies = "S";
o.derivative_direction = "radial";
derivativeGeometry = stagnant_cap_refined_geometry(o);
derivativeSolve = solve_stagnant_cap(o,struct('kind','force', ...
    'Fg',[0 0 -1]),derivativeGeometry);
assert(derivativeSolve.forceResidual<1e-10 && ...
    all(isfinite(derivativeSolve.U)) && ...
    isfinite(norm(derivativeSolve.torque)));
fprintf('solve_stagnant_cap self-test passed.\n');
end
