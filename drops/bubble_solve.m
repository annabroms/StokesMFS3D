function solution = bubble_solve(opt,Fg,geometry)
%BUBBLE_SOLVE Solve force-driven exterior stagnant-cap bubble flow.
if nargin == 0
    bubble_self_test();
    solution = [];
    return
end
opt = bubble_options(opt);
if nargin<3 || isempty(geometry)
    geometry = bubble_geometry(opt);
end
rows = 3*geometry.collocationCount;
columns = geometry.coefficientCount;
thin = min(rows,columns);
inputMB = 8*rows*columns/1e6;
factorsMB = 8*(rows*thin+columns*thin+thin^2)/1e6;
fprintf(['Bubble preflight: %d distinct source sites, %d cone families, ', ...
    '%d azimuths per collocation ring; matrix %d x %d. ', ...
    'Input %.1f MB, known dense storage >= %.1f MB ', ...
    'before assembly arrays and SVD workspace.\n'], ...
    geometry.sourceCount,geometry.coneFamilyCount, ...
    geometry.ringAzimuthCount,rows,columns,inputMB, ...
    2*inputMB+factorsMB);
system = bubble_system(geometry,opt,Fg);
Uref = abs(Fg(3))/stagnant_cap_exact_drag( ...
    opt.phi,0,opt.muOut,opt.radius);
tractionScale = opt.muOut*Uref/opt.radius;
rowScale = [ones(system.normalRows+system.capTangentRows,1)/Uref; ...
    ones(system.mobileTractionRows,1)/tractionScale];
A = rowScale.*system.matrix;
b = rowScale.*system.rhs;
assert(isequal(size(A),[rows,columns]));
fprintf(['Bubble: %d base + %d cone sources = %d points, %d scalar unknowns; ', ...
    '%d cap + %d mobile collocation points.\n'], ...
    geometry.baseSourceCount,geometry.coneSourceCount, ...
    geometry.sourceCount,columns,geometry.capCount,geometry.mobileCount);
fprintf(['Boundary rows: normal %d, cap tangent %d, mobile traction %d. ', ...
    'SVD matrix %d x %d; input %.1f MB, factors %.1f MB, ', ...
    'known dense storage >= %.1f MB, excluding workspace.\n'], ...
    system.normalRows,system.capTangentRows, ...
    system.mobileTractionRows,rows,columns,inputMB,factorsMB, ...
    2*inputMB+factorsMB);
[rightOverSingular,left,singularValues] = ...
    getPseudoFactors(A,opt.svdTolerance,0);
assert(~isempty(singularValues),'SVD discarded every singular value.');
alpha = rightOverSingular*(left.'*b);
lambda = alpha-system.forceCompletion*(system.forceMap*alpha) ...
    +system.lambda0;
U = -opt.completionScale*system.forceMap*alpha;
FH = -system.forceMap*lambda;
torque = -system.momentMap(4:6,:)*lambda;
solution.lambda = lambda;
solution.alpha = alpha;
solution.U = U;
solution.FH = FH;
solution.Fg = Fg(:);
solution.torque = torque;
solution.forceResidual = norm(FH+Fg(:))/norm(Fg);
solution.transverseDrift = norm(U(1:2))/max(abs(U(3)),eps);
solution.basis = system.basis;
solution.grid = geometry;
solution.Uref = Uref;
solution.tractionScale = tractionScale;
solution.matrixSize = [rows,columns];
solution.retainedRank = numel(singularValues);
solution.svdTolerance = opt.svdTolerance;
solution.svdInputMB = inputMB;
solution.svdKnownDenseMB = 2*inputMB+factorsMB;
solution.fitRMS = norm(A*alpha-b)/sqrt(rows);
solution.fitMaximum = norm(A*alpha-b,inf);
fprintf('Retained numerical rank %d/%d; force residual %.3g.\n', ...
    solution.retainedRank,min(rows,columns),solution.forceResidual);
end

function bubble_self_test()
% Check force balance, finite axial drift, endpoints, and a cone-source solve.
for phi = [0,pi/3,pi]
    o = bubble_options(struct('sourceRequest',50,'phi',phi, ...
        'ringLevels',2,'svdTolerance',1e-11));
    g = bubble_geometry(o);
    s = bubble_solve(o,[0;0;-1],g);
    assert(s.forceResidual<1e-10 && all(isfinite(s.U)));
    assert(isequal(s.matrixSize,[3*g.collocationCount, ...
        g.coefficientCount]));
    referenceUz = -1/stagnant_cap_exact_drag(phi,0,o.muOut,o.radius);
    speedError = abs(s.U(3)-referenceUz)/abs(referenceUz);
    if phi==0 || phi==pi
        assert(speedError<5e-3);
    else
        assert(speedError<0.1);
    end
end
o = bubble_options(struct('sourceRequest',50,'phi',pi/3, ...
    'n_level',1,'n_az_source',6,'n_az_colloc',8));
s = bubble_solve(o,[0;0;-1]);
assert(numel(s.lambda)==s.grid.coefficientCount && ...
    s.grid.coneSourceCount==6 && s.forceResidual<1e-10);
mixed = bubble_options(struct('sourceRequest',50,'phi',pi/3, ...
    'n_level',1,'n_az_source',6,'n_az_colloc',8, ...
    'refinementFamilies',["S","T_t","T_n","T_n2"]));
mixedSolution = bubble_solve(mixed,[0;0;-1]);
assert(mixedSolution.forceResidual<1e-10 && ...
    numel(mixedSolution.lambda)==mixedSolution.grid.coefficientCount && ...
    numel(mixedSolution.basis)==5);
fprintf('bubble_solve self-test passed.\n');
end
