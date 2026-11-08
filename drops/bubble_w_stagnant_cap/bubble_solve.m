function solution = bubble_solve(opt,Fg,geometry,pointWeights,capMode)
%BUBBLE_SOLVE Solve force-driven exterior stagnant-cap bubble flow.
% useAreaWeights applies sqrt(local area/mean area) to each boundary row.
% useColumnScaling normalizes the resulting columns to unit Euclidean norm.
% Recover physical coefficients before the force and translation maps.
% The relative SVD cutoff is applied to the scaled matrix.
% Optional capMode="cartesian" changes cap rows, retaining force completion.
if nargin == 0
    bubble_self_test();
    solution = [];
    return
end
if nargin<5, capMode = "split"; end
capMode = string(capMode);
assert(isscalar(capMode) && ismember(capMode,["split","cartesian"]));
opt = bubble_options(opt);
if nargin<3 || isempty(geometry)
    geometry = bubble_geometry(opt);
end
rows = 3*geometry.collocationCount;
columns = geometry.coefficientCount;
memoryEstimate = bubble_svd_memory(rows,columns);
assert(memoryEstimate.estimatedGB<=opt.maxSvdMemoryGB, ...
    'bubble_solve:MemoryBudget', ...
    'Estimated dense solve %.3g GB exceeds budget %.3g GB; reduce points or sources.', ...
    memoryEstimate.estimatedGB,opt.maxSvdMemoryGB);
if nargin<4 || isempty(pointWeights)
    pointWeights = ones(geometry.collocationCount,3);
end
validateattributes(pointWeights,{'numeric'}, ...
    {'size',[geometry.collocationCount,3],'real','finite','positive'});
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
system = bubble_system(geometry,opt,Fg,capMode);
Uref = abs(Fg(3))/stagnant_cap_exact_drag( ...
    opt.phi,0,opt.muOut,opt.radius);
tractionScale = opt.muOut*Uref/opt.radius;
rowScale = [ones(system.capVelocityRows+system.normalRows+system.capTangentRows,1)/Uref; ...
    ones(system.mobileTractionRows,1)/tractionScale];
physicalRowScale = rowScale;
% Columns of pointWeights correspond to normal, cap velocity, mobile traction.
conditionScale = sqrt(pointWeights).*opt.conditionWeights;
if opt.useConditionCountScaling
    counts = [geometry.collocationCount,geometry.capCount,geometry.mobileCount];
    conditionScale = conditionScale./sqrt(max(counts,1));
end
if capMode=="split"
    conditionRowScale = [conditionScale(:,1); ...
        repelem(conditionScale(geometry.cap,2),2); ...
        repelem(conditionScale(geometry.mobile,3),2)];
else
    capNormal = conditionScale(geometry.cap,1);
    capTangent = conditionScale(geometry.cap,2);
    assert(all(abs(capNormal-capTangent)<=1e-13*max(capNormal,capTangent)), ...
        'bubble_solve:AnisotropicCapWeights', ...
        'Cartesian cap mode needs equal effective normal and tangential weights at each cap point.');
    conditionRowScale = [repelem(capNormal,3); ...
        conditionScale(geometry.mobile,1); ...
        repelem(conditionScale(geometry.mobile,3),2)];
end
rowScale = rowScale.*conditionRowScale;
collocationAreas = [];
if opt.useAreaWeights
    collocationAreas = collocation_areas(geometry.normals,opt.radius);
    pointWeights = sqrt(collocationAreas/mean(collocationAreas));
    if capMode=="split"
        % Legacy rows: all normals, cap tangents, then mobile traction.
        areaRowWeights = [pointWeights;repelem(pointWeights(geometry.cap),2); ...
            repelem(pointWeights(geometry.mobile),2)];
    else
        areaRowWeights = [repelem(pointWeights(geometry.cap),3); ...
            pointWeights(geometry.mobile);repelem(pointWeights(geometry.mobile),2)];
    end
    rowScale = rowScale.*areaRowWeights;
end
A = rowScale.*system.matrix;
b = rowScale.*system.rhs;
columnScale = ones(columns,1);
if opt.useColumnScaling
    columnNorms = vecnorm(A,2,1).';
    nonzero = columnNorms>0;
    columnScale(nonzero) = 1./columnNorms(nonzero);
    A = A.*columnScale.';
end
assert(isequal(size(A),[rows,columns]));
fprintf('Column scaling %d; local-area weighting %d; relative SVD cutoff %.3g.\n', ...
    opt.useColumnScaling,opt.useAreaWeights,opt.svdTolerance);
fprintf(['Bubble: %d base + %d cone sources = %d points, %d scalar unknowns; ', ...
    '%d cap + %d mobile collocation points.\n'], ...
    geometry.baseSourceCount,geometry.coneSourceCount, ...
    geometry.sourceCount,columns,geometry.capCount,geometry.mobileCount);
fprintf(['Boundary mode %s; row blocks %d/%d/%d. ', ...
    'SVD matrix %d x %d; input %.1f MB, factors %.1f MB, ', ...
    'known dense storage >= %.1f MB, excluding workspace.\n'], ...
    capMode,system.blockCounts,rows,columns,inputMB,factorsMB, ...
    2*inputMB+factorsMB);
[rightOverSingular,left,singularValues] = ...
    getPseudoFactors(A,opt.svdTolerance,0);
assert(~isempty(singularValues),'SVD discarded every singular value.');
scaledAlpha = rightOverSingular*(left.'*b);
alpha = columnScale.*scaledAlpha;
lambdaFree = alpha-system.forceCompletion*(system.forceMap*alpha);
lambda = lambdaFree+system.lambda0;
U = -opt.completionScale*system.forceMap*alpha;
FH = -system.forceMap*lambda;
torque = -system.momentMap(4:6,:)*lambda;
solution.lambda = lambda;
solution.lambdaFree = lambdaFree;
solution.lambda0 = system.lambda0;
solution.capMode = capMode;
solution.blockCounts = system.blockCounts;
solution.alpha = alpha;
solution.U = U;
solution.FH = FH;
solution.Fg = Fg(:);
solution.torque = torque;
solution.forceResidual = norm(FH+Fg(:))/norm(Fg);
solution.forceDefect = FH+Fg(:);
solution.forceFreeResidual = norm(system.forceMap*lambdaFree)/norm(Fg);
solution.completionForceResidual = norm(system.forceMap*system.lambda0-Fg(:))/norm(Fg);
solution.transverseDrift = norm(U(1:2))/max(abs(U(3)),eps);
solution.basis = system.basis;
solution.grid = geometry;
solution.Uref = Uref;
solution.tractionScale = tractionScale;
solution.matrixSize = [rows,columns];
solution.retainedRank = numel(singularValues);
solution.svdTolerance = opt.svdTolerance;
solution.sigmaMax = singularValues(1);
solution.svdAbsoluteThreshold = opt.svdTolerance*singularValues(1);
solution.smallestRetainedSingularValue = singularValues(end);
solution.svdInputMB = inputMB;
solution.svdKnownDenseMB = 2*inputMB+factorsMB;
solution.svdMemoryEstimateGB = memoryEstimate.estimatedGB;
solution.conditionWeights = opt.conditionWeights;
solution.useConditionCountScaling = opt.useConditionCountScaling;
solution.useColumnScaling = opt.useColumnScaling;
solution.useAreaWeights = opt.useAreaWeights;
solution.columnScale = columnScale;
solution.rowScale = rowScale;
solution.collocationAreas = collocationAreas;
fitResidual = system.matrix*alpha-system.rhs;
solution.fitRMS = norm(physicalRowScale.*fitResidual)/sqrt(rows);
solution.fitMaximum = norm(physicalRowScale.*fitResidual,inf);
solution.weightedFitRMS = norm(rowScale.*fitResidual)/sqrt(rows);
solution.weightedFitMaximum = norm(rowScale.*fitResidual,inf);
fprintf('sigmaMax %.6g; absolute SVD cutoff %.6g.\n', ...
    solution.sigmaMax,solution.svdAbsoluteThreshold);
fprintf('Retained numerical rank %d/%d; force residual %.3g.\n', ...
    solution.retainedRank,min(rows,columns),solution.forceResidual);
fprintf('Force-free residual %.3g; completion force residual %.3g.\n', ...
    solution.forceFreeResidual,solution.completionForceResidual);
end

function areas = collocation_areas(unit,radius)
% Lump spherical triangle areas onto nodes of the combined collocation cloud.
unit = unit./vecnorm(unit,2,2);
triangles = convhull(unit);
u = unit(triangles(:,1),:);
v = unit(triangles(:,2),:);
w = unit(triangles(:,3),:);
numerator = abs(sum(u.*cross(v,w,2),2));
denominator = 1+sum(u.*v,2)+sum(v.*w,2)+sum(w.*u,2);
triangleAreas = 2*radius^2*atan2(numerator,denominator);
areas = accumarray(triangles(:),repmat(triangleAreas/3,3,1), ...
    [size(unit,1),1]);
assert(all(isfinite(areas)) && all(areas>0) && ...
    abs(sum(areas)/(4*pi*radius^2)-1)<1e-10, ...
    'bubble_solve:BadAreas','Collocation areas must cover the whole sphere.');
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
% A regular octahedron has six equal nodal areas for every sphere radius.
octahedron = [eye(3);-eye(3)];
areas = collocation_areas(octahedron,2);
assert(max(abs(areas-16*pi/6))<1e-12);
% Closely graded rings must still give positive areas summing to 4*pi*a^2.
clustered = bubble_options(struct('sourceRequest',12,'n_level',5, ...
    'n_az_source',6,'n_az_colloc',8,'collocClusterRatio',0.5));
clusteredGeometry = bubble_geometry(clustered);
collocation_areas(clusteredGeometry.normals,clustered.radius);
% Check every flag combination against a direct pseudoinverse on a tiny cone.
small = bubble_options(struct('sourceRequest',12,'phi',pi/4, ...
    'n_level',1,'n_az_source',6,'n_az_colloc',8, ...
    'refinementFamilies',["S","T_t","T_n","T_n2","R"]));
g = bubble_geometry(small);
system = bubble_system(g,small,[0;0;-1]);
for useArea = [false true]
    for useColumns = [false true]
        small.useAreaWeights = useArea;
        small.useColumnScaling = useColumns;
        s = bubble_solve(small,[0;0;-1],g);
        A = (s.rowScale.*system.matrix).*s.columnScale.';
        b = s.rowScale.*system.rhs;
        expected = s.columnScale.*(pinv(A,small.svdTolerance*norm(A,2))*b);
        assert(norm(s.alpha-expected)<1e-9*max(1,norm(expected)));
        assert(norm(s.U+small.completionScale*system.forceMap*s.alpha)<1e-12);
        assert(s.forceResidual<1e-9 && all(isfinite(s.lambda)));
        if useColumns
            norms = vecnorm(A,2,1);
            assert(max(abs(norms(norms>0)-1))<1e-12);
        end
        if useArea
            assert(abs(sum(s.collocationAreas)-4*pi)<1e-10);
            weights = sqrt(s.collocationAreas/mean(s.collocationAreas));
            residual = system.matrix*s.alpha-system.rhs;
            normal = residual(1:system.normalRows)/s.Uref;
            capTangent = residual(system.normalRows+ ...
                (1:system.capTangentRows))/s.Uref;
            mobileTangent = residual(system.normalRows+system.capTangentRows+ ...
                (1:system.mobileTractionRows))/s.tractionScale;
            expectedNorm = norm([weights.*normal; ...
                repelem(weights(g.cap),2).*capTangent; ...
                repelem(weights(g.mobile),2).*mobileTangent]);
            assert(abs(s.weightedFitRMS*sqrt(size(A,1))-expectedNorm)<1e-12);
        else
            assert(isempty(s.collocationAreas) && s.fitRMS==s.weightedFitRMS);
        end
    end
end
% Equivalent cap rows must preserve physical fields with force completion.
for phi = [0 pi/3 pi]
    for tolerance = [1e-3 1e-9]
        options = bubble_options(struct('sourceRequest',12,'phi',phi, ...
            'svdTolerance',tolerance,'n_level',double(phi==pi/3), ...
            'n_az_source',6,'n_az_colloc',8));
        geometry = bubble_geometry(options);
        split = bubble_solve(options,[0;0;-1],geometry);
        cart = bubble_solve(options,[0;0;-1],geometry,[],"cartesian");
        assert(isequal(cart.matrixSize,split.matrixSize) && cart.retainedRank==split.retainedRank);
        assert(cart.forceResidual<1e-10 && cart.forceFreeResidual<1e-10 && ...
            cart.completionForceResidual<1e-12);
        assert(norm(cart.U-split.U)/split.Uref<1e-7);
        fieldSplit = bubble_evaluate(split,geometry.points,options,geometry.normals);
        fieldCart = bubble_evaluate(cart,geometry.points,options,geometry.normals);
        assert(max(vecnorm(fieldCart.velocity-fieldSplit.velocity,2,2))/split.Uref<1e-7);
        assert(max(vecnorm(fieldCart.tangentTraction-fieldSplit.tangentTraction,2,2))/ ...
            split.tractionScale<1e-7);
    end
end
options = bubble_options(struct('sourceRequest',12,'phi',pi/3,'n_level',1, ...
    'n_az_source',6,'n_az_colloc',8));
geometry = bubble_geometry(options);
for useArea = [false true]
    for useColumns = [false true]
        options.useAreaWeights = useArea;
        options.useColumnScaling = useColumns;
        split = bubble_solve(options,[0;0;-1],geometry);
        cart = bubble_solve(options,[0;0;-1],geometry,[],"cartesian");
        assert(norm(split.U-cart.U)/split.Uref<1e-7 && cart.forceResidual<1e-10);
    end
end
options.useConditionCountScaling = true;
try
    bubble_solve(options,[0;0;-1],geometry,[],"cartesian");
    error('bubble_solve:SelfTest','Anisotropic effective cap weights accepted.');
catch ME
    assert(strcmp(ME.identifier,'bubble_solve:AnisotropicCapWeights'));
end
fprintf('bubble_solve self-test passed (split and Cartesian force completion).\n');
end
