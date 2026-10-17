function result = sep15_sphere_laplace_spheroid_sweep(delta,opt)
%SEP15_SPHERE_LAPLACE_SPHEROID_SWEEP Sweep chopped-spheroid MFS sources only.
%   result = sep15_sphere_laplace_spheroid_sweep
%   result = sep15_sphere_laplace_spheroid_sweep(0.01,opt)
%   opt = struct('nodeCounts',[100,200,400],'betas',[0.1,0.3,0.5], ...
%       'widthScales',[0.5,1],'visualise',false);
%   result = sep15_sphere_laplace_spheroid_sweep(0.01,opt);
%
% Two spheres have potentials [1;-1].  Basic sources and collocation nodes
% use repository spherical designs.  Each enhancing parent is a mapped
% spherical design, chopped at the basic proxy.  No enhancing sphere caps
% are tested.  All reported residuals and coefficient sizes use max norms.
% One parent/cap is reflected onto the second sphere.  Its dimensions depend
% on the continuous proxy radius, not on the nearest sampled proxy source.
% Only the basic solution and selected best enhanced solutions are retained;
% result.summary contains the metrics for every parameter combination.
% With visualise=true, logarithmic comparison plots update during the sweep,
% rather than appearing only after all cases have finished.
%
% Verified low-strength setup for radius 0.63 and delta 0.01, potentials [1;-1]:
%   opt = struct('nodeCounts',200,'betas',0.3,'widthScales',0.5, ...
%       'basicSourceCount',1000,'basicCollocationCount',1600, ...
%       'collocationPatchLayers',24);
%   result = sep15_sphere_laplace_spheroid_sweep(0.01,opt);
% Mirrored-grid helper check: max |c_j| about 0.85, independent max residual
% about 7.3e-5 on a finer grid with additional independent contact rings.
% No coefficient normalization is used.  More sources need not reduce their
% magnitudes: nearly dependent columns can cause large cancelling charges.

if nargin < 1 || isempty(delta)
    delta = 0.01;
end
if nargin < 2 || isempty(opt)
    opt = struct();
end

% EDIT THESE SWEEP CONTROLS.  Counts are per parent surface, BEFORE chopping.
sweepControls.radius = 0.63;
sweepControls.nodeCounts = [100,200,400,640];
sweepControls.betas = [0.1,0.2,0.3,0.5,0.7,0.9];
sweepControls.widthScales = [0.2,0.5,1,1];
sweepControls.residualTarget = 1e-3;   % independent relative max residual
sweepControls.coefficientTarget = 1;  % strict <= 1, not merely "order one"

% Beta moves the outer tip from the nearest accumulation point towards the
% body surface.  Width scale multiplies both transverse axes: < 1 narrows,
% > 1 widens, and 1 recovers the original confocal parent.  Rescaling width
% changes the foci.  The continuous proxy radius can require extra stretching
% if the requested tip would not reach it; sampled proxy sources need not be
% enclosed.  Effective beta is identical on both spheres.  Invalid parents
% are skipped.  The parent and retained sources are mirrored point by point.

% FIXED GRID CONTROLS: these are held constant throughout the sweep.
% ENHANCING BASIS: "SLP" = monopoles; "SLP+DLP" adds an independent
% parent-normal dipole per cap node.  Basic proxy nodes always remain SLP.
gridControls.enhancingLaplaceBasis = "SLP";
% Positive physical length L multiplying DLP columns; moment = L*d*n.
% Fixed L gives fitted coefficients in common charge units for max norms.
gridControls.laplaceDipoleLengthScale = 1;
% Fixed switches for BOTH the basic solve and every enhanced sweep case.
% true: normalize columns in the 2-norm AFTER row weighting; undo on coefficients.
gridControls.columnScaling = false;
% true: multiply equations AND data by sqrt(local area on the COMPLETE grid),
% accounting for merged patch density; this is area-weighted least squares.
gridControls.surfaceElementPreconditioning = false;
% Cutoffs apply to the scaled matrix, but sweep targets/plots use physical
% coefficient maxima and unweighted inf-norm residuals, never scaled values.
gridControls.basicGridDistribution = "sphericaldesign";
gridControls.basicSourceCount = 700;
gridControls.basicCollocationCount = 700*1.2;
gridControls.enhancingSourceDistribution = "sphericaldesign";
gridControls.enhancingSourceShape = "spheroid";
gridControls.sourceOffset = 0.125;
gridControls.collocationPatchLayers = 18;
gridControls.collocationExponent = 3;
gridControls.collocationPatchAngle = 0.5;
gridControls.collocationMinAzimuth = 8;
% Multiply the circumference-based azimuthal count on ALL rings.  Increase
% layers for radial refinement; increase this factor for azimuthal refinement.
% Try layers 36, oversampling 2, and keep the exponent/angle unchanged first.
gridControls.collocationAzimuthOversampling = 1;
gridControls.basicSvdCutoff = 1e-13;
gridControls.enhancedSvdCutoff = 1e-13;
gridControls.testResolution = [96,60];
% Show live coefficient/residual curves and the final one-particle plots.
gridControls.visualise = true;
% Confirm exact sizes before the basic solve and EACH enhanced solve.
% false bypasses prompts only for an explicitly pre-approved sweep.
gridControls.confirmBeforeSolve = true;

% opt can override the editable sweep and fixed-grid defaults above.
names = fieldnames(sweepControls);
for k = 1:numel(names)
    name = names{k};
    if isfield(opt,name)
        sweepControls.(name) = opt.(name);
        opt = rmfield(opt,name);
    end
end
names = fieldnames(gridControls);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opt,name)
        opt.(name) = gridControls.(name);
    end
end

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'geometry','spheres'));
addpath(fullfile(projectDirectory,'laplace'));
opt = lap3d_complete_grid_options(opt);
if opt.enhancingSourceShape ~= "spheroid" || ...
        opt.enhancingSourceDistribution ~= "sphericaldesign" || ...
        opt.basicGridDistribution ~= "sphericaldesign"
    error('sep15_sphere_laplace_spheroid_sweep:WrongGrid', ...
        'This sweep uses spherical designs and chopped spheroids only.');
end
if ~isnumeric(delta) || ~isscalar(delta) || ~isreal(delta) || ...
        ~isfinite(delta) || delta <= 0 || ...
        ~isnumeric(sweepControls.radius) || ~isscalar(sweepControls.radius) || ...
        ~isreal(sweepControls.radius) || ~isfinite(sweepControls.radius) || ...
        sweepControls.radius <= opt.sourceOffset
    error('sep15_sphere_laplace_spheroid_sweep:BadGeometry', ...
        'delta must be positive; radius must be greater than sourceOffset.');
end
names = {'nodeCounts','betas','widthScales'};
for k = 1:numel(names)
    name = names{k};
    values = sweepControls.(name);
    if ~isnumeric(values) || ~isvector(values) || isempty(values) || ...
            ~isreal(values) || any(~isfinite(values)) || any(values <= 0)
        error('sep15_sphere_laplace_spheroid_sweep:BadSweepValues', ...
            '%s must be a nonempty vector of positive finite numbers.',name);
    end
    sweepControls.(name) = unique(values(:).');
end
if any(sweepControls.nodeCounts < 3) || ...
        any(sweepControls.nodeCounts ~= round(sweepControls.nodeCounts)) || ...
        any(sweepControls.betas >= 1)
    error('sep15_sphere_laplace_spheroid_sweep:BadSweepValues', ...
        'Node counts must be integers >= 3, and beta must lie between 0 and 1.');
end
names = {'residualTarget','coefficientTarget'};
for k = 1:numel(names)
    value = sweepControls.(names{k});
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value <= 0
        error('sep15_sphere_laplace_spheroid_sweep:BadTarget', ...
            '%s must be a positive finite scalar.',names{k});
    end
end

setupStart = tic;
pair = lap3d_sphere_pair(sweepControls.radius,delta);
basicGrid = lap3d_build_two_body_grids(pair,opt);

% Strengthen the independent boundary check near contact.  These test rings
% differ from collocation rings in layer count, grading, and outer angle.
% The contact node is omitted here because it is already a collocation node.
testPatchOptions = opt;
testPatchOptions.collocationPatchLayers = 2*opt.collocationPatchLayers+1;
testPatchOptions.collocationExponent = opt.collocationExponent+0.5;
testPatchOptions.collocationPatchAngle = 0.97*opt.collocationPatchAngle;
testPatchOptions.collocationMinAzimuth = opt.collocationMinAzimuth+2;
testPatchMask = false(size(basicGrid.testPoints,1),1);
bodies = {pair.E1,pair.E2};
contacts = {pair.p1,pair.p2};
for body = 1:2
    testPatch = collocation_surface_patch( ...
        bodies{body},contacts{body},testPatchOptions);
    testPatch = testPatch(2:end,:);
    basicGrid.testPoints = [basicGrid.testPoints;testPatch];
    basicGrid.testBody = [basicGrid.testBody;body*ones(size(testPatch,1),1)];
    testPatchMask = [testPatchMask;true(size(testPatch,1),1)]; %#ok<AGROW>
end
setupTime = toc(setupStart);
bodyPotentials = [1;-1];
if ~mfs_confirm_system_sizes([size(basicGrid.collocation,1), ...
        size(basicGrid.sources,1)],"basic Laplace sweep",opt.confirmBeforeSolve)
    result = [];
    return
end
basicSolution = solve_laplace_capacitance_mfs( ...
    basicGrid,bodyPotentials,opt.basicSvdCutoff);
fprintf('\nChopped-spheroid sweep: radius %.4g, delta %.4g, potentials [1;-1]\n', ...
    sweepControls.radius,delta);
fprintf('Enhancing basis %s; DLP column length scale %g\n', ...
    opt.enhancingLaplaceBasis,opt.laplaceDipoleLengthScale);
fprintf('Column scaling %s; sqrt(surface element) row weighting %s\n', ...
    string(opt.columnScaling),string(opt.surfaceElementPreconditioning));
if opt.columnScaling || opt.surfaceElementPreconditioning
    fprintf('SVD cutoffs/ranks refer to the scaled matrix; plotted maxima stay physical.\n');
end
fprintf('Basic sources/collocation per body %d/%d; fixed cutoff %.1e\n', ...
    size(basicGrid.sourceBlocks{1},1),size(basicGrid.collocationBlocks{1},1), ...
    opt.enhancedSvdCutoff);
fprintf('Basic max residual: collocation %.3e, independent %.3e; max |c_j| %.3e\n', ...
    basicSolution.collocationMaxResidual,basicSolution.testMaxResidual, ...
    basicSolution.coefficientNorm);

[nodeCounts,betas,widthScales] = ndgrid( ...
    sweepControls.nodeCounts,sweepControls.betas,sweepControls.widthScales);
rowTemplate = struct('requestedNodes',NaN,'parentNodes',NaN, ...
    'beta',NaN,'widthScale',NaN,'effectiveBeta1',NaN,'effectiveBeta2',NaN, ...
    'capNodes1',NaN,'capNodes2',NaN,'majorRadius',NaN,'minorRadius',NaN, ...
    'enhancingCollocation1',NaN,'enhancingCollocation2',NaN, ...
    'totalEnhancingCollocation',NaN, ...
    'totalSources',NaN,'totalUnknowns',NaN,'dipoleUnknowns',NaN, ...
    'totalCollocation',NaN,'retainedRank',NaN, ...
    'coefficientInf',NaN,'capCoefficientInf',NaN,'proxyCoefficientInf',NaN, ...
    'collocationMaxResidual',NaN,'testMaxResidual',NaN, ...
    'nearContactMaxResidual',NaN,'charge1',NaN,'charge2',NaN, ...
    'relativeChargeChangeInf',NaN,'gridTime',NaN,'assemblyTime',NaN, ...
    'solveTime',NaN,'testTime',NaN,'valid',false,'message',"");
rows = repmat(rowTemplate,numel(nodeCounts),1);
bestAccurate = struct();
bestAccuracy = struct();
lowestCoefficient = struct();
smallestAccurateCoefficient = inf;
smallestCoefficient = inf;
smallestResidual = inf;
geometryErrors = { ...
    'lap3d_build_two_body_grids:SourceParentOutsideBody', ...
    'lap3d_build_two_body_grids:AccumulationNotEnclosed', ...
    'lap3d_build_two_body_grids:ProxySegmentNotEnclosed', ...
    'lap3d_build_two_body_grids:EmptySourceCap', ...
    'lap3d_build_two_body_grids:OblateSourceParent', ...
    'lap3d_build_two_body_grids:TooFewCollocationNodes'};
sweepFigure = [];
if opt.visualise
    sweepFigure = lap3d_plot_spheroid_sweep( ...
        rows,sweepControls,delta,basicSolution);
end
fprintf(['  requested  beta  width  cap1/cap2  patch total    max |c_j|', ...
    '    colloc max      test max\n']);
for k = 1:numel(nodeCounts)
    caseOptions = opt;
    caseOptions.enhancingSourceCount = nodeCounts(k);
    caseOptions.sourceBeta = betas(k);
    caseOptions.sourceWidthScale = widthScales(k);
    rows(k).requestedNodes = nodeCounts(k);
    rows(k).beta = betas(k);
    rows(k).widthScale = widthScales(k);
    gridStart = tic;
    try
        [~,caseGrid] = lap3d_build_two_body_grids(pair,caseOptions,basicGrid);
    catch exception
        if ~any(strcmp(exception.identifier,geometryErrors))
            rethrow(exception)
        end
        rows(k).message = string(exception.message);
        fprintf('  %9d %5.2f %6.2f  invalid: %s\n', ...
            nodeCounts(k),betas(k),widthScales(k),exception.identifier);
        if opt.visualise && isgraphics(sweepFigure)
            lap3d_plot_spheroid_sweep( ...
                rows,sweepControls,delta,basicSolution,sweepFigure);
        end
        continue
    end
    rows(k).gridTime = toc(gridStart);
    if ~mfs_confirm_system_sizes([size(caseGrid.collocation,1), ...
            size(caseGrid.sources,1)+numel(caseGrid.dipoleSourceIndices)], ...
            sprintf('enhanced Laplace sweep case %d',k),opt.confirmBeforeSolve)
        rows(k).message = "Cancelled: remaining sweep cases not solved.";
        break
    end
    solution = solve_laplace_capacitance_mfs( ...
        caseGrid,bodyPotentials,opt.enhancedSvdCutoff);
    metadata1 = caseGrid.sourceCapMetadata{1};
    metadata2 = caseGrid.sourceCapMetadata{2};
    capMask = caseGrid.enhancingSourceMask;
    rows(k).parentNodes = metadata1.parentNodeCount;
    rows(k).effectiveBeta1 = metadata1.effectiveBeta;
    rows(k).effectiveBeta2 = metadata2.effectiveBeta;
    rows(k).capNodes1 = size(caseGrid.sourceCapBlocks{1},1);
    rows(k).capNodes2 = size(caseGrid.sourceCapBlocks{2},1);
    rows(k).enhancingCollocation1 = caseGrid.enhancingCollocationCounts(1);
    rows(k).enhancingCollocation2 = caseGrid.enhancingCollocationCounts(2);
    rows(k).totalEnhancingCollocation = caseGrid.totalEnhancingCollocation;
    rows(k).majorRadius = metadata1.parentAxes(3);
    rows(k).minorRadius = metadata1.parentAxes(1);
    rows(k).totalSources = size(caseGrid.sources,1);
    rows(k).totalUnknowns = solution.matrixSize(2);
    rows(k).dipoleUnknowns = numel(caseGrid.dipoleSourceIndices);
    rows(k).totalCollocation = solution.matrixSize(1);
    rows(k).retainedRank = solution.retainedRank;
    rows(k).coefficientInf = solution.coefficientNorm;
    rows(k).capCoefficientInf = norm(solution.sourceCoefficientMax(capMask),inf);
    rows(k).proxyCoefficientInf = norm(solution.sourceCoefficientMax(~capMask),inf);
    rows(k).collocationMaxResidual = solution.collocationMaxResidual;
    rows(k).testMaxResidual = solution.testMaxResidual;
    rows(k).nearContactMaxResidual = norm(solution.testError(testPatchMask),inf);
    rows(k).charge1 = solution.charges(1);
    rows(k).charge2 = solution.charges(2);
    rows(k).relativeChargeChangeInf = norm( ...
        solution.charges-basicSolution.charges,inf)/ ...
        max(norm(solution.charges,inf),eps);
    rows(k).assemblyTime = solution.assemblyTime;
    rows(k).solveTime = solution.solveTime;
    rows(k).testTime = solution.testEvaluationTime;
    rows(k).valid = true;
    fprintf('  %9d %5.2f %6.2f %4d/%-5d %11d %12.3e %13.3e %13.3e\n', ...
        nodeCounts(k),betas(k),widthScales(k), ...
        rows(k).capNodes1,rows(k).capNodes2, ...
        rows(k).totalEnhancingCollocation, ...
        rows(k).coefficientInf,rows(k).collocationMaxResidual, ...
        rows(k).testMaxResidual);
    if opt.visualise && isgraphics(sweepFigure)
        lap3d_plot_spheroid_sweep( ...
            rows,sweepControls,delta,basicSolution,sweepFigure);
    end

    selected.options = caseOptions;
    selected.metrics = rows(k);
    selected.grid = caseGrid;
    selected.solution = solution;
    if solution.coefficientNorm < smallestCoefficient
        smallestCoefficient = solution.coefficientNorm;
        lowestCoefficient = selected;
    end
    if solution.testMaxResidual < smallestResidual
        smallestResidual = solution.testMaxResidual;
        bestAccuracy = selected;
    end
    if solution.testMaxResidual <= sweepControls.residualTarget && ...
            solution.coefficientNorm < smallestAccurateCoefficient
        smallestAccurateCoefficient = solution.coefficientNorm;
        bestAccurate = selected;
    end
end

summary = struct2table(rows);
summary.meetsResidualTarget = summary.valid & ...
    summary.testMaxResidual <= sweepControls.residualTarget;
summary.hasOrderOneCoefficients = summary.meetsResidualTarget & ...
    summary.coefficientInf <= sweepControls.coefficientTarget;
fprintf('\nAccurate cases ranked by max coefficient (residual target %.1e):\n', ...
    sweepControls.residualTarget);
accurateRows = sortrows(summary(summary.meetsResidualTarget,:), ...
    'coefficientInf');
if isempty(accurateRows)
    fprintf('  No tested setup meets the residual target.\n');
else
    disp(accurateRows(1:min(5,height(accurateRows)), ...
        {'requestedNodes','parentNodes','capNodes1','capNodes2','beta', ...
         'widthScale','effectiveBeta1','effectiveBeta2', ...
         'totalEnhancingCollocation','coefficientInf','testMaxResidual'}));
end
fprintf('Setups meeting BOTH targets (residual <= %.1e, max coefficient <= %.3g): %d\n', ...
    sweepControls.residualTarget,sweepControls.coefficientTarget, ...
    nnz(summary.hasOrderOneCoefficients));

result.delta = delta;
result.geometry = pair;
result.bodyPotentials = bodyPotentials;
result.options = opt;
result.sweepControls = sweepControls;
result.independentPatchOptions = testPatchOptions;
result.setupTime = setupTime;
result.summary = summary;
result.basic.grid = basicGrid;
result.basic.solution = basicSolution;
result.bestAccurate = bestAccurate;
result.bestAccuracy = bestAccuracy;
result.lowestCoefficient = lowestCoefficient;
result.orderOneFound = any(summary.hasOrderOneCoefficients);
% All source sweeps keep the same basic collocation grid and patch controls.
result.enhancingCollocationCounts = [NaN;NaN];
result.totalEnhancingCollocation = NaN;
firstValid = find(summary.valid,1);
if ~isempty(firstValid)
    result.enhancingCollocationCounts = [summary.enhancingCollocation1(firstValid); ...
        summary.enhancingCollocation2(firstValid)];
    result.totalEnhancingCollocation = summary.totalEnhancingCollocation(firstValid);
    fprintf('Enhancing collocation counts: %d + %d = %d total, excluding basic nodes.\n', ...
        result.enhancingCollocationCounts,result.totalEnhancingCollocation);
end

if opt.visualise
    % Reuse the established one-particle grid/residual/source-magnitude plots.
    chosen = bestAccurate;
    if isempty(fieldnames(chosen))
        chosen = bestAccuracy;
    end
    if ~isempty(fieldnames(chosen))
        comparisonResult.caseName = "sphere spheroid sweep";
        comparisonResult.geometry = pair;
        comparisonResult.delta = delta;
        comparisonResult.basic = result.basic;
        comparisonResult.enhanced.grid = chosen.grid;
        comparisonResult.enhanced.solution = chosen.solution;
        lap3d_plot_two_body_experiment(comparisonResult);
    end
end
end
