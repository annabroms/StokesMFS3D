function result = sep10_sphere_laplace_grid_experiment(radius,delta,opt)
%SEP10_SPHERE_LAPLACE_GRID_EXPERIMENT Compare basic and enhanced sphere grids.
%
%   result = sep10_sphere_laplace_grid_experiment
%   result = sep10_sphere_laplace_grid_experiment(0.7,0.05)
%   result = sep10_sphere_laplace_grid_experiment(0.7,0.05,opt)
%
% Two equal spheres are placed on the x axis with closest gap delta.  The
% basic proxy and base collocation nodes are spherical designs.  Additional
% sources are spherical-design nodes mapped to the parent surface and then
% chopped at the basic proxy.  The Laplace solver and plotting are shared
% with the ellipsoid experiment.  The radius and grid controls are editable
% below; opt can override the grid controls for an individual invocation.
% opt.comparisonSphereRadiusFactors optionally sets the two comparison
% radii as fractions of the accumulation-point distance to the body.

if nargin < 1 || isempty(radius)
    radius = 0.63;
end
if nargin < 2 || isempty(delta)
    delta = 0.01;
end
if nargin < 3 || isempty(opt)
    opt = struct();
end

% Grid controls.  These are the main lines to edit.
% Counts are PER BODY.  get_sphdesign selects the largest available design
% not exceeding the requested count.  Enhancing counts refer to the whole
% parent surface BEFORE chopping; the actual retained counts are printed.
% ENHANCING BASIS: "SLP" = monopoles; "SLP+DLP" = independently fitted
% monopole and parent-normal dipole per enhancing node.  Basic nodes stay SLP.
% Combined mode adds one unknown per cap node, also on comparison spheres.
gridControls.enhancingLaplaceBasis = "SLP";
% Positive physical length L multiplying DLP columns; moment = L*d*n.
% Keep fixed to compare max fitted coefficients in common charge units.
gridControls.laplaceDipoleLengthScale = 1;
% Apply these independently to basic, spheroid and BOTH enhancing sphere solves.
% true: unit 2-norm columns AFTER row weighting; recover physical coefficients.
gridControls.columnScaling = false;
% true: weight equations AND data by sqrt(local area per collocation node),
% including the density of the merged contact patch (no overlap double count).
% This is area-weighted least squares, not minimization of the max residual.
gridControls.surfaceElementPreconditioning = false;
% SVD cutoffs apply to the scaled matrix. Magnitude plots and reported
% inf-norm residuals use physical coefficients and unweighted potentials.
gridControls.basicGridDistribution = "sphericaldesign";
gridControls.basicSourceCount = 700;
gridControls.basicCollocationCount = 700*1.2;
gridControls.enhancingSourceDistribution = "sphericaldesign";
gridControls.enhancingSourceCount = 640;

% Separate cross-sectional/axial counts require the "tensor" alternative.
% [Nu,Nv] = points around each ring / Gauss-Legendre layers along the FULL
% parent long axis.  Counts become even; layers cluster towards both tips.
% "quasiuniform" reduces azimuthal counts on the smaller near-tip rings.
% For "sphericaldesign", this is ignored: use enhancingSourceCount instead.
gridControls.fineSourceResolution = [32,20];

% Chopping depth = sourceOffset + sourceChopOffset (physical length units).
% 0 cuts at the continuous basic proxy; positive values extend the cap
% farther INTO the proxy, negative values stop short and leave a gap.
% For example, +0.02 retains sources 0.02 deeper than the original cut.
% Only the cut changes, not beta, width, or the parent contact-facing tip.
% Increase the design count (or Nv in tensor mode) to resolve the cut edge.
gridControls.sourceChopOffset = 0;

% INDEPENDENT ACCURACY GRID ONLY: these nodes do not enter either solve.
gridControls.testResolution = [128,80]; % full-body [azimuth, axial] resolution
gridControls.testPatchLayers = 48; % independently graded contact rings; 0 disables
gridControls.testPatchAzimuthOversampling = 2; % adaptive azimuthal ring counts
% true waits for y/yes after showing system sizes; false is for an already-
% approved automated run.  Each additional comparison surface is confirmed.
gridControls.confirmBeforeSolve = true;

% Spheroid beta: put the contact-facing tip this fraction of the distance
% from the nearest accumulation point towards the body's radial surface.
% 0 < beta < 1; 0.3 means 30% of that distance.  Increasing beta stretches
% the spheroid and moves its tip closer to the boundary.  It may be enlarged
% further if needed to reach the continuous proxy radius R-sourceOffset.
% Sampled proxy sources need not be enclosed.  One parent/cap is reflected
% onto the second sphere, so both effective betas and source counts match.
% A parent surface leaving the body is rejected.  The two comparison-sphere
% radii are independent of beta.
gridControls.sourceBeta = 0.3;

% Multiply both transverse spheroid axes: < 1 is narrower, > 1 is wider.
% The beta-controlled tip is unchanged by width scaling.  The proxy-radius
% bound can require stretching, independently of width or sampled sources.
% Width scaling changes the foci; the accumulation points must
% remain enclosed, and the parent must stay prolate and inside the body.
gridControls.sourceWidthScale = 1;

% Collocation patch: extra boundary nodes around each closest contact point.
% Number of concentric rings, plus one node at the contact point.  More
% layers make radial spacing finer and normally add more collocation nodes.
gridControls.collocationPatchLayers = 18;

% Positive exponential grading strength.  Larger values cluster the inner
% rings more tightly towards contact; small positive values are nearly uniform.
% Ring angles: rho_j = patchAngle*(exp(alpha*j/L)-1)/(exp(alpha)-1),
% where alpha is this exponent, L is the layer count, and j = 1,...,L.
gridControls.collocationExponent = 3;

% Patch half-angle in RADIANS on the parameter sphere, not a physical length.
% 0.5 rad is about 28.6 degrees.  On a sphere of radius R, the outer ring is
% R*patchAngle in surface arc length from contact.  Larger angles cover more
% of the body; the final ring is always at this angle.  Require 0 < angle < pi.
gridControls.collocationPatchAngle = 0.5;

% Minimum azimuthal nodes on every ring (integer >= 4).  Actual ring counts
% adapt to circumference/radial spacing and are rounded up to even integers.
% Increase this minimum to give the small rings more azimuthal resolution.
gridControls.collocationMinAzimuth = 8;

% Multiply the circumference-based azimuthal counts on ALL rings.  A factor
% of 2 roughly doubles those counts; the minimum above is only a floor.
% For a denser patch, try layers 36 and oversampling 2 before changing the
% exponent.  Always compare the independent max residual, not only collocation.
gridControls.collocationAzimuthOversampling = 1;

% Values supplied in opt take precedence over the editable defaults above.
controlNames = fieldnames(gridControls);
for k = 1:numel(controlNames)
    name = controlNames{k};
    if ~isfield(opt,name)
        opt.(name) = gridControls.(name);
    end
end

if ~isscalar(radius) || ~isreal(radius) || ~isfinite(radius) || radius <= 0
    error('sep10_sphere_laplace_grid_experiment:BadRadius', ...
        'radius must be a positive finite scalar.');
end
if ~isscalar(delta) || ~isreal(delta) || ~isfinite(delta) || delta <= 0
    error('sep10_sphere_laplace_grid_experiment:BadGap', ...
        'delta must be a positive finite scalar.');
end

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'geometry','spheres'));
addpath(fullfile(projectDirectory,'laplace'));

comparisonRadiusFactors = [0.3,0.6];
if isfield(opt,'comparisonSphereRadiusFactors')
    comparisonRadiusFactors = opt.comparisonSphereRadiusFactors;
    opt = rmfield(opt,'comparisonSphereRadiusFactors');
end
if ~isnumeric(comparisonRadiusFactors) || ...
        numel(comparisonRadiusFactors) ~= 2 || ...
        any(~isfinite(comparisonRadiusFactors)) || ...
        any(comparisonRadiusFactors <= 0) || ...
        any(comparisonRadiusFactors >= 1) || ...
        comparisonRadiusFactors(1) == comparisonRadiusFactors(2)
    error('sep10_sphere_laplace_grid_experiment:BadSphereRadiusFactors', ...
        ['comparisonSphereRadiusFactors must contain two distinct ', ...
         'numbers strictly between zero and one.']);
end
comparisonRadiusFactors = comparisonRadiusFactors(:).';

opt = lap3d_complete_grid_options(opt);
if opt.sourceOffset >= radius
    error('sep10_sphere_laplace_grid_experiment:BadSourceOffset', ...
        'sourceOffset must be smaller than the sphere radius.');
end

pairStart = tic;
pair = lap3d_sphere_pair(radius,delta);
pairBuildTime = toc(pairStart);

result = lap3d_run_two_body_grid_experiment(pair,"sphere",opt);
if isempty(result)
    return % declined the basic/enhanced solves
end
result.sphereRadius = radius;
result.pairBuildTime = pairBuildTime;
result.geometryTime = pairBuildTime+result.gridBuildTime;
result.options.comparisonSphereRadiusFactors = comparisonRadiusFactors;

methodTemplate = struct('name',"",'sourceRadius',NaN, ...
    'radiusFactor',NaN,'grid',struct(),'solution',struct(), ...
    'comparisonToBasic',struct());
methods = repmat(methodTemplate,1,3);
methods(1).name = "chopped spheroid";
methods(1).grid = result.enhanced.grid;
methods(1).solution = result.enhanced.solution;
methods(1).comparisonToBasic = result.comparison;

accumulationDistance = norm(pair.p1-pair.Pminus(:,1));
for k = 1:2
    sphereOpt = opt;
    sphereOpt.enhancingSourceShape = "sphere";
    sphereOpt.sourceSphereRadius = ...
        comparisonRadiusFactors(k)*accumulationDistance;
    [~,sphereGrid] = lap3d_build_two_body_grids(pair,sphereOpt);
    if ~mfs_confirm_system_sizes([size(sphereGrid.collocation,1), ...
            size(sphereGrid.sources,1)+numel(sphereGrid.dipoleSourceIndices)], ...
            sprintf('comparison sphere radius %.6g', ...
            sphereOpt.sourceSphereRadius),opt.confirmBeforeSolve)
        result = [];
        return
    end
    sphereSolution = solve_laplace_capacitance_mfs( ...
        sphereGrid,result.bodyPotentials,sphereOpt.enhancedSvdCutoff);

    comparison.collocationMaxResidualRatio = ...
        sphereSolution.collocationMaxResidual/ ...
        max(result.basic.solution.collocationMaxResidual,eps);
    comparison.testMaxResidualRatio = sphereSolution.testMaxResidual/ ...
        max(result.basic.solution.testMaxResidual,eps);
    comparison.relativeChargeDifference = norm( ...
        sphereSolution.charges-result.basic.solution.charges,inf)/ ...
        max(norm(sphereSolution.charges),eps);

    methods(k+1).name = sprintf('sphere, radius %.4g', ...
        sphereOpt.sourceSphereRadius);
    methods(k+1).sourceRadius = sphereOpt.sourceSphereRadius;
    methods(k+1).radiusFactor = comparisonRadiusFactors(k);
    methods(k+1).grid = sphereGrid;
    methods(k+1).solution = sphereSolution;
    methods(k+1).comparisonToBasic = comparison;
end
result.sourceShapeComparisons = methods;

fprintf('\nEnhancing source-surface comparison (inf-norm residual ratios to basic):\n');
fprintf('  %-25s %12s %12s %12s %12s %8s\n', ...
    'source surface','colloc ratio','test ratio','max |c_j|','added/body','rank');
for k = 1:numel(methods)
    comparison = methods(k).comparisonToBasic;
    addedSourceCount = size(methods(k).grid.sourceCapBlocks{1},1);
    fprintf('  %-25s %12.3e %12.3e %12.3e %12d %8d\n', ...
        methods(k).name,comparison.collocationMaxResidualRatio, ...
        comparison.testMaxResidualRatio, ...
        methods(k).solution.coefficientNorm,addedSourceCount, ...
        methods(k).solution.retainedRank);
end

if opt.visualise
    lap3d_plot_two_body_experiment(result);
end
end
