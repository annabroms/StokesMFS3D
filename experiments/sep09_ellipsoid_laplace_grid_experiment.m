function result = sep09_ellipsoid_laplace_grid_experiment(caseName,delta,opt)
%SEP09_ELLIPSOID_LAPLACE_GRID_EXPERIMENT Compare basic and enhanced MFS grids.
%
%   result = sep09_ellipsoid_laplace_grid_experiment
%   result = sep09_ellipsoid_laplace_grid_experiment("tip",0.05)
%   result = sep09_ellipsoid_laplace_grid_experiment("side",0.05,opt)
%   result = sep09_ellipsoid_laplace_grid_experiment("oblique",0.05,opt)
%   opt.semiAxes = [1,1,1];  % optional equal-axis sphere check
%
% Exactly one two-ellipsoid geometry is treated per call.  Both solves use
% the exterior Laplace Dirichlet data [1;-1].  Shared grid construction,
% solution, validation, and plotting live in the project laplace folder.
% Without arguments, use the editable geometry/grid controls below.  This
% geometry is independent of the reflection script's editable defaults.
% opt overrides the geometry/grid controls below.

% GEOMETRY CHOICES: "tip" = major axes along the separation direction;
% "side" = parallel major axes perpendicular to separation; "oblique" =
% differently rotated bodies.  Exactly one pair per call; "all" is rejected.
% Edit the no-argument defaults below, or pass caseName and delta directly.
if nargin < 1 || isempty(caseName)
    caseName = "side";
end
% delta is the positive closest surface gap, in the same units as semiAxes.
% Smaller gaps usually need more contact nodes and careful SVD truncation.
if nargin < 2 || isempty(delta)
    delta = 0.01;
end
% opt is a scalar struct.  opt.FIELD overrides the corresponding control
% below; opt.semiAxes overrides the axes.  Unknown option names are rejected.
if nargin < 3 || isempty(opt)
    opt = struct();
end

% MAIN GEOMETRY CONTROL.  Set equal axes here for a sphere check, or pass
% opt.semiAxes.  The first axis is the major axis; require a >= b = c > 0.
% [1;1;1] gives unit spheres; [1;0.7;0.7] gives prolate ellipsoids.
semiAxes = [1;0.9;0.7];

% MAIN GRID CONTROLS.  Counts are per body.  For spherical-design grids the
% repository selects the largest design not exceeding each requested count.
% enhancingSourceCount is ignored in tensor/quasiuniform mode: use [Nu,Nv].
% ENHANCING LAPLACE BASIS: "SLP" = monopoles only; "SLP+DLP" adds one
% independent parent-normal dipole at EACH enhancing node.  Basic proxy
% nodes always remain monopoles.  Combined mode adds one unknown per cap node.
gridControls.enhancingLaplaceBasis = "SLP+DLP";
% Positive physical length L multiplying DLP columns.  Fitted SLP/DLP
% coefficients then share charge units; physical dipole moment is L*d*n.
% Keep L fixed when comparing coefficient max norms; no quadrature weights.
gridControls.laplaceDipoleLengthScale = 1;

% Independent preconditioning switches, applied to BOTH basic/enhanced solves.
% true: normalize each matrix column in the 2-norm AFTER any row weighting.
% The scaling is undone on coefficients before magnitude plots and reporting.
gridControls.columnScaling = false;
% true: multiply each equation and its data by sqrt(local surface area).
% Areas account for the density of the COMPLETE merged collocation grid.
% This fits an area-weighted residual, not the maximum pointwise residual.
gridControls.surfaceElementPreconditioning = true;
% SVD cutoffs below act on the resulting scaled matrix. Reported residuals
% stay unweighted and in the inf norm; false/false preserves the old solve.
% BASIC GRID CHOICES: "sphericaldesign" maps repository sphere designs;
% "ellipsoidquadrature" uses getEllipsoidGrids and coarseSourceResolution.
gridControls.basicGridDistribution = "ellipsoidquadrature"; %"sphericaldesign";
% ENHANCING GRID CHOICES: "tensor" = a fixed azimuth count per parent ring;
% "quasiuniform" = tapered ring counts near tips; "sphericaldesign" = an
% affine-mapped whole-parent design, with no independent ring/layer counts.
gridControls.enhancingSourceDistribution =  "quasiuniform"; %"tensor";
% Requested basic sources PER BODY (integer >= 3); sphericaldesign only.
% Actual available design counts may be smaller.  Does not set cap density.
gridControls.basicSourceCount = 700;
% Requested base boundary nodes PER BODY (integer >= 3); sphericaldesign
% only.  These enter the solve; keep collocation more numerous than sources.
gridControls.basicCollocationCount = 1400;
gridControls.basicCollocationCount = 1000;
% Requested WHOLE-parent design count PER BODY, before chopping (integer
% >= 3).  Ignored for tensor/quasiuniform; actual cap counts are printed.
gridControls.enhancingSourceCount = 1500;

% BASIC PROXY SOURCE RESOLUTION ONLY: [azimuth scale, Gauss-Legendre layers].
% Two integers >= 4; ignored when basicGridDistribution = "sphericaldesign".
% basicSourceCount/basicCollocationCount are ignored in quadrature mode.
gridControls.coarseSourceResolution = [28,18];
gridControls.coarseSourceResolution = [35,28];

% BASIC COLLOCATION RESOLUTION ONLY: independently editable [Nu,Nv].
% Active when basicGridDistribution = "ellipsoidquadrature".
% Nu controls nodes AROUND each ring (fewer on small rings near the poles);
% Nv controls Gauss-Legendre layers from one pole to the other.  Both are
% integers >= 4, rounded to even.  The total is NOT simply Nu*Nv.
% Increase Nu to refine around rings, Nv to add layers, or both to refine
% the whole boundary.  Changes neither proxy/cap sources nor test nodes.
% [70,56] preserves the previously automatic grid for sources [35,28].
% This explicit choice no longer follows changes to coarseSourceResolution.
% Set [] only if you want getEllipsoidGrids to choose its default again.
gridControls.basicCollocationResolution = [70,56];
gridControls.basicCollocationResolution = [round(1.5*35),round(1.5*28)];


% ENHANCING SURFACE CHOICES: "spheroid" = stretched, chopped prolate parent;
% "sphere" = a chopped sphere centred at the nearest accumulation point.
% The basic proxy is retained in either case.  Default is the spheroid.
gridControls.enhancingSourceShape = "spheroid";
% For shape "sphere", enter a positive radius enclosing all accumulation
% points while staying inside the body.  For "spheroid", this MUST be [].
gridControls.sourceSphereRadius = [];

% CHOPPED-SPHEROID SOURCE RESOLUTION, NOT the basic proxy resolution.
% Tensor mode makes the cross-sectional count below directly adjustable.
% [Nu,Nv] then means Nu points AROUND each cross-sectional ring and Nv
% Gauss-Legendre layers ALONG the complete parent spheroid's long axis.
% Layers cluster towards the tips; counts are rounded to even integers.
% The parent has Nu*Nv nodes before chopping.  "quasiuniform" instead
% decreases the azimuthal count on the small rings near the tips.
% This resolution is ignored for "sphericaldesign"; use the count above.
% Both entries must be integers >= 4.  Increasing the first resolves each
% cross section; increasing the second adds axial layers before chopping.
gridControls.fineSourceResolution = [52,30]; %[32,20];

% Quasiuniform only: add azimuths on all rings or raise the small-ring floor.
% Tensor mode ignores these two controls and uses exactly Nu per parent ring.
% ExtraAzimuthNodes: integer >= 0, added to the circumference-scaled count
% before rounding/floor.  Larger = denser ring sampling.
gridControls.enhancingSourceExtraAzimuthNodes = 0;
% MinAzimuthNodes: integer >= 2; raises the soft floor on small near-tip rings.
gridControls.enhancingSourceMinAzimuthNodes = 8;

% SOLVE GRID: basicCollocationResolution sets the basic quadrature boundary;
% basicCollocationCount sets it only when using spherical designs instead.
% These controls add the near-contact collocation patch used in the solve.
% Layers: positive integer; more rings refine radially and also add azimuths.
gridControls.collocationPatchLayers = 25; % number of graded rings
% Exponent: positive; small values are nearly uniform, larger values cluster
% towards contact.  rho_j = angle*expm1(exponent*j/L)/expm1(exponent).
gridControls.collocationExponent = 5; % larger = stronger contact clustering
% Angle: 0 < angle < pi, on the parameter sphere.  Larger = wider patch,
% not automatically finer spacing.  0.5 rad is approximately 28.6 degrees.
gridControls.collocationPatchAngle = 0.5; %0.5; % patch half-angle in radians
% Minimum: integer >= 4.  Actual counts adapt to circumference/radial spacing
% and become even.  Raising this especially refines the small inner rings.
gridControls.collocationMinAzimuth = 8; % minimum nodes around each ring
% Oversampling: positive; 1 = standard, 2 approximately doubles adaptive
% azimuthal counts.  Increase this for cross-sectional boundary resolution.
gridControls.collocationAzimuthOversampling = 1; % multiply adaptive ring counts

% ACCURACY-CHECK GRID ONLY: these nodes never enter either solve.
% [Nu,Nv] sets the full-body evaluation grid via repository quadrature;
% it is quasiuniform in azimuth, so Nu is an approximate maximum per ring.
% The original [64,40] is increased to [128,80] for a finer global check.
% Entries are integers >= 4, rounded to even counts.  Larger = finer testing
% only, not a more accurate fitted solution; it may reveal larger error peaks.
gridControls.testResolution = [128,80];
% Add independently graded near-contact evaluation rings: more rings and
% azimuths than collocation, with different radii to avoid testing only fit nodes.
% Increase both controls until the measured max residual stops changing.
% Test rings use exponent+0.5, angle*0.97 and minimum azimuth+2, rather than
% simply rechecking the fitted rings.  Layers is an integer >= 0.
gridControls.testPatchLayers = 48; % 0 disables extra near-contact test rings
% Positive factor for independent test-ring azimuths; ignored when layers=0.
gridControls.testPatchAzimuthOversampling = 1;

% SOURCE DUPLICATES: false skips duplicate checks/removal on the enhancing
% sources, including overlap with basic proxy nodes.  true removes near-
% duplicates from the cap; neither choice changes basic sources.  Collocation nodes
% are still deduplicated.  Can also be overridden with opt.removeDuplicates.
gridControls.removeDuplicates = false;

% A deeper proxy smooths the basic source representation.  Offset 0.25
% still encloses the [1,0.7,0.7] focal segment with about 0.03586 clearance.
% It also leaves a nonempty enhancing cap in the equal-axis sphere check.
% Positive inward NORMAL displacement in physical length units (not beta).
% Too large an offset can lose focal enclosure or source interiority.
gridControls.sourceOffset = 0.25;
% Beta: 0 < beta < 1; fraction extending the parent tip from the nearest
% accumulation point towards its radial body boundary.  Larger approaches
% the boundary.  Proxy enclosure may stretch further (effectiveBeta).
% Used for shape "spheroid" only; "sphere" uses sourceSphereRadius instead.
gridControls.sourceBeta = 0.7; % tip fraction from accumulation to boundary

% CAP EXTENT TOWARDS THE BASIC PROXY (physical length units, not beta).
% Retain parent nodes whose distance to the body surface is less than
% sourceOffset + sourceChopOffset.  0 cuts at the continuous basic proxy;
% +0.02 extends 0.02 deeper into it; -0.02 leaves a 0.02 clearance gap.
% This changes only the cut, not the parent shape or contact-facing tip.
% A discrete grid need not have a node exactly at the cut: increase the
% design count, or Nv in tensor mode, for a better-resolved cap edge.
% Require sourceOffset + sourceChopOffset > 0.  Excessive proxy overlap
% can increase redundancy and source coefficients, rather than accuracy.
gridControls.sourceChopOffset = 0.05;

% SPHEROID WIDTH: multiply both transverse axes by this factor.
% 1 = unchanged, 0.8 = narrower, 1.2 = wider.  The ellipsoid proxy-enclosure
% rule may also stretch the axial direction, changing effective beta.
% The parent must remain prolate, inside the body, and enclose accumulation.
% Positive factor, used for "spheroid" only; "sphere" ignores this setting.
gridControls.sourceWidthScale = 0.8;

% SVD TRUNCATION: retain singular values sigma >= cutoff*sigma_max.
% Larger cutoffs discard more modes, usually reducing coefficient growth
% at the cost of boundary accuracy.  The two solves are independently tunable.
% Both cutoffs must lie in (0,1); typical trials are 1e-13, 1e-11 and 1e-9.
% Basic cutoff affects only the solve on the original proxy/boundary grids.
gridControls.basicSvdCutoff = 1e-12;
% Enhanced cutoff affects only the3 solve with added cap/patch nodes.
gridControls.enhancedSvdCutoff = 1e-9; % try 1e-11 for stronger truncation

% PLOTTING: true shows the comparison, grids and independently scaled
% source/residual scatter plots; false skips figures.  No automatic saving.
gridControls.visualise = true;

% CONSENT BEFORE SOLVING: true shows exact basic/enhanced matrix dimensions
% and matrix-only memory, then waits for y/yes.  Anything else cancels with [].
% Geometry/grids are built to determine sizes; kernels/SVD/evaluations wait.
% false skips the prompt ONLY for an already-approved automated run.
gridControls.confirmBeforeSolve = true;

caseName = lower(string(caseName));

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'geometry','spheres'));
addpath(fullfile(projectDirectory,'laplace'));

[pair,opt,pairBuildTime] = prepare_ellipsoid_grid_experiment( ...
    caseName,delta,semiAxes,gridControls,opt);

result = lap3d_run_two_body_grid_experiment(pair,caseName,opt);
if isempty(result)
    return % user declined; do not plot or access missing solution fields
end
result.options.semiAxes = pair.E1.axes;
result.pairBuildTime = pairBuildTime;
result.geometryTime = pairBuildTime+result.gridBuildTime;
if opt.visualise
    lap3d_plot_two_body_experiment(result);
end
end
