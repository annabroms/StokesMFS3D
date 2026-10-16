function result = sep15_ellipsoid_stokes_grid_experiment(caseName,delta,opt)
%SEP15_ELLIPSOID_STOKES_GRID_EXPERIMENT Stokes counterpart of the Laplace test.
%   result = sep15_ellipsoid_stokes_grid_experiment
%   opt.enhancingFamilies = 'SRTD';
%   opt.semiAxes = [1,1,1];  % equal-axis sphere check
%   result = sep15_ellipsoid_stokes_grid_experiment("tip",0.01,opt);
% To match a completed Laplace run exactly, including edited grid settings:
%   result = sep15_ellipsoid_stokes_grid_experiment( ...
%       laplaceResult.caseName,laplaceResult.delta,laplaceResult.options);
% One configuration per call.  Default motions translate towards each other
% along the gap normal with zero rotation.  Viscosity is 1; velocity decays
% at infinity.  Geometry, sources, collocation and evaluation reuse Laplace
% helpers; Stokes solving/reporting/plotting are shared with the sphere test.
% Dense Stokes systems have three unknowns per node per selected family.
% Start with smaller counts when trying many families.  No files are saved.
% Geometry/grid defaults below match the current ellipsoid Laplace script.
% This is a settings copy, not automatic synchronization of future edits.

% GEOMETRY CHOICES: "tip" = separation along the major axes; "side" =
% separation across parallel major axes; "oblique" = differently rotated
% bodies.  One pair per call; "all" is rejected.  Arguments override defaults.
if nargin < 1 || isempty(caseName)
    caseName = "tip";
end
% Positive closest surface gap, in the same length units as semiAxes.
% Smaller gaps usually require finer contact grids and careful truncation.
if nargin < 2 || isempty(delta)
    delta = 0.01;
end
% Scalar option struct: opt.FIELD overrides the controls below.
% Unknown names are rejected.  Passing laplaceResult.options matches its grids.
if nargin < 3 || isempty(opt)
    opt = struct();
end
% Positive physical-body axes [a;b;c]; opt.semiAxes overrides this.
% [1;1;1] gives unit spheres; three distinct axes give a triaxial ellipsoid.
% The accumulation model uses all three axes; the separate focal-segment
% diagnostic remains specific to prolate spheroids [a;b;b].
semiAxes = [1;0.9;0.7];

% BASIC GRID CHOICES: "sphericaldesign" maps repository sphere designs;
% "ellipsoidquadrature" uses getEllipsoidGrids and coarseSourceResolution.
% These are independent of the enhancing source distribution below.
gridControls.basicGridDistribution = "ellipsoidquadrature";
% Requested original proxy sources PER BODY, integer >= 3; design mode only.
% Available design counts may be smaller.  Basic nodes always carry S only.
gridControls.basicSourceCount = 700;
% Requested base boundary nodes PER BODY, integer >= 3; design mode only.
% These enter the least-squares solve, unlike testResolution/testPatch nodes.
gridControls.basicCollocationCount = 1000;
% Basic PROXY SOURCE quadrature: [azimuth scale, Gauss-Legendre layers],
% integers >= 4.  Ignored in design mode; the two basic design counts are
% ignored in quadrature mode.  Boundary resolution is independent below.
gridControls.coarseSourceResolution = [35,28];
% BASIC COLLOCATION RESOLUTION ONLY: independently editable [Nu,Nv].
% To use this control, set basicGridDistribution = "ellipsoidquadrature".
% Nu controls nodes AROUND each ring (fewer near the poles); Nv controls
% Gauss-Legendre layers from pole to pole.  Integers >= 4, rounded to even.
% The total is NOT Nu*Nv.  Increase either entry, or both, for finer fitting.
% This matches the Laplace script's last assignment, before even rounding.
% For example, use [80,60] for a finer boundary without changing sources.
% [] keeps repository oversampling; ignored in spherical-design mode, where
% basicCollocationCount above is the active boundary-node control instead.
gridControls.basicCollocationResolution = [round(1.5*35),round(1.5*28)];

% SOURCE SURFACE CHOICES: "spheroid" = stretched chopped prolate parent;
% "sphere" = chopped sphere centred at the nearest accumulation point.
% Both retain the original proxy.  Sphere mode ignores beta/width below.
gridControls.enhancingSourceShape = "spheroid";
% Leave [] for "spheroid".  For "sphere", enter a positive radius enclosing
% every accumulation point while keeping the parent inside the body.
gridControls.sourceSphereRadius = [];
% SOURCE DISTRIBUTIONS: "tensor" = fixed azimuths per parent ring;
% "quasiuniform" = fewer azimuths on small rings near tips;
% "sphericaldesign" = whole-parent design, controlled by enhancingSourceCount.
gridControls.enhancingSourceDistribution = "quasiuniform"; % tapered ring counts
% [nodes AROUND each cross section, axial layers along the FULL parent].
% Chopping removes deep layers/nodes.  Never changes the basic proxy grid.
% Integers >= 4, rounded to even.  Tensor has Nu*Nv parent nodes; quasiuniform
% uses Nu as an azimuthal scale, not an exact count.  Axial layers cluster
% towards both tips.  Ignored in sphericaldesign mode.
gridControls.fineSourceResolution = [42,30];
% Integer >= 0, added to the circumference-scaled azimuth count before
% rounding/floor.  Increase to refine all quasiuniform rings; otherwise ignored.
gridControls.enhancingSourceExtraAzimuthNodes = 0; % quasiuniform only
% Integer >= 2: soft azimuthal floor on small rings; quasiuniform only.
gridControls.enhancingSourceMinAzimuthNodes = 8; % quasiuniform small-ring floor
% Integer >= 3, PER BODY for the WHOLE parent BEFORE chopping.  Available
% designs may be smaller; ignored in tensor/quasiuniform mode.
gridControls.enhancingSourceCount = 1200; % used only for "sphericaldesign"
% true removes near-duplicate enhancing sources, including overlap with the
% basic proxy; false keeps every chopped source.  Basic sources stay untouched;
% collocation duplicates are always removed, independently of this flag.
gridControls.removeDuplicates = false; % source nodes only
% Positive inward NORMAL displacement in physical length units.  Larger
% offsets put proxies deeper; excessive offsets can lose focal enclosure.
% Sources must remain strictly inside the physical bodies.
gridControls.sourceOffset = 0.25;
% 0 < beta < 1: requested fraction of the nearest-accumulation-to-radial-
% boundary gap extending the spheroid tip.  Larger moves towards the boundary.
% Sampled-proxy enclosure can stretch further, giving effectiveBeta != beta
% and different values on nonspherical bodies.  Ignored for shape "sphere".
gridControls.sourceBeta = 0.2; % fraction moving the parent tip towards boundary
% Positive transverse-axis multiplier: 1 unchanged, 0.8 narrower, 1.2 wider.
% Changes the foci; proxy enclosure may stretch the axial direction too.
% Parent must remain prolate, enclose accumulation and stay inside the body.
% Ignored for shape "sphere", whose radius is set independently above.
gridControls.sourceWidthScale = 0.8; % transverse-axis multiplier
% Signed PHYSICAL length beyond the original proxy chopping depth:
% 0 cuts at the proxy; positive extends the cap deeper INTO it; negative
% stops short and leaves a gap.  Require sourceOffset+sourceChopOffset > 0.
% Changes only the cut, not the tip/parent.  Overlap may increase redundancy.
gridControls.sourceChopOffset = 0.05; % signed depth beyond the proxy cut

% COLLOCATION PATCH: physical boundary nodes entering the enhanced solve.
% Positive integer number of graded rings, plus one contact node.  More
% layers refine radial spacing and normally increase azimuthal counts too.
gridControls.collocationPatchLayers = 25;
% Positive exponent; small = nearly uniform, larger = clustered at contact.
% rho_j = angle*expm1(exponent*j/L)/expm1(exponent), j=1,...,L.
gridControls.collocationExponent = 5;
% 0 < angle < pi RADIANS on the parameter sphere.  Larger covers a wider
% patch, not automatically finer spacing; 0.5 rad is about 28.6 degrees.
gridControls.collocationPatchAngle = 0.5;
% Integer >= 4, azimuthal floor per ring; actual counts adapt to circumference
% and radial spacing and become even.  Raising this refines small inner rings.
gridControls.collocationMinAzimuth = 8;
% Positive azimuthal factor: 1 standard; 2 approximately doubles adaptive
% ring counts.  The total velocity equations must exceed ALL family unknowns.
gridControls.collocationAzimuthOversampling = 1;
% SVD cutoffs in (0,1), relative to the system's largest singular value.
% Larger discards more modes, usually lowering coefficients but potentially
% worsening boundary accuracy; try 1e-13, 1e-11, 1e-9.  Basic solve only:
gridControls.basicSvdCutoff = 1e-12;
% Independent cutoff for the enhanced solve, including all selected families.
gridControls.enhancedSvdCutoff = 1e-12;

% INDEPENDENT ACCURACY CHECK ONLY: no effect on either fitted solution.
% [Nu,Nv], integers >= 4, rounded to even: whole-body repository quadrature.
% Nu is an approximate maximum per ring (quasiuniform), Nv the layer count.
% Increase both until the measured inf-norm residual stabilizes.
gridControls.testResolution = [128,80]; % full-body [azimuth, axial] resolution
% Integer >= 0; 0 disables extra test rings.  Uses exponent+0.5, angle*0.97,
% minimum azimuth+2 to avoid merely rechecking the fitted rings.  Contact is
% checked too.  Increasing this may reveal larger previously missed peaks.
gridControls.testPatchLayers = 48; % independent graded contact rings; 0 disables
% Positive independent test-ring azimuth factor; larger = denser checking.
% Ignored when testPatchLayers=0.  Both solves use exactly the same test nodes.
gridControls.testPatchAzimuthOversampling = 2;
% true shows body-1 source/residual plots, marking nearest accumulation black;
% false skips figures.  Each colour scale is individual; no automatic saving.
gridControls.visualise = true;
% true previews both dense matrix sizes, including THREE components per
% node and ALL enhancing families, and waits for y/yes before kernels/SVDs.
% Geometry/grids are built first for exact sizes.  Other answers return [].
% false is an explicit bypass ONLY for an already-approved automated run.
gridControls.confirmBeforeSolve = true;

% STOKES CONTROLS.  Basic nodes always use S (Stokeslets).
% At enhancing nodes choose any combination: S = Stokeslet, R = rotlet,
% T = stresslet with parent normal, D = potential dipole; 3 unknowns each.
% Any nonempty combination of DISTINCT letters: e.g. 'S', 'SD', 'SRTD';
% also accepts ["S","D"].  More families increase dense-SVD cost and can
% require more collocation equations.  Letter order is canonicalized.
stokesControls.enhancingFamilies = 'SRD'; % for example 'SD' or 'SRTD'
% Default translation is resolved along the geometry's gap normal below;
% pass opt.translations (two rows) to prescribe another rigid translation.
% This zero matrix is a placeholder, NOT stationary particles: it is replaced
% by [gapNormal';-gapNormal'] unless opt.translations is supplied.
% For explicit Cartesian velocities, uncomment/edit the line below:
% opt.translations = [1,0,0;-1,0,0]; % rows = bodies, columns = x,y,z
stokesControls.translations = zeros(2,3);
% Angular velocities: finite 2-by-3 matrix, rows = bodies, columns = x,y,z.
% zeros means no rotation; positive directions follow the right-hand rule
% about each body centre.  Edit here or supply opt.angularVelocities.
stokesControls.angularVelocities = zeros(2,3);
% Positive length L: multiply R/T columns by L and D columns by L^2, putting
% fitted coefficients in common force units.  Changing L changes coefficient
% max norms and possibly the truncated solution; physical strengths are also
% returned.  Independent of source geometry; not a trick to make strengths small.
stokesControls.singularityLengthScale = semiAxes(1); % R/T scale L; D scale L^2
% Positive integer number of evaluation targets per batch.  Smaller uses less
% temporary memory, larger may be faster.  Does NOT change grid density/accuracy.
stokesControls.evaluationBatchSize = 256; % memory control, not node resolution

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory,fullfile(projectDirectory,'common','dense'), ...
    fullfile(projectDirectory,'geometry','spheres'), ...
    fullfile(projectDirectory,'geometry','ellipsoids'),fullfile(projectDirectory,'laplace'));
assert(isstruct(opt) && isscalar(opt),'opt must be a scalar struct.');
translationsSpecified = isfield(opt,'translations');
[stokesControls,opt] = complete_stokes_two_body_controls(opt,stokesControls);
caseName = lower(string(caseName));
[pair,opt,pairBuildTime] = prepare_ellipsoid_grid_experiment( ...
    caseName,delta,semiAxes,gridControls,opt);
if ~translationsSpecified
    stokesControls.translations = [pair.n.';-pair.n.'];
end
result = run_stokes_two_body_grid_experiment(pair,caseName,opt,stokesControls);
if isempty(result)
    return % declined: no solution or figures to process
end
result.options.semiAxes = pair.E1.axes;
result.pairBuildTime = pairBuildTime;
result.geometryTime = pairBuildTime+result.gridBuildTime;
if opt.visualise
    plot_stokes_two_body_mfs(result);
end
end
