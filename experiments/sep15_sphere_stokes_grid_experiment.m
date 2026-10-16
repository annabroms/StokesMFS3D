function result = sep15_sphere_stokes_grid_experiment(radius,delta,opt)
%SEP15_SPHERE_STOKES_GRID_EXPERIMENT Basic/enhanced two-sphere resistance.
%
%   result = sep15_sphere_stokes_grid_experiment
%   opt.enhancingFamilies = 'SRTD';  % also accepts ["S","R","T","D"]
%   opt.enhancingSourceCount = 400;
%   opt.sourceBeta = 0.5;
%   opt.visualise = false;
%   result = sep15_sphere_stokes_grid_experiment(0.63,0.01,opt);
%
% This is the Stokes counterpart of the sphere Laplace grid experiment,
% NOT a scalar capacitance problem.  The spheres translate towards each
% other with velocities [1,0,0] and [-1,0,0], with zero rotation by default.
% Viscosity is 1 and velocity decays at infinity.  Only the chopped spheroid
% is used for enhancement.  Geometry and grids reuse the existing helpers.
% Figures show body 1 only; the nearest accumulation point is black.
% No files are saved.  Default dense SVDs are substantially larger than
% their Laplace counterparts; begin with smaller counts for a quick test.

if nargin < 1 || isempty(radius)
    radius = 0.63;
end
if nargin < 2 || isempty(delta)
    delta = 0.01;
end
if nargin < 3 || isempty(opt)
    opt = struct();
end

% MAIN STOKES CONTROLS: edit the enhancing combination here.
% S = Stokeslet; R = rotlet; T = stresslet with the parent surface normal;
% D = potential dipole.  Basic proxy nodes ALWAYS carry Stokeslets only.
% Each selected family adds THREE unknowns per retained enhancing node.
stokesControls.enhancingFamilies = 'SRTD';
stokesControls.translations = [1,0,0;-1,0,0];
stokesControls.angularVelocities = zeros(2,3);

% R/T columns are scaled by L and D columns by L^2.  With this fixed scale,
% all fitted coefficients have force units.  Changing L changes their
% max norm and the truncated-SVD solution; physical strengths are returned
% separately.  This is NOT a normalization to make coefficients small.
stokesControls.singularityLengthScale = radius;
stokesControls.evaluationBatchSize = 256;  % memory control, not grid density

% MAIN GRID CONTROLS, matching the sphere Laplace setup.  Counts are per
% body.  The repository picks the largest spherical design <= each request.
% Enhancing count is for the WHOLE parent spheroid BEFORE chopping.  The
% retained count, number of unknowns and total extra collocation are printed.
gridControls.basicGridDistribution = "sphericaldesign";
gridControls.basicSourceCount = 700;
gridControls.basicCollocationCount = 700*1.2;
gridControls.enhancingSourceDistribution = "sphericaldesign";
gridControls.enhancingSourceShape = "spheroid";
gridControls.enhancingSourceCount = 640;
gridControls.fineSourceResolution = [32,20]; % tensor alternative: [around, axial]
gridControls.sourceChopOffset = 0; % signed cap depth beyond the proxy cut
gridControls.sourceOffset = 0.125;  % proxy radius = radius-sourceOffset

% EASY BETA CONTROL: fraction of the accumulation-to-boundary distance by
% which the contact-facing spheroid tip moves outwards.  Require 0<beta<1.
% The continuous proxy radius may require further stretching.  Sampled
% proxy nodes need not be enclosed.  Both caps are exact reflected copies.
gridControls.sourceBeta = 0.3;
gridControls.sourceWidthScale = 1;  % scales transverse axes, not tip position

% Patch rings have rho_j=angle*expm1(alpha*j/L)/expm1(alpha).
gridControls.collocationPatchLayers = 18; % L: more rings/refinement
gridControls.collocationExponent = 3;     % alpha: stronger contact clustering
gridControls.collocationPatchAngle = 0.5; % half-angle in radians
gridControls.collocationMinAzimuth = 8;   % azimuthal count floor on each ring
gridControls.collocationAzimuthOversampling = 1; % multiplies adaptive counts
% For more collocation equations, increase layers and/or oversampling.
% The system must have more equations than ALL selected family unknowns.
gridControls.basicSvdCutoff = 1e-13;
gridControls.enhancedSvdCutoff = 1e-13;
gridControls.testResolution = [128,80]; % independent full-body grid only
gridControls.testPatchLayers = 48; % independent contact evaluation rings
gridControls.testPatchAzimuthOversampling = 2;
gridControls.visualise = true;
% true shows exact system sizes and waits for y/yes before kernels/SVDs.
% false bypasses confirmation only for an already-approved automated run.
gridControls.confirmBeforeSolve = true;

% Normalize options once at the entry point; opt overrides editable values.
validateattributes(radius,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(delta,{'numeric'},{'real','finite','scalar','positive'});
[stokesControls,opt] = complete_stokes_two_body_controls(opt,stokesControls);
gridNames = fieldnames(gridControls);
for k = 1:numel(gridNames)
    name = gridNames{k};
    if ~isfield(opt,name)
        opt.(name) = gridControls.(name);
    end
end

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory,fullfile(projectDirectory,'common','dense'), ...
    fullfile(projectDirectory,'geometry','spheres'), ...
    fullfile(projectDirectory,'geometry','ellipsoids'), ...
    fullfile(projectDirectory,'laplace'));
opt = lap3d_complete_grid_options(opt);
assert(opt.sourceOffset < radius,'sourceOffset must be smaller than radius.');
assert(opt.enhancingSourceShape == "spheroid", ...
    'This experiment uses a chopped spheroid only.');

pairStart = tic;
pair = lap3d_sphere_pair(radius,delta);
pairBuildTime = toc(pairStart);
result = run_stokes_two_body_grid_experiment(pair,"sphere",opt,stokesControls);
if isempty(result)
    return % declined: no numerical solutions or figures to process
end
result.sphereRadius = radius;
result.pairBuildTime = pairBuildTime;
result.geometryTime = pairBuildTime+result.gridBuildTime;
if opt.visualise
    plot_stokes_two_body_mfs(result);
end
end
