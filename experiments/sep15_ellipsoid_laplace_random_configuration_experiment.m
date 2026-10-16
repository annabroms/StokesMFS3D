function summary = sep15_ellipsoid_laplace_random_configuration_experiment(nruns,delta,opt)
%SEP15_ELLIPSOID_LAPLACE_RANDOM_CONFIGURATION_EXPERIMENT Fixed-gap random pairs.
%   summary = sep15_ellipsoid_laplace_random_configuration_experiment;
%   summary = sep15_ellipsoid_laplace_random_configuration_experiment(10,0.01);
%   opt.randomSeed = 7; opt.basicCollocationResolution = [80,60];
%   summary = sep15_ellipsoid_laplace_random_configuration_experiment(10,0.01,opt);
% Each pair uses independent random rotations and a random gap direction.
% Facing support points certify the exact closest gap delta.  There is no
% random global translation: the gap midpoint is at the origin for plotting.
% Reuse the reflection accumulation model, basic/enhanced grids, and dense
% Laplace solves with potentials [1;-1].  No helpers or kernels are copied.
% Accumulation locations are osculating-sphere predictions, not an exact
% ellipsoid singularity construction.  The existing prolate focal-segment
% diagnostic does not certify the focal set of a triaxial ellipsoid.
% One configuration figure per attempted solve, plus a sweep summary figure.
% Relative residual = ||u-u_prescribed||_inf / ||[1;-1]||_inf; coefficient
% size = max(abs(c)).  Test residuals use the independent finer check grid.
% Residuals and coefficients use max norms over BOTH particles.  Failed cap
% checks are recorded, not repaired or resampled.  No files are saved.

% NUMBER OF CONFIGURATIONS: positive integer; change 10 here for a longer run.
if nargin < 1 || isempty(nruns)
    nruns = 20;
end
% FIXED GAP: positive closest surface distance, in physical length units.
if nargin < 2 || isempty(delta)
    delta = 0.01;
end
% OPTIONS: opt.FIELD overrides controls below; opt.semiAxes overrides axes.
if nargin < 3 || isempty(opt)
    opt = struct();
end
% GEOMETRY: positive axes [a;b;c]; equal axes give the sphere limit.
semiAxes = [1;0.9;0.7];
% RANDOM SEED: nonnegative integer; change for a different repeatable sample.
% Restores the caller's RNG state on exit.  Also accepts opt.randomSeed.
randomSeed = 1;

% GRID DEFAULTS copied from the current sep09 ellipsoid Laplace experiment.
% This is a snapshot: later edits in that experiment are not auto-imported.
% ENHANCING BASIS: "SLP" = monopoles; "SLP+DLP" adds independently fitted
% parent-normal dipoles on the cap.  Basic proxy nodes always stay SLP.
% Combined mode adds one unknown per enhancing node in every configuration.
gridControls.enhancingLaplaceBasis = "SLP";
% Positive physical length L for DLP columns; dipole moment = L*d*n.
% Keep fixed when comparing fitted coefficient max norms in charge units.
gridControls.laplaceDipoleLengthScale = 1;
% Independent switches for BOTH basic/enhanced solves in every configuration.
% true: normalize columns in the 2-norm AFTER row weighting; undo on coefficients.
gridControls.columnScaling = false;
% true: weight equations AND data by sqrt(local area on the COMPLETE grid),
% accounting for merged patch density; this fits an area-weighted residual.
gridControls.surfaceElementPreconditioning = false;
% SVD cutoffs apply to the scaled matrix. Coefficient magnitudes and max
% residuals are always reported in physical, unweighted units (not scaled).
% BASIC DISTRIBUTION: quadrature uses the two resolutions below; designs
% use basicSourceCount/basicCollocationCount instead (counts PER BODY).
gridControls.basicGridDistribution = "ellipsoidquadrature";
% Proxy source [azimuth scale, GL layers]; integers >= 4, rounded to even.
gridControls.coarseSourceResolution = [35,28];
% Independent BASIC boundary [azimuth scale, GL layers]; integers >= 4.
% More azimuths resolve rings, more layers resolve pole-to-pole spacing.
% Ring counts taper near poles, so total is not Nu*Nv.  [] uses default.
gridControls.basicCollocationResolution = [round(1.5*35),round(1.5*28)];
% Requested design sources per body, integer >= 3; design mode only.
gridControls.basicSourceCount = 700;
% Requested design boundary nodes per body, integer >= 3; design mode only.
gridControls.basicCollocationCount = 1000;

% ENHANCING SHAPE: "spheroid" is the chopped stretched prolate parent;
% "sphere" instead uses a positive sourceSphereRadius around accumulation.
gridControls.enhancingSourceShape = "spheroid";
% Leave [] for spheroids; set a positive enclosing radius for sphere mode.
gridControls.sourceSphereRadius = [];
% DISTRIBUTION: "quasiuniform" tapers ring counts; "tensor" fixes counts;
% "sphericaldesign" affine-maps a whole-parent repository design.
gridControls.enhancingSourceDistribution = "quasiuniform";
% Enhancing [azimuth scale, axial GL layers] BEFORE chopping, not basic.
% Tensor fixes the first entry per ring; quasiuniform tapers it near tips.
gridControls.fineSourceResolution = [42,30];
% Extra ring azimuths, integer >= 0; quasiuniform only.
gridControls.enhancingSourceExtraAzimuthNodes = 0;
% Small-ring azimuthal floor, integer >= 2; quasiuniform only.
gridControls.enhancingSourceMinAzimuthNodes = 8;
% Requested whole-parent design count per body, integer >= 3; designs only.
gridControls.enhancingSourceCount = 1200;
% Positive inward NORMAL proxy displacement in physical length units.
gridControls.sourceOffset = 0.25;
% Requested tip fraction towards radial boundary, 0 < beta < 1.
% Nonspherical proxy enclosure can increase the effective beta.
gridControls.sourceBeta = 0.2;
% Positive multiplier of both transverse spheroid axes; 1 = unchanged.
gridControls.sourceWidthScale = 0.8;
% Signed physical depth beyond proxy cut; positive overlaps the proxy.
% Require sourceOffset+sourceChopOffset > 0; changes cut, not parent tip.
gridControls.sourceChopOffset = 0.05;
% false keeps duplicate enhancing sources; collocation is always deduplicated.
gridControls.removeDuplicates = false;

% COLLOCATION PATCH: positive integer count of graded rings around contact.
gridControls.collocationPatchLayers = 25;
% Positive exponent; larger clusters rings more strongly at contact.
gridControls.collocationExponent = 5;
% Patch half-angle in parameter-sphere radians, between 0 and pi.
gridControls.collocationPatchAngle = 0.5;
% Minimum azimuths per ring, integer >= 4; adaptive counts become even.
gridControls.collocationMinAzimuth = 8;
% Positive adaptive azimuth multiplier; 2 roughly doubles ring counts.
gridControls.collocationAzimuthOversampling = 1;

% INDEPENDENT CHECK ONLY: whole-body [azimuth scale, GL layers], >= 4.
% These nodes do not enter the fit; refine until max residual stabilizes.
gridControls.testResolution = [128,80];
% Independent graded contact rings, integer >= 0; 0 disables them.
gridControls.testPatchLayers = 48;
% Positive azimuth multiplier for independent contact checking only.
gridControls.testPatchAzimuthOversampling = 2;
% Separate relative SVD cutoffs in (0,1); larger discards more modes.
gridControls.basicSvdCutoff = 1e-12;
gridControls.enhancedSvdCutoff = 1e-12;
% true makes one figure per configuration and the final comparison figure.
gridControls.visualise = true;
% true asks y/yes for EACH configuration's exact basic/enhanced matrix sizes.
% No kernel or SVD starts before consent; declining stops the remaining runs.
% false is a bypass only for an explicitly pre-approved automated sweep.
gridControls.confirmBeforeSolve = false;

validateattributes(nruns,{'numeric'},{'real','finite','scalar','integer','positive'});
assert(isstruct(opt) && isscalar(opt),'opt must be a scalar struct.');
if isfield(opt,'randomSeed')
    randomSeed = opt.randomSeed;
    opt = rmfield(opt,'randomSeed');
end
validateattributes(randomSeed,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',0,'<=',2^32-1});

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory,fullfile(projectDirectory,'geometry','ellipsoids'), ...
    fullfile(projectDirectory,'geometry','spheres'),fullfile(projectDirectory,'laplace'));
originalRng = rng;
rngCleanup = onCleanup(@() rng(originalRng));
rng(randomSeed,'twister');

% Invariant: randomized geometry changes neither the kernels, data, nor
% prescribed grid controls.  Columns below always mean [basic,enhanced].
configurations = cell(nruns,1);
runs = cell(nruns,1);
status = repmat("not attempted",nruns,1);
failureMessages = strings(nruns,1);
gaps = nan(nruns,1);
collocationMaxRelativeResidual = nan(nruns,2);
evaluationMaxRelativeResidual = nan(nruns,2);
nearContactMaxRelativeResidual = nan(nruns,2);
coefficientMax = nan(nruns,2);
cancelled = false;

% opt is completed once on the first call to the shared preparation routine.
for runIndex = 1:nruns
    fprintf('\nRandom configuration %d/%d, seed %d\n',runIndex,nruns,randomSeed);
    if runIndex == 1
        [pair,opt,pairBuildTime] = prepare_ellipsoid_grid_experiment( ...
            "random",delta,semiAxes,gridControls,opt);
        semiAxes = pair.E1.axes;
        if opt.visualise
            [meshX,meshY,meshZ] = ellipsoid_matlab(0,0,0, ...
                semiAxes(1),semiAxes(2),semiAxes(3),40);
            meshPoints = [meshX(:),meshY(:),meshZ(:)].';
        end
    else
        pairStart = tic;
        pair = sep09_ellipsoid_reflection_accumulation_experiment( ...
            "random",delta,false,semiAxes);
        pairBuildTime = toc(pairStart);
    end
    pair.label = sprintf('Random ellipsoid configuration %d',runIndex);
    configurations{runIndex} = pair;
    gaps(runIndex) = pair.gap;
    assert(abs(pair.gap-delta) < 200*eps*max(1,delta), ...
        'Random configuration did not preserve the requested gap.');

    try
        runResult = lap3d_run_two_body_grid_experiment( ...
            pair,sprintf('random_%02d',runIndex),opt);
        if isempty(runResult)
            status(runIndex) = "cancelled";
            cancelled = true;
            break
        end
        runResult.options.semiAxes = semiAxes;
        runResult.pairBuildTime = pairBuildTime;
        runResult.geometryTime = pairBuildTime+runResult.gridBuildTime;
        runs{runIndex} = runResult;
        solutions = {runResult.basic.solution,runResult.enhanced.solution};
        for method = 1:2
            collocationMaxRelativeResidual(runIndex,method) = ...
                solutions{method}.collocationMaxResidual;
            evaluationMaxRelativeResidual(runIndex,method) = ...
                solutions{method}.testMaxResidual;
            nearContactMaxRelativeResidual(runIndex,method) = ...
                solutions{method}.testNearContactMaxResidual;
            coefficientMax(runIndex,method) = solutions{method}.coefficientNorm;
        end
        status(runIndex) = "solved";
    catch exception
        % Keep failed geometries in the sample: no silent parameter changes.
        status(runIndex) = "failed";
        failureMessages(runIndex) = string(exception.identifier)+": "+ ...
            string(exception.message);
        warning('random_ellipsoid_laplace:ConfigurationFailed', ...
            'Configuration %d failed: %s',runIndex,failureMessages(runIndex));
    end

    if opt.visualise
        figure('Name',sprintf('Random Laplace configuration %d',runIndex));
        ax = axes;
        hold(ax,'on')
        bodies = {pair.E1,pair.E2};
        contacts = {pair.p1,pair.p2};
        accumulation = {pair.Pminus,pair.Pplus};
        particleColors = [0.20,0.45,0.85;0.85,0.40,0.20];
        for body = 1:2
            worldPoints = bodies{body}.c+bodies{body}.Q*meshPoints;
            surf(ax,reshape(worldPoints(1,:),size(meshX)), ...
                reshape(worldPoints(2,:),size(meshY)), ...
                reshape(worldPoints(3,:),size(meshZ)), ...
                'FaceColor',particleColors(body,:),'FaceAlpha',0.45, ...
                'EdgeColor','none');
            [~,nearestIndex] = min(vecnorm(accumulation{body}-contacts{body},2,1));
            nearestPoint = accumulation{body}(:,nearestIndex);
            scatter3(ax,nearestPoint(1),nearestPoint(2),nearestPoint(3), ...
                65,'k','filled');
        end
        plot3(ax,[pair.p1(1),pair.p2(1)],[pair.p1(2),pair.p2(2)], ...
            [pair.p1(3),pair.p2(3)],'k-o','MarkerFaceColor',[1,0.85,0], ...
            'LineWidth',1.2);
        axis(ax,'equal')
        grid(ax,'on')
        view(ax,3)
        xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
        configurationTitle = sprintf('Configuration %d/%d, gap %.6g, seed %d, %s', ...
            runIndex,nruns,pair.gap,randomSeed,opt.enhancingLaplaceBasis);
        if status(runIndex) == "solved"
            title(ax,{configurationTitle, ...
                sprintf('Basic: max rel. residual col/test %.2e / %.2e; max |c| %.2e', ...
                collocationMaxRelativeResidual(runIndex,1), ...
                evaluationMaxRelativeResidual(runIndex,1),coefficientMax(runIndex,1)), ...
                sprintf('Enhanced: max rel. residual col/test %.2e / %.2e; max |c| %.2e', ...
                collocationMaxRelativeResidual(runIndex,2), ...
                evaluationMaxRelativeResidual(runIndex,2),coefficientMax(runIndex,2)), ...
                'All maxima over both particles; nearest accumulation points in black'});
        else
            title(ax,{configurationTitle,'FAILED: see summary.failureMessages'}, ...
                'Color',[0.8,0,0]);
        end
        drawnow
    end
end

runNumbers = (1:nruns).';
summary.nruns = nruns;
summary.delta = delta;
summary.semiAxes = semiAxes;
summary.randomSeed = randomSeed;
summary.options = opt;
summary.methodNames = ["basic","enhanced"];
summary.configurations = configurations;
summary.runs = runs;
summary.status = status;
summary.failureMessages = failureMessages;
summary.cancelled = cancelled;
summary.collocationMaxRelativeResidual = collocationMaxRelativeResidual;
summary.evaluationMaxRelativeResidual = evaluationMaxRelativeResidual;
summary.nearContactMaxRelativeResidual = nearContactMaxRelativeResidual;
summary.coefficientMax = coefficientMax;
summary.table = table(runNumbers,gaps,status, ...
    collocationMaxRelativeResidual(:,1),collocationMaxRelativeResidual(:,2), ...
    evaluationMaxRelativeResidual(:,1),evaluationMaxRelativeResidual(:,2), ...
    coefficientMax(:,1),coefficientMax(:,2),'VariableNames', ...
    {'Run','Gap','Status','BasicCollocationResidual','EnhancedCollocationResidual', ...
    'BasicEvaluationResidual','EnhancedEvaluationResidual', ...
    'BasicCoefficientMax','EnhancedCoefficientMax'});
fprintf('\nRandom configurations: %d solved, %d failed, %d cancelled.\n', ...
    nnz(status == "solved"),nnz(status == "failed"),nnz(status == "cancelled"));
disp(summary.table)

if opt.visualise && any(status == "solved")
    figure('Name','Random ellipsoid Laplace configuration comparison');
    layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(layout);
    semilogy(ax,runNumbers,evaluationMaxRelativeResidual,'o-');
    hold(ax,'on')
    semilogy(ax,runNumbers,collocationMaxRelativeResidual,'x--');
    xlabel(ax,'configuration'); ylabel(ax,'max relative residual');
    grid(ax,'on'); xticks(ax,runNumbers);
    legend(ax,{'basic evaluation','enhanced evaluation', ...
        'basic collocation','enhanced collocation'},'Location','best');
    title(ax,'Max relative residual over both particles')
    ax = nexttile(layout);
    semilogy(ax,runNumbers,coefficientMax,'o-');
    xlabel(ax,'configuration'); ylabel(ax,'max |c_j|');
    grid(ax,'on'); xticks(ax,runNumbers);
    legend(ax,{'basic','enhanced'},'Location','best');
    title(ax,'Coefficient max norm over all sources')
    title(layout,sprintf('Fixed gap %.3g, seed %d; %d/%d configurations solved', ...
        delta,randomSeed,nnz(status == "solved"),nruns))
end
end
