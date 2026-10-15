% OCT15_2026_TWO_DROP_AXIAL_PRESSURE_GAP_SWEEP
% Two identical active drops have their common direction e along the line
% of centres. Test whether BOTH translations point toward the other drop,
% and measure exterior pressure at the midpoint. The fixed-sphere model
% leaves interior pressure gauges undetermined, so the plane plot shows
% exterior pressure only. Run this script to perform the gap sweep.

clear; close all

%% Physical and numerical parameters
dropRadius = 1;
innerProxyRadius = 0.6;
outerProxyRadius = 2;
S1 = 2;                           % same active-stress strength on both drops
interiorViscosity = 10;            % ratio to exterior viscosity mu_0 = 1
activeDirection = [1 0 0];       % same e_i, along the centre-to-centre axis
gapList = [1 0.5 0.2 0.1 0.05 0.02 0.01 0.001];
representativeGap = 0.01;
requestedN = 600;                 % sources per proxy surface and per drop
collocationToSourceRatio = 1.2;
validationResolution = 1000;     % minimum independent points per drop
oneBodySvdTolerance = 1e-13;
krylovTolerance = 1e-10;
maxIterations = 150;              % unrestarted GMRES iteration limit
tractionRowWeight = 1;
scale_columns = false;
interactionBatchSize = 128;       % targets per temporary GMRES kernel block
validationBatchSize = 128;        % targets per temporary validation block
pressureBatchSize = 128;          % targets per temporary pressure kernel block
planeGridResolution = [401 241]; % [x,z] points for pressure cross-section
planePaddingRadii = 1.5;
view_solution = true;

%% Paths and spherical-design grids
experimentDirectory = fileparts(mfilename('fullpath'));
dropsDirectory = fileparts(experimentDirectory);
repositoryDirectory = fileparts(dropsDirectory);
addpath(dropsDirectory,fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'common','eval'), ...
    fullfile(repositoryDirectory,'common','gmres'), ...
    fullfile(repositoryDirectory,'geometry','spheres'));

gapList = sort(gapList);
if ~ismember(representativeGap,gapList)
    error('The representative gap must appear in gapList.');
end
catalogueFile = fullfile(repositoryDirectory,'geometry','spheres', ...
    'spherical_designs','filelist.txt');
catalogueID = fopen(catalogueFile,'r');
if catalogueID < 0
    error('Cannot read the spherical-design catalogue.');
end
catalogueData = textscan(catalogueID,'%*[^.].%d');
fclose(catalogueID);
availableCounts = catalogueData{1};
[sourceDirections,~] = get_sphdesign(requestedN);
N = size(sourceDirections,1);
minimumM = ceil(collocationToSourceRatio*N);
collocationIndex = find(availableCounts>=minimumM,1);
if isempty(collocationIndex)
    error('No spherical design supplies the required collocation nodes.');
end
[collocationDirections,~] = get_sphdesign(availableCounts(collocationIndex));
M = size(collocationDirections,1);
validationIndex = find(availableCounts>=max(validationResolution,M+1),1);
if isempty(validationIndex)
    error('No spherical design supplies the independent validation nodes.');
end
[validationDirections,~] = get_sphdesign(availableCounts(validationIndex));
% Rotate the independent design so that it has no collocation points.
azimuth = 0.37;
elevation = 0.23;
rotationZ = [cos(azimuth) -sin(azimuth) 0; ...
             sin(azimuth) cos(azimuth) 0;0 0 1];
rotationY = [cos(elevation) 0 sin(elevation); ...
             0 1 0;-sin(elevation) 0 cos(elevation)];
validationDirections = (rotationZ*rotationY*validationDirections.').';
validationCount = size(validationDirections,1);
fprintf('Per drop: requested N=%d, actual N=%d; required M>=%d, actual M=%d (M/N=%.3f).\n', ...
    requestedN,N,minimumM,M,M/N);
fprintf('Independent validation: requested >=%d, actual %d points per drop.\n', ...
    validationResolution,validationCount);

%% Solve settings and physical residual scales
activeDirection = activeDirection/norm(activeDirection);
activeDirections = repmat(activeDirection,2,1);
S1ByDrop = [S1;S1];
if S1 == 0
    velocityScale = 1;
    tractionScale = 1;
    residualScaleLabel = 'unit-scaled';
else
    velocityScale = abs(S1)*dropRadius; % mu_0 = 1
    tractionScale = abs(S1);
    residualScaleLabel = 'relative';
end
solverOpt.interiorViscosity = interiorViscosity;
solverOpt.completionScale = 1/(8*pi*dropRadius);
solverOpt.velocityScale = velocityScale;
solverOpt.tractionScale = tractionScale;
solverOpt.tractionRowWeight = tractionRowWeight;
solverOpt.oneBodySvdTolerance = oneBodySvdTolerance;
solverOpt.krylovTolerance = krylovTolerance;
solverOpt.maxIterations = maxIterations;
solverOpt.interactionBatchSize = interactionBatchSize;
solverOpt.scale_columns = scale_columns;

%% Gap sweep: motion, midpoint pressure, and independent-grid residuals
caseCount = numel(gapList);
translationX = zeros(caseCount,2);
transverseSpeed = zeros(caseCount,1);
gapRate = zeros(caseCount,1);
midpointPressure = zeros(caseCount,1);
relativeBoundaryMaxima = zeros(caseCount,3);
gmresIterations = zeros(caseCount,1);
gmresResiduals = zeros(caseCount,1);
originalResiduals = zeros(caseCount,1);
setupTimes = zeros(caseCount,1);
solveTimes = zeros(caseCount,1);
oneBodyRanks = zeros(caseCount,1);
newFactorCounts = zeros(caseCount,1);
factorBank = {};

for caseIndex = 1:caseCount
    delta = gapList(caseIndex);
    centreDistance = 2*dropRadius+delta;
    centres = [-centreDistance/2 0 0;centreDistance/2 0 0];
    innerSources = cell(1,2);
    outerSources = cell(1,2);
    surfacePoints = cell(1,2);
    normals = cell(1,2);
    validationPoints = cell(1,2);
    validationNormals = cell(1,2);
    for body = 1:2
        innerSources{body} = centres(body,:)+ ...
            innerProxyRadius*sourceDirections;
        outerSources{body} = centres(body,:)+ ...
            outerProxyRadius*sourceDirections;
        surfacePoints{body} = centres(body,:)+ ...
            dropRadius*collocationDirections;
        normals{body} = collocationDirections;
        validationPoints{body} = centres(body,:)+ ...
            dropRadius*validationDirections;
        validationNormals{body} = validationDirections;
    end

    % Each body has [alpha(3N);gamma(3N)] coefficients. The recovered
    % exterior strengths lambda=P*alpha are force free, and U comes from
    % alpha's completion. Identical translated drops reuse one SVD factor.
    [translations,~,gamma,lambda,diagnostics,factorBank] = ...
        solve_drops_preconditioned(innerSources,outerSources, ...
        surfacePoints,normals,S1ByDrop,activeDirections, ...
        solverOpt,factorBank);
    translationX(caseIndex,:) = translations(:,1).';
    transverseSpeed(caseIndex) = max(vecnorm(translations(:,2:3),2,2));
    gapRate(caseIndex) = translations(2,1)-translations(1,1);

    % The midpoint is x=(0,0,0), strictly in the exterior for delta>0.
    % p_0 tends to zero at infinity; no interior pressure gauge enters.
    midpointPressure(caseIndex) = drop_exterior_pressure( ...
        [0 0 0],innerSources,lambda,pressureBatchSize);

    maximumPhysicalErrors = drop_validate_surface_conditions( ...
        lambda,gamma,translations,innerSources,outerSources, ...
        validationPoints,validationNormals,activeDirections, ...
        S1ByDrop,interiorViscosity,validationBatchSize);
    relativeBoundaryMaxima(caseIndex,:) = max( ...
        maximumPhysicalErrors./[velocityScale velocityScale tractionScale], ...
        [],1);
    gmresIterations(caseIndex) = diagnostics.iterations;
    gmresResiduals(caseIndex) = diagnostics.regularizedResidual;
    originalResiduals(caseIndex) = diagnostics.originalResidual;
    setupTimes(caseIndex) = diagnostics.setupTime;
    solveTimes(caseIndex) = diagnostics.solveTime;
    oneBodyRanks(caseIndex) = diagnostics.retainedRanks(1);
    newFactorCounts(caseIndex) = diagnostics.newFactorCount;

    fprintf('delta=%7.4g: U1x=% .5e, U2x=% .5e, d(delta)/dt=% .5e, p_mid=% .5e; GMRES %d, rel BC [%.2e %.2e %.2e].\n', ...
        delta,translationX(caseIndex,1),translationX(caseIndex,2), ...
        gapRate(caseIndex),midpointPressure(caseIndex), ...
        gmresIterations(caseIndex),relativeBoundaryMaxima(caseIndex,:));

    if delta == representativeGap
        representativeCentres = centres;
        representativeSources = innerSources;
        representativeStrengths = lambda;
    end
end

%% Signed motion test and pressure summary
% With c_1 to the left of c_2, both drops move toward one another exactly
% when U_1x>0 and U_2x<0. A negative d(delta)/dt means the gap closes.
leftTowardSpeed = translationX(:,1);
rightTowardSpeed = -translationX(:,2);
bothPointToward = leftTowardSpeed>0 & rightTowardSpeed>0;
fprintf('\nMotion toward one another at %d of %d gaps; gap closes at %d of %d gaps.\n', ...
    nnz(bothPointToward),caseCount,nnz(gapRate<0),caseCount);
fprintf('Boundary residuals are %s: velocity scale %.3g, traction scale %.3g.\n', ...
    residualScaleLabel,velocityScale,tractionScale);
fprintf('Pressure is physical exterior pressure with p_0 -> 0 at infinity; midpoint magnitude is |p_mid|.\n');
disp(table(gapList(:),leftTowardSpeed,rightTowardSpeed,gapRate, ...
    midpointPressure,abs(midpointPressure),bothPointToward, ...
    transverseSpeed,gmresIterations,gmresResiduals,originalResiduals, ...
    relativeBoundaryMaxima(:,1),relativeBoundaryMaxima(:,2), ...
    relativeBoundaryMaxima(:,3),oneBodyRanks,newFactorCounts, ...
    setupTimes,solveTimes, ...
    'VariableNames',{'Gap','U1Toward','U2Toward','GapRate', ...
    'MidpointPressure','PressureMagnitude','BothPointToward', ...
    'MaxTransverseSpeed','GMRESIterations','GMRESResidual', ...
    'OriginalResidual','ContinuityRel','NormalRel','TractionRel', ...
    'OneBodyRank','NewFactors','SetupSeconds','SolveSeconds'}));

%% Gap plots and exterior pressure in the x-z plane
if view_solution
    figure('Name','Axial drop motion and midpoint pressure versus gap', ...
        'Color','w');
    tiledlayout(1,2);
    nexttile;
    plot(gapList,leftTowardSpeed,'o-','LineWidth',1.3); hold on;
    plot(gapList,rightTowardSpeed,'s-','LineWidth',1.3);
    yline(0,'k:');
    set(gca,'XScale','log');
    axis tight; grid on;
    xlabel('$\delta$','Interpreter','latex');
    ylabel('Speed toward the other drop','Interpreter','latex');
    title('Positive values mean approach','Interpreter','latex');
    legend({'$U_{1x}$','$-U_{2x}$'}, ...
        'Interpreter','latex','Location','best');

    nexttile;
    plot(gapList,abs(midpointPressure),'o-','LineWidth',1.3);
    set(gca,'XScale','log');
    axis tight; grid on;
    xlabel('$\delta$','Interpreter','latex');
    ylabel('$|p_0(0,0,0)|$','Interpreter','latex');
    title('Exterior midpoint pressure magnitude','Interpreter','latex');

    centreDistance = 2*dropRadius+representativeGap;
    xExtent = centreDistance/2+dropRadius+ ...
        planePaddingRadii*dropRadius;
    zExtent = (1+planePaddingRadii)*dropRadius;
    xGrid = linspace(-xExtent,xExtent,planeGridResolution(1));
    zGrid = linspace(-zExtent,zExtent,planeGridResolution(2));
    [gridX,gridZ] = meshgrid(xGrid,zGrid);
    planeTargets = [gridX(:),zeros(numel(gridX),1),gridZ(:)];
    exterior = true(size(planeTargets,1),1);
    for body = 1:2
        exterior = exterior & ...
            sum((planeTargets-representativeCentres(body,:)).^2,2) > ...
            dropRadius^2;
    end
    planePressure = nan(size(planeTargets,1),1);
    planePressure(exterior) = drop_exterior_pressure( ...
        planeTargets(exterior,:),representativeSources, ...
        representativeStrengths,pressureBatchSize);

    figure('Name','Signed exterior pressure through two active drops', ...
        'Color','w');
    pressureImage = imagesc(xGrid,zGrid,reshape(planePressure,size(gridX)));
    set(pressureImage,'AlphaData',reshape(exterior,size(gridX)));
    set(gca,'YDir','normal');
    hold on;
    theta = linspace(0,2*pi,301);
    for body = 1:2
        plot(representativeCentres(body,1)+dropRadius*cos(theta), ...
            dropRadius*sin(theta),'k-','LineWidth',1.4);
    end
    plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',4);
    maximumPressure = max(abs(planePressure(exterior)));
    if maximumPressure > 0
        clim([-maximumPressure maximumPressure]);
    end
    % A gray zero distinguishes near-zero exterior pressure from the white
    % masked drop interiors; blue and red retain the pressure sign.
    colourPositions = linspace(0,1,257).';
    pressureColours = interp1([0 0.5 1], ...
        [0.10 0.27 0.60;0.70 0.70 0.70;0.65 0.07 0.12], ...
        colourPositions);
    colormap(pressureColours);
    pressureColourbar = colorbar;
    pressureColourbar.Label.String = '$p_0$';
    pressureColourbar.Label.Interpreter = 'latex';
    axis equal tight;
    xlim([-xExtent xExtent]);
    ylim([-zExtent zExtent]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    title(sprintf('Exterior pressure at $\\delta=%.3g$; interiors masked', ...
        representativeGap),'Interpreter','latex');
end
