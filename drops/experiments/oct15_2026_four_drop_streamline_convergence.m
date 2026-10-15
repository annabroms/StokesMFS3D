% OCT15_2026_FOUR_DROP_STREAMLINE_CONVERGENCE  Compare two MFS resolutions.
% Solve exactly the same four-drop problem with requested N=400 and 800.
% Evaluate both fields at the SAME x-z points, trace from the SAME seeds,
% and compare path positions at equal arc length from each seed. The runs
% are sequential, so their large one-body SVD factors are not held at once.
% The instantaneous lab-frame lines may cross moving interfaces; they are
% streamlines of the fixed-time field, not fluid particle trajectories.

clear;

%% Problem and solver parameters
P = 4;
radii = [0.8;1.0;1.2;0.9];
interiorViscosities = [1.5;2;3;4]; % exterior viscosity mu_0 = 1
gap = 1;
innerProxyRadiusRatio = 0.6;
outerProxyRadiusRatio = 2;
activeDirection = [0 0 1];
randomSeed = 15;
S1Minimum = 0.5;
S1Maximum = 1.5;
requestedResolutions = [400 800]; % source nodes per proxy surface
collocationToSourceRatio = 1.2;
validationResolution = 1000;
oneBodySvdTolerance = 1e-13;
krylovTolerance = 1e-10;
maxIterations = 150;             % total unrestarted GMRES iterations
tractionRowWeight = 1;
scale_columns = false;
interactionBatchSize = 128;      % targets per temporary interaction block
validationBatchSize = 128;       % targets per temporary validation block
view_solution = true;            % draw overlays and field-difference maps

%% Shared x-z field samples and streamline settings
plotPaddingInRadii = 1.5;
plotGridCountX = 161;
plotGridCountZ = 101;
globalSeedsPerDrop = 4;
interiorSeedsPerDrop = 6;
globalSeedRadiusRatio = 1.25;
interiorSeedRadiusRatio = 0.45;
maxArcLength = 20*max(radii);     % in each direction from a seed
maxStep = 0.04*min(radii);        % ode45 maximum arclength step
minimumSpeed = 1e-6;             % stagnation stop threshold
odeRelativeTolerance = 1e-6;
odeAbsoluteTolerance = 1e-8*min(radii);
streamlineComparisonSamples = 101; % matched arclength samples/direction

%% Repository routines and spherical-design catalogue
experimentDirectory = fileparts(mfilename('fullpath'));
dropsDirectory = fileparts(experimentDirectory);
repositoryDirectory = fileparts(dropsDirectory);
addpath(dropsDirectory,fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'common','gmres'), ...
    fullfile(repositoryDirectory,'geometry','spheres'));
catalogueFile = fullfile(repositoryDirectory,'geometry','spheres', ...
    'spherical_designs','filelist.txt');
catalogueID = fopen(catalogueFile,'r');
if catalogueID < 0
    error('Cannot read the spherical-design catalogue.');
end
catalogueData = textscan(catalogueID,'%*[^.].%d');
fclose(catalogueID);
availableCounts = catalogueData{1};

%% Fixed geometry, active stresses, validation points, and comparison seeds
centreX = zeros(P,1);
for body = 2:P
    centreX(body) = centreX(body-1)+radii(body-1)+gap+radii(body);
end
centreX = centreX-mean(centreX);
centres = [centreX zeros(P,2)];
rng(randomSeed,'twister');
S1 = S1Minimum+(S1Maximum-S1Minimum)*randperm(P).'/(P+1);
activeDirections = repmat(activeDirection/norm(activeDirection),P,1);

validationIndex = find(availableCounts>=validationResolution,1);
if isempty(validationIndex)
    error('No spherical design supplies the requested validation nodes.');
end
[validationDirections,~] = get_sphdesign(availableCounts(validationIndex));
azimuth = 0.37;
elevation = 0.23;
rotationZ = [cos(azimuth) -sin(azimuth) 0; ...
             sin(azimuth) cos(azimuth) 0;0 0 1];
rotationY = [cos(elevation) 0 sin(elevation); ...
             0 1 0;-sin(elevation) 0 cos(elevation)];
validationDirections = (rotationZ*rotationY*validationDirections.').';
validationCount = size(validationDirections,1);
validationPoints = cell(1,P);
validationNormals = cell(1,P);
for body = 1:P
    validationPoints{body} = centres(body,:)+radii(body)*validationDirections;
    validationNormals{body} = validationDirections;
end

xLimits = [min(centreX-radii),max(centreX+radii)] + ...
    plotPaddingInRadii*max(radii)*[-1 1];
zLimits = max(radii)*(1+plotPaddingInRadii)*[-1 1];
xGrid = linspace(xLimits(1),xLimits(2),plotGridCountX);
zGrid = linspace(zLimits(1),zLimits(2),plotGridCountZ);
[gridX,gridZ] = meshgrid(xGrid,zGrid);
gridPoints = [gridX(:),gridZ(:)];
interiorMask = false(size(gridPoints,1),1);
bodyMasks = cell(1,P);
globalSeeds = zeros(P*globalSeedsPerDrop,2);
interiorSeeds = zeros(P*interiorSeedsPerDrop,2);
for body = 1:P
    centreXZ = centres(body,[1 3]);
    bodyMasks{body} = sum((gridPoints-centreXZ).^2,2)<radii(body)^2;
    interiorMask = interiorMask | bodyMasks{body};
    globalAngles = 2*pi*(0:globalSeedsPerDrop-1)'/globalSeedsPerDrop + ...
        pi/(2*globalSeedsPerDrop);
    globalRows = (body-1)*globalSeedsPerDrop+(1:globalSeedsPerDrop);
    globalSeeds(globalRows,:) = centreXZ + ...
        globalSeedRadiusRatio*radii(body)* ...
        [cos(globalAngles),sin(globalAngles)];
    interiorAngles = 2*pi*(0:interiorSeedsPerDrop-1)'/interiorSeedsPerDrop + ...
        pi/(2*interiorSeedsPerDrop);
    interiorRows = (body-1)*interiorSeedsPerDrop+(1:interiorSeedsPerDrop);
    interiorSeeds(interiorRows,:) = centreXZ + ...
        interiorSeedRadiusRatio*radii(body)* ...
        [cos(interiorAngles),sin(interiorAngles)];
end
globalClearance = @(xz) min([xz(1)-xLimits(1),xLimits(2)-xz(1), ...
    xz(2)-zLimits(1),zLimits(2)-xz(2)]);
odeTolerances = [odeRelativeTolerance,odeAbsoluteTolerance];

%% Solve each resolution, retaining only fields, paths, and diagnostics
runCount = numel(requestedResolutions);
actualN = zeros(runCount,1);
actualM = zeros(runCount,1);
iterationCounts = zeros(runCount,1);
solverResiduals = zeros(runCount,1);
relativeBoundaryMaxima = zeros(runCount,3);
translationsByResolution = cell(runCount,1);
globalFields = cell(runCount,1);
interiorRelativeFields = cell(runCount,1);
globalPaths = cell(runCount,1);
interiorPaths = cell(runCount,1);
bodyScales = [abs(S1).*radii,abs(S1).*radii,abs(S1)];
fprintf('P=%d; gap %.3g; random seed %d; independent validation %d nodes/drop.\n', ...
    P,gap,randomSeed,validationCount);
for resolution = 1:runCount
    requestedN = requestedResolutions(resolution);
    [sourceDirections,~] = get_sphdesign(requestedN);
    N = size(sourceDirections,1);
    minimumM = ceil(collocationToSourceRatio*N);
    collocationIndex = find(availableCounts>=minimumM,1);
    if isempty(collocationIndex)
        error('No spherical design supplies M>=ceil(1.2*N).');
    end
    [collocationDirections,~] = ...
        get_sphdesign(availableCounts(collocationIndex));
    M = size(collocationDirections,1);
    actualN(resolution) = N;
    actualM(resolution) = M;

    % Upper estimate for the retained self matrix, compact left and right
    % singular vectors if the block has full column rank. SVD workspace,
    % kernel blocks, and MATLAB copies require additional memory.
    rowCount = 6*M;
    columnCount = 6*N;
    factorStorageGiB = P*8*(2*rowCount*columnCount + ...
        columnCount^2)/2^30;
    fprintf('\nRequested N=%d, actual N=%d; required M>=%d, actual M=%d (M/N=%.3f).\n', ...
        requestedN,N,minimumM,M,M/N);
    fprintf('Approximate retained one-body matrices and SVD vectors: %.2f GiB for %d drops; SVD workspace is additional.\n', ...
        factorStorageGiB,P);

    innerSources = cell(1,P);
    outerSources = cell(1,P);
    surfacePoints = cell(1,P);
    normals = cell(1,P);
    for body = 1:P
        innerSources{body} = centres(body,:)+ ...
            innerProxyRadiusRatio*radii(body)*sourceDirections;
        outerSources{body} = centres(body,:)+ ...
            outerProxyRadiusRatio*radii(body)*sourceDirections;
        surfacePoints{body} = centres(body,:)+ ...
            radii(body)*collocationDirections;
        normals{body} = collocationDirections;
    end

    % Coefficients within each drop are [alpha(3N);gamma(3N)]. The solver
    % uses one local SVD per distinct radius/viscosity and unrestarted GMRES.
    solverOpt.interiorViscosity = interiorViscosities;
    solverOpt.completionScale = 1./(8*pi*radii);
    solverOpt.velocityScale = max(abs(S1).*radii);
    solverOpt.tractionScale = max(abs(S1));
    solverOpt.tractionRowWeight = tractionRowWeight;
    solverOpt.oneBodySvdTolerance = oneBodySvdTolerance;
    solverOpt.krylovTolerance = krylovTolerance;
    solverOpt.maxIterations = maxIterations;
    solverOpt.interactionBatchSize = interactionBatchSize;
    solverOpt.scale_columns = scale_columns;
    [translations,~,gamma,lambda,diagnostics] = ...
        solve_drops_preconditioned(innerSources,outerSources, ...
        surfacePoints,normals,S1,activeDirections,solverOpt);
    translationsByResolution{resolution} = translations;
    iterationCounts(resolution) = diagnostics.iterations;
    solverResiduals(resolution) = diagnostics.regularizedResidual;
    fprintf('GMRES %d iterations; ||H eta-b||/||b|| %.3e; ||A z-b||/||b|| %.3e.\n', ...
        diagnostics.iterations,diagnostics.regularizedResidual, ...
        diagnostics.originalResidual);
    fprintf('Retained one-body ranks: %s; setup %.2fs, solve %.2fs.\n', ...
        mat2str(diagnostics.retainedRanks.'), ...
        diagnostics.setupTime,diagnostics.solveTime);

    maximumPhysicalErrors = drop_validate_surface_conditions( ...
        lambda,gamma,translations,innerSources,outerSources, ...
        validationPoints,validationNormals,activeDirections,S1, ...
        interiorViscosities,validationBatchSize);
    relativeBoundaryMaxima(resolution,:) = ...
        max(maximumPhysicalErrors./bodyScales,[],1);
    fprintf('Independent-grid relative maxima: continuity %.3e; normal velocity %.3e; tangential traction %.3e.\n', ...
        relativeBoundaryMaxima(resolution,:));

    % Every global-grid point selects the proper exterior or interior
    % representation. Interior fields below subtract the matching U_i.
    globalVelocity = @(xz) drop_plane_velocity(xz,0,centres,radii, ...
        innerSources,outerSources,lambda,gamma, ...
        interiorViscosities,interactionBatchSize);
    globalFields{resolution} = globalVelocity(gridPoints);
    interiorRelativeFields{resolution} = nan(size(gridPoints));
    interiorPaths{resolution} = cell(P*interiorSeedsPerDrop,1);
    for body = 1:P
        centreXZ = centres(body,[1 3]);
        interiorVelocity = @(xz) drop_plane_velocity(xz,body, ...
            centres,radii,innerSources,outerSources,lambda,gamma, ...
            interiorViscosities,interactionBatchSize)- ...
            translations(body,[1 3]);
        interiorRelativeFields{resolution}(bodyMasks{body},:) = ...
            interiorVelocity(gridPoints(bodyMasks{body},:));
        interiorClearance = @(xz) radii(body)-norm(xz-centreXZ);
        rows = (body-1)*interiorSeedsPerDrop+(1:interiorSeedsPerDrop);
        interiorPaths{resolution}(rows) = get_streamline( ...
            interiorSeeds(rows,:),interiorVelocity,interiorClearance, ...
            maxArcLength,maxStep,minimumSpeed,odeTolerances);
    end
    globalPaths{resolution} = get_streamline( ...
        globalSeeds,globalVelocity,globalClearance, ...
        maxArcLength,maxStep,minimumSpeed,odeTolerances);
end

%% Field and streamline differences at fixed physical sample locations
globalDifference = vecnorm(globalFields{2}-globalFields{1},2,2);
interiorDifference = vecnorm( ...
    interiorRelativeFields{2}(interiorMask,:)- ...
    interiorRelativeFields{1}(interiorMask,:),2,2);
globalVelocityScale = max(vecnorm(globalFields{2},2,2));
interiorVelocityScale = max(vecnorm( ...
    interiorRelativeFields{2}(interiorMask,:),2,2));
[globalMaximumPathDifference,globalRmsPathDifference,globalLengthDifference] = ...
    compare_streamline_paths(globalPaths{1},globalPaths{2}, ...
    globalSeeds,streamlineComparisonSamples);
[interiorMaximumPathDifference,interiorRmsPathDifference,interiorLengthDifference] = ...
    compare_streamline_paths(interiorPaths{1},interiorPaths{2}, ...
    interiorSeeds,streamlineComparisonSamples);
translationDifference = vecnorm( ...
    translationsByResolution{2}-translationsByResolution{1},2,2);
translationScale = max(vecnorm(translationsByResolution{2},2,2));

fprintf('\nComparison: actual N=%d versus N=%d; same geometry, S_i, seeds, and sample grids.\n', ...
    actualN(1),actualN(2));
fprintf('Global field: max |du| %.3e, RMS %.3e; relative to max fine-grid |u|=%.3e: %.3e, %.3e.\n', ...
    max(globalDifference),sqrt(mean(globalDifference.^2)), ...
    globalVelocityScale,max(globalDifference)/max(globalVelocityScale,realmin), ...
    sqrt(mean(globalDifference.^2))/max(globalVelocityScale,realmin));
fprintf('Interior moving-frame field: max |du_rel| %.3e, RMS %.3e; relative to max fine-grid |u_rel|=%.3e: %.3e, %.3e.\n', ...
    max(interiorDifference),sqrt(mean(interiorDifference.^2)), ...
    interiorVelocityScale, ...
    max(interiorDifference)/max(interiorVelocityScale,realmin), ...
    sqrt(mean(interiorDifference.^2))/max(interiorVelocityScale,realmin));
fprintf('Translations: largest |U_fine-U_coarse| %.3e; relative to max fine-grid |U|=%.3e: %.3e.\n', ...
    max(translationDifference),translationScale, ...
    max(translationDifference)/max(translationScale,realmin));
fprintf('Global streamlines: max matched separation %.3e (%.3e / min radius), max per-line RMS %.3e, max arc-length mismatch %.3e.\n', ...
    max(globalMaximumPathDifference), ...
    max(globalMaximumPathDifference)/min(radii), ...
    max(globalRmsPathDifference),max(globalLengthDifference));
fprintf('Interior streamlines: max matched separation %.3e (%.3e / min radius), max per-line RMS %.3e, max arc-length mismatch %.3e.\n', ...
    max(interiorMaximumPathDifference), ...
    max(interiorMaximumPathDifference)/min(radii), ...
    max(interiorRmsPathDifference),max(interiorLengthDifference));
fprintf('Path distances compare equal arc length from each seed only over the portion traced at both resolutions.\n');
fprintf('Single-point paths (coarse/fine): global %d/%d, interior %d/%d; their traced-length differences remain reported.\n', ...
    nnz(cellfun(@(path) size(path,1)==1,globalPaths{1})), ...
    nnz(cellfun(@(path) size(path,1)==1,globalPaths{2})), ...
    nnz(cellfun(@(path) size(path,1)==1,interiorPaths{1})), ...
    nnz(cellfun(@(path) size(path,1)==1,interiorPaths{2})));
globalMaximumByDrop = zeros(P,1);
interiorMaximumByDrop = zeros(P,1);
globalLengthByDrop = zeros(P,1);
interiorLengthByDrop = zeros(P,1);
for body = 1:P
    globalRows = (body-1)*globalSeedsPerDrop+(1:globalSeedsPerDrop);
    interiorRows = (body-1)*interiorSeedsPerDrop+(1:interiorSeedsPerDrop);
    globalMaximumByDrop(body) = max(globalMaximumPathDifference(globalRows));
    interiorMaximumByDrop(body) = max( ...
        interiorMaximumPathDifference(interiorRows));
    globalLengthByDrop(body) = max(globalLengthDifference(globalRows));
    interiorLengthByDrop(body) = max(interiorLengthDifference(interiorRows));
end
fprintf('Per-drop maxima across streamlines seeded near or inside each drop:\n');
disp(table((1:P).',globalMaximumByDrop,interiorMaximumByDrop, ...
    globalLengthByDrop,interiorLengthByDrop, ...
    'VariableNames',{'Drop','GlobalPathSeparation', ...
    'InteriorPathSeparation','GlobalLengthDifference', ...
    'InteriorLengthDifference'}));
disp(table(requestedResolutions.',actualN,actualM,iterationCounts, ...
    solverResiduals,relativeBoundaryMaxima(:,1), ...
    relativeBoundaryMaxima(:,2),relativeBoundaryMaxima(:,3), ...
    'VariableNames',{'RequestedN','ActualN','ActualM','GMRESIterations', ...
    'GMRESResidual','ContinuityRel','NormalRel','TractionRel'}));

%% Overlay matched streamlines and show pointwise field differences
if view_solution
    theta = linspace(0,2*pi,201);
    figure('Name','Global-frame streamline convergence','Color','w');
    imagesc(xGrid,zGrid,reshape(vecnorm(globalFields{2},2,2),size(gridX)));
    set(gca,'YDir','normal'); hold on;
    for body = 1:P
        plot(centreX(body)+radii(body)*cos(theta), ...
            radii(body)*sin(theta),'w-','LineWidth',1.5, ...
            'HandleVisibility','off');
    end
    coarseHandle = plot(nan,nan,'--','Color',[0.85 0.15 0.15], ...
        'LineWidth',1.2);
    fineHandle = plot(nan,nan,'-','Color',[0.05 0.05 0.05], ...
        'LineWidth',1);
    for line = 1:numel(globalPaths{1})
        coarsePath = globalPaths{1}{line};
        finePath = globalPaths{2}{line};
        plot(coarsePath(:,1),coarsePath(:,2),'--', ...
            'Color',[0.85 0.15 0.15],'LineWidth',1.2, ...
            'HandleVisibility','off');
        plot(finePath(:,1),finePath(:,2),'-', ...
            'Color',[0.05 0.05 0.05],'LineWidth',1, ...
            'HandleVisibility','off');
    end
    axis equal tight; xlim(xLimits); ylim(zLimits);
    colormap(turbo); colorbar;
    xlabel('$x$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    title('Instantaneous lab-frame streamlines; colour is fine-grid $|\mathbf{u}|$', ...
        'Interpreter','latex');
    legend([coarseHandle,fineHandle], ...
        {sprintf('Actual $N=%d$',actualN(1)), ...
         sprintf('Actual $N=%d$',actualN(2))}, ...
        'Interpreter','latex','Location','best');

    figure('Name','Interior moving-frame streamline convergence','Color','w');
    fineInteriorSpeed = vecnorm(interiorRelativeFields{2},2,2);
    speedImage = imagesc(xGrid,zGrid, ...
        reshape(fineInteriorSpeed,size(gridX)));
    set(speedImage,'AlphaData',reshape(interiorMask,size(gridX)));
    set(gca,'YDir','normal'); hold on;
    for body = 1:P
        plot(centreX(body)+radii(body)*cos(theta), ...
            radii(body)*sin(theta),'-','Color',[0.15 0.15 0.15], ...
            'LineWidth',1.3,'HandleVisibility','off');
    end
    coarseHandle = plot(nan,nan,'--','Color',[0.85 0.15 0.15], ...
        'LineWidth',1.2);
    fineHandle = plot(nan,nan,'-','Color',[0.05 0.05 0.05], ...
        'LineWidth',1);
    for line = 1:numel(interiorPaths{1})
        coarsePath = interiorPaths{1}{line};
        finePath = interiorPaths{2}{line};
        plot(coarsePath(:,1),coarsePath(:,2),'--', ...
            'Color',[0.85 0.15 0.15],'LineWidth',1.2, ...
            'HandleVisibility','off');
        plot(finePath(:,1),finePath(:,2),'-', ...
            'Color',[0.05 0.05 0.05],'LineWidth',1, ...
            'HandleVisibility','off');
    end
    axis equal tight; xlim(xLimits); ylim(zLimits);
    colormap(turbo); colorbar;
    xlabel('$x$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    title('Drop-relative interior streamlines; colour is fine-grid $|\mathbf{u}_i-\mathbf{U}_i|$', ...
        'Interpreter','latex');
    legend([coarseHandle,fineHandle], ...
        {sprintf('Actual $N=%d$',actualN(1)), ...
         sprintf('Actual $N=%d$',actualN(2))}, ...
        'Interpreter','latex','Location','best');

    figure('Name','Velocity-field convergence on the common x-z grid', ...
        'Color','w');
    tiledlayout(1,2);
    nexttile;
    imagesc(xGrid,zGrid,reshape(globalDifference,size(gridX)));
    set(gca,'YDir','normal'); hold on;
    for body = 1:P
        plot(centreX(body)+radii(body)*cos(theta), ...
            radii(body)*sin(theta),'w-','LineWidth',1.2);
    end
    axis equal tight; xlim(xLimits); ylim(zLimits);
    colorbar; colormap(turbo);
    xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
    title('$|\mathbf{u}_{N_2}-\mathbf{u}_{N_1}|$ in global frame', ...
        'Interpreter','latex');
    nexttile;
    interiorDifferenceMap = nan(size(gridPoints,1),1);
    interiorDifferenceMap(interiorMask) = interiorDifference;
    differenceImage = imagesc(xGrid,zGrid, ...
        reshape(interiorDifferenceMap,size(gridX)));
    set(differenceImage,'AlphaData',reshape(interiorMask,size(gridX)));
    set(gca,'YDir','normal'); hold on;
    for body = 1:P
        plot(centreX(body)+radii(body)*cos(theta), ...
            radii(body)*sin(theta),'k-','LineWidth',1.2);
    end
    axis equal tight; xlim(xLimits); ylim(zLimits);
    colorbar; colormap(turbo);
    xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
    title('$|(\mathbf{u}_i-\mathbf{U}_i)_{N_2}-(\mathbf{u}_i-\mathbf{U}_i)_{N_1}|$', ...
        'Interpreter','latex');
end
