% OCT15_2026_FOUR_DROP_STREAMLINES  Flow around four unequal active drops.
% The drops lie on the x axis, so the x-z section through their centres is
% a symmetry plane for the common z-directed activity. This demo makes two
% views of the solved MFS field: fixed-time lab-frame streamlines that may
% cross drop interfaces, and interior streamlines relative to each U_i.
% A lab-frame streamline is not the path of a fluid particle through the
% moving interface. Run this script to perform the solve and draw the views.

clear;

%% Physical and numerical parameters
P = 4;
radii = [0.8;1.0;1.2;0.9];
interiorViscosities = [1.5;2;3;4]; % exterior viscosity mu_0 = 1
gap = 1;                           % exact gap between neighbouring spheres
innerProxyRadiusRatio = 0.6;      % inner proxy radius is 0.6*a_i
outerProxyRadiusRatio = 2;        % outer proxy radius is 2*a_i
activeDirection = [1 0 0];        % common direction of active surface stress
randomSeed = 15;                  % reproduces the distinct S_i
S1Minimum = 0.5;
S1Maximum = 1.5;
requestedN = 400;                % source nodes on EACH proxy surface
collocationToSourceRatio = 1.2;
validationResolution = 1000;     % minimum independent nodes per surface
oneBodySvdTolerance = 1e-13;
krylovTolerance = 1e-10;
maxIterations = 150;              % unrestarted GMRES iteration limit
tractionRowWeight = 1;
scale_columns = false;
interactionBatchSize = 128;       % targets per temporary kernel block
validationBatchSize = 128;        % targets per validation kernel block
view_solution = true;             % switch both streamline figures on/off

% Streamline parameters. ode45 traces distance along each curve because
% get_streamline integrates u/|u|. These are lengths in the same units as a.
plotPaddingInRadii = 1.5;         % plot box beyond the outermost surfaces
plotGridCountX = 401;
plotGridCountZ = 221;
globalSeedsPerDrop = 8;           % seeds outside each drop
interiorSeedsPerDrop = 8;         % seeds inside each drop
globalSeedRadiusRatio = 1.25;
interiorSeedRadiusRatio = 0.45;
maxArcLength = 20*max(radii);     % in EACH tracing direction
maxStep = 0.04*min(radii);        % largest adaptive ode45 arclength step
minimumSpeed = 1e-6;             % stop nearly stagnant streamlines
odeRelativeTolerance = 1e-6;
odeAbsoluteTolerance = 1e-8*min(radii);

%% Repository routines and available spherical designs
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
% Rotating the validation design separates it from the collocation grid.
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
fprintf('Independent validation: requested >=%d, actual %d per drop.\n', ...
    validationResolution,validationCount);

%% Four centres with exact neighbour gap, then per-drop surfaces and sources
centreX = zeros(P,1);
for body = 2:P
    centreX(body) = centreX(body-1)+radii(body-1)+gap+radii(body);
end
centreX = centreX-mean(centreX);
centres = [centreX zeros(P,2)];
rng(randomSeed,'twister');
S1 = S1Minimum+(S1Maximum-S1Minimum)*randperm(P).'/(P+1);
activeDirection = activeDirection/norm(activeDirection);
activeDirections = repmat(activeDirection,P,1);

innerSources = cell(1,P);
outerSources = cell(1,P);
surfacePoints = cell(1,P);
normals = cell(1,P);
validationPoints = cell(1,P);
validationNormals = cell(1,P);
for body = 1:P
    innerSources{body} = centres(body,:)+ ...
        innerProxyRadiusRatio*radii(body)*sourceDirections;
    outerSources{body} = centres(body,:)+ ...
        outerProxyRadiusRatio*radii(body)*sourceDirections;
    surfacePoints{body} = centres(body,:)+radii(body)*collocationDirections;
    normals{body} = collocationDirections;
    validationPoints{body} = centres(body,:)+radii(body)*validationDirections;
    validationNormals{body} = validationDirections;
end

%% One-body setup, unrestarted GMRES solve, physical coefficients
% The only options struct carries solver controls. Each completion has its
% own beta_i=1/(8*pi*a_i), so U_i=-beta_i*sum_n alpha_i,n.
tractionScale = max(abs(S1));
velocityScale = max(abs(S1).*radii); % mu_0=1
solverOpt.interiorViscosity = interiorViscosities;
solverOpt.completionScale = 1./(8*pi*radii);
solverOpt.velocityScale = velocityScale;
solverOpt.tractionScale = tractionScale;
solverOpt.tractionRowWeight = tractionRowWeight;
solverOpt.oneBodySvdTolerance = oneBodySvdTolerance;
solverOpt.krylovTolerance = krylovTolerance;
solverOpt.maxIterations = maxIterations;
solverOpt.interactionBatchSize = interactionBatchSize;
solverOpt.scale_columns = scale_columns;
[translations,~,gamma,lambda,diagnostics] = solve_drops_preconditioned( ...
    innerSources,outerSources,surfacePoints,normals,S1, ...
    activeDirections,solverOpt);

speeds = vecnorm(translations,2,2);
fprintf('\nP=%d, exact neighbour gap %.3g, S_i random seed %d.\n',P,gap,randomSeed);
fprintf('Unrestarted GMRES: %d iterations, target %.1e, achieved ||H eta-b||/||b||=%.3e.\n', ...
    diagnostics.iterations,krylovTolerance,diagnostics.regularizedResidual);
fprintf('Original scaled residual ||A z-b||/||b||=%.3e; setup %.3fs, solve %.3fs.\n', ...
    diagnostics.originalResidual,diagnostics.setupTime,diagnostics.solveTime);
fprintf('One-body ranks: %s; new factorizations: %d.\n', ...
    mat2str(diagnostics.retainedRanks.'),diagnostics.newFactorCount);
disp(table((1:P).',radii,interiorViscosities,S1, ...
    translations(:,1),translations(:,2),translations(:,3),speeds, ...
    'VariableNames',{'Body','Radius','ViscosityRatio','S_i', ...
    'Ux','Uy','Uz','TranslationSpeed'}));

%% Boundary conditions at a rotated, independent spherical design
% Physical conditions: |u_0-u_i|, |n dot (u_0-U_i)|, and the magnitude of
% the two tangential traction-jump defects. Normalize velocity by
% |S_i|*a_i/mu_0 and traction by |S_i| for each drop.
maximumPhysicalErrors = drop_validate_surface_conditions( ...
    lambda,gamma,translations,innerSources,outerSources, ...
    validationPoints,validationNormals,activeDirections,S1, ...
    interiorViscosities,validationBatchSize);
bodyScales = [abs(S1).*radii,abs(S1).*radii,abs(S1)];
maximumRelativeErrors = maximumPhysicalErrors./bodyScales;
fprintf('\nIndependent-grid maxima, with velocity scale |S_i|*a_i and traction scale |S_i|:\n');
disp(table((1:P).',maximumPhysicalErrors(:,1),maximumPhysicalErrors(:,2), ...
    maximumPhysicalErrors(:,3),maximumRelativeErrors(:,1), ...
    maximumRelativeErrors(:,2),maximumRelativeErrors(:,3), ...
    'VariableNames',{'Body','ContinuityAbs','NormalAbs', ...
    'TractionAbs','ContinuityRel','NormalRel','TractionRel'}));
fprintf('Across all drops, relative maxima: continuity %.3e, normal %.3e, traction %.3e.\n', ...
    max(maximumRelativeErrors,[],1));

%% Two x-z plane views of the solved velocity field
if view_solution
    odeTolerances = [odeRelativeTolerance,odeAbsoluteTolerance];
    xLimits = [min(centreX-radii),max(centreX+radii)] + ...
        plotPaddingInRadii*max(radii)*[-1 1];
    zLimits = max(radii)*(1+plotPaddingInRadii)*[-1 1];
    xGrid = linspace(xLimits(1),xLimits(2),plotGridCountX);
    zGrid = linspace(zLimits(1),zLimits(2),plotGridCountZ);
    [gridX,gridZ] = meshgrid(xGrid,zGrid);
    gridPoints = [gridX(:),gridZ(:)];
    theta = linspace(0,2*pi,201);

    % The lab-frame field is selected by physical region at EVERY target,
    % including ode45's intermediate evaluations near an interface.
    globalVelocity = @(xz) drop_plane_velocity(xz,0,centres,radii, ...
        innerSources,outerSources,lambda,gamma, ...
        interiorViscosities,interactionBatchSize);
    globalClearance = @(xz) min([xz(1)-xLimits(1), ...
        xLimits(2)-xz(1),xz(2)-zLimits(1),zLimits(2)-xz(2)]);
    globalSpeed = vecnorm(globalVelocity(gridPoints),2,2);
    globalSeeds = zeros(P*globalSeedsPerDrop,2);
    for body = 1:P
        seedAngles = 2*pi*(0:globalSeedsPerDrop-1)'/globalSeedsPerDrop + ...
            pi/(2*globalSeedsPerDrop);
        selected = (body-1)*globalSeedsPerDrop+(1:globalSeedsPerDrop);
        globalSeeds(selected,:) = centres(body,[1 3]) + ...
            globalSeedRadiusRatio*radii(body)* ...
            [cos(seedAngles),sin(seedAngles)];
    end
    globalPaths = get_streamline(globalSeeds,globalVelocity,globalClearance, ...
        maxArcLength,maxStep,minimumSpeed,odeTolerances);

    figure('Name','Four drops: instantaneous global-frame streamlines', ...
        'Color','w');
    imagesc(xGrid,zGrid,reshape(globalSpeed,size(gridX)));
    set(gca,'YDir','normal');
    hold on;
    for body = 1:P
        plot(centreX(body)+radii(body)*cos(theta), ...
            radii(body)*sin(theta),'w-','LineWidth',1.5);
    end
    for line = 1:numel(globalPaths)
        path = globalPaths{line};
        plot(path(:,1),path(:,2),'-','Color',[0.05 0.05 0.05], ...
            'LineWidth',0.9);
    end
    axis equal;
    xlim(xLimits);
    ylim(zLimits);
    axis tight;
    colorbar;
    colormap(turbo);
    xlabel('$x$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    title('Instantaneous lab-frame streamlines, speed $|\mathbf{u}|$', ...
        'Interpreter','latex');

    % Within drop i, subtract its own translation. The signed distance
    % a_i-|x-c_i| is an ode45 terminal event: these lines stop at the
    % interface. No exterior field appears in this moving-frame view.
    interiorSpeed = nan(size(gridPoints,1),1);
    interiorPaths = cell(P,1);
    for body = 1:P
        centreXZ = centres(body,[1 3]);
        inside = sum((gridPoints-centreXZ).^2,2)<radii(body)^2;
        interiorVelocity = @(xz) drop_plane_velocity(xz,body, ...
            centres,radii,innerSources,outerSources,lambda,gamma, ...
            interiorViscosities,interactionBatchSize)- ...
            translations(body,[1 3]);
        interiorSpeed(inside) = vecnorm( ...
            interiorVelocity(gridPoints(inside,:)),2,2);
        interiorClearance = @(xz) radii(body)-norm(xz-centreXZ);
        seedAngles = 2*pi*(0:interiorSeedsPerDrop-1)'/ ...
            interiorSeedsPerDrop + pi/(2*interiorSeedsPerDrop);
        seeds = centreXZ+interiorSeedRadiusRatio*radii(body)* ...
            [cos(seedAngles),sin(seedAngles)];
        interiorPaths{body} = get_streamline( ...
            seeds,interiorVelocity,interiorClearance, ...
            maxArcLength,maxStep,minimumSpeed,odeTolerances);
    end

    figure('Name','Four drops: moving-frame interior streamlines', ...
        'Color','w');
    speedImage = imagesc(xGrid,zGrid,reshape(interiorSpeed,size(gridX)));
    set(speedImage,'AlphaData',reshape(isfinite(interiorSpeed),size(gridX)));
    set(gca,'YDir','normal');
    hold on;
    for body = 1:P
        plot(centreX(body)+radii(body)*cos(theta), ...
            radii(body)*sin(theta),'-','Color',[0.15 0.15 0.15], ...
            'LineWidth',1.4);
        for line = 1:numel(interiorPaths{body})
            path = interiorPaths{body}{line};
            plot(path(:,1),path(:,2),'-','Color',[0.08 0.08 0.08], ...
                'LineWidth',0.9);
        end
    end
    axis equal;
    xlim(xLimits);
    ylim(zLimits);
    axis tight;
    colorbar;
    colormap(turbo);
    xlabel('$x$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    title('Interior streamlines in each drop frame, speed $|\mathbf{u}_i-\mathbf{U}_i|$', ...
        'Interpreter','latex');
end
