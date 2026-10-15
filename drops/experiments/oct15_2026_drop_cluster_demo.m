% OCT15_2026_DROP_CLUSTER_DEMO  MFS motion of an active drop cluster.
% This is a single solve, in the style of demo/demo_spheres.m. The centres
% come from grow_cluster; each sphere has at least one neighbour at the
% requested surface gap. One-body preconditioned, unrestarted GMRES solves
% the rectangular MFS system without assembling a global dense matrix.
%
% Each body's 6N unknowns are [alpha; gamma], with Cartesian components
% interleaved at each source. Its 6M equations are [velocity continuity;
% normal velocity; two tangential traction components]. The recovered
% exterior Stokeslet strengths are force free, and U_i is the translation.

clear;

%% Parameters
P = 20;
dropRadius = 1;
gap = 1;                         % physical surface-to-surface distance
randomSeed = 15;                 % reproduces both the cluster and S_i
S1Minimum = 0.5;                 % distinct levels are assigned randomly
S1Maximum = 1.5;
activeDirection = [0 0 1];      % common direction of the active traction
interiorViscosity = 2;          % exterior viscosity is one
innerProxyRadius = 0.6;
outerProxyRadius = 2;
requestedN = 400;               % requested source points per drop
collocationToSourceRatio = 1.2;
validationResolution = 1000;    % at least this many independent nodes/drop
validationPointsPerBatch = 128; % targets per evaluation chunk; not a grid size
oneBodySvdTolerance = 1e-13;
krylovTolerance = 1e-10;
maxIterations = 150;            % total unrestarted GMRES iterations
tractionRowWeight = 1;
scale_columns = false;
view_solution = true;           % show velocity and residual figures
maximumArrowLengthInRadii = 2; % display scale, common to all arrows
if view_solution
    close all;
end

%% Repository paths and spherical-design grids
experimentDirectory = fileparts(mfilename('fullpath'));
dropsDirectory = fileparts(experimentDirectory);
repositoryDirectory = fileparts(dropsDirectory);
addpath(dropsDirectory,fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'common','gmres'), ...
    fullfile(repositoryDirectory,'geometry','spheres'));

% get_sphdesign returns the largest available design with at most its
% requested count. For M, select the first available design meeting 1.2N.
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
sourceCount = size(sourceDirections,1);
minimumM = ceil(collocationToSourceRatio*sourceCount);
collocationIndex = find(availableCounts >= minimumM,1);
if isempty(collocationIndex)
    error('No spherical design supplies the requested collocation count.');
end
[collocationDirections,~] = ...
    get_sphdesign(availableCounts(collocationIndex));
collocationCount = size(collocationDirections,1);
validationIndex = find(availableCounts >= ...
    max(validationResolution,collocationCount+1),1);
if isempty(validationIndex)
    error('No spherical design supplies the independent validation count.');
end
[validationDirections,~] = ...
    get_sphdesign(availableCounts(validationIndex));
% Rotate this design so it shares no nodes with the collocation design.
azimuth = 0.37;
elevation = 0.23;
rotationZ = [cos(azimuth) -sin(azimuth) 0; ...
             sin(azimuth) cos(azimuth) 0;0 0 1];
rotationY = [cos(elevation) 0 sin(elevation); ...
             0 1 0;-sin(elevation) 0 cos(elevation)];
validationDirections = (rotationZ*rotationY*validationDirections.').';
validationCount = size(validationDirections,1);

fprintf('Per drop: requested N=%d, actual N=%d; requested M>=%d, actual M=%d (M/N=%.3f).\n', ...
    requestedN,sourceCount,minimumM,collocationCount, ...
    collocationCount/sourceCount);
fprintf('Independent validation: requested >=%d, actual %d points per drop.\n', ...
    validationResolution,validationCount);

%% Cluster and body-specific active traction
rng(randomSeed,'twister');
% grow_cluster assumes unit-radius spheres. Rescale both its gap input and
% its returned centres to give the physical gap above for any dropRadius.
[unitCentres,~] = grow_cluster(P,gap/dropRadius,3);
centres = dropRadius*unitCentres;
minimumActualGap = inf;
for body = 1:P-1
    neighbourGaps = vecnorm(centres(body+1:end,:)-centres(body,:),2,2) ...
        -2*dropRadius;
    minimumActualGap = min(minimumActualGap,min(neighbourGaps));
end
% A random permutation of distinct levels guarantees different S_i values.
S1 = S1Minimum+(S1Maximum-S1Minimum)* ...
    randperm(P).'/(P+1);
activeDirection = activeDirection/norm(activeDirection);
activeDirections = repmat(activeDirection,P,1);

innerSources = cell(1,P);
outerSources = cell(1,P);
surfacePoints = cell(1,P);
normals = cell(1,P);
validationPoints = cell(1,P);
validationNormals = cell(1,P);
for body = 1:P
    innerSources{body} = centres(body,:)+innerProxyRadius*sourceDirections;
    outerSources{body} = centres(body,:)+outerProxyRadius*sourceDirections;
    surfacePoints{body} = centres(body,:)+dropRadius*collocationDirections;
    normals{body} = collocationDirections;
    validationPoints{body} = ...
        centres(body,:)+dropRadius*validationDirections;
    validationNormals{body} = validationDirections;
end

%% One-body setup, unrestarted GMRES, and recovery
% The existing solver accepts one small options struct. The geometry,
% strengths, and resulting velocities remain explicit arrays and cells.
% A shared scale keeps identical translated drops eligible for one SVD.
tractionScale = max(abs(S1));
velocityScale = tractionScale*dropRadius;
completionScale = 1/(8*pi*dropRadius);
solverOpt.interiorViscosity = interiorViscosity;
solverOpt.completionScale = completionScale;
solverOpt.velocityScale = velocityScale;
solverOpt.tractionScale = tractionScale;
solverOpt.tractionRowWeight = tractionRowWeight;
solverOpt.oneBodySvdTolerance = oneBodySvdTolerance;
solverOpt.krylovTolerance = krylovTolerance;
solverOpt.maxIterations = maxIterations;
solverOpt.scale_columns = scale_columns;
% Process all M targets per source body in one direct-interaction batch.
% This controls temporary memory only; it does not change N or M.
solverOpt.interactionBatchSize = collocationCount;

[translations,~,gamma,lambda,diagnostics] = solve_drops_preconditioned( ...
    innerSources,outerSources,surfacePoints,normals,S1, ...
    activeDirections,solverOpt);
% All spheres have the same relative grids and viscosity. Reusing one SVD
% is therefore an expected structural invariant of this experiment.
assert(diagnostics.newFactorCount == 1, ...
    'Congruent drops should share one one-body SVD factorization.');

speeds = vecnorm(translations,2,2);
fprintf('\nP=%d drops; requested gap %.6g, smallest generated gap %.6g; random seed %d.\n', ...
    P,gap,minimumActualGap,randomSeed);
fprintf('Unrestarted GMRES iterations: %d / %d (target %.3e).\n', ...
    diagnostics.iterations,maxIterations,krylovTolerance);
fprintf('Preconditioned residual ||H eta-b||/||b||: %.3e.\n', ...
    diagnostics.regularizedResidual);
fprintf('Original scaled residual ||A z-b||/||b||: %.3e.\n', ...
    diagnostics.originalResidual);
fprintf('One-body rank: %d / %d; setup %.3fs, solve %.3fs.\n', ...
    diagnostics.retainedRanks(1),6*sourceCount, ...
    diagnostics.setupTime,diagnostics.solveTime);
fprintf('Largest coefficient component: %.3e.\n', ...
    diagnostics.maxCoefficientComponent);
disp(table((1:P).',centres(:,1),centres(:,2),centres(:,3), ...
    S1,translations(:,1),translations(:,2),translations(:,3),speeds, ...
    'VariableNames',{'Body','x','y','z','S_i','Ux','Uy','Uz','Speed'}));

%% Surface-condition violations at new, independent nodes
% These are pointwise physical residuals, evaluated off the solve grid.
% Velocity continuity uses |u_out-u_in|. Normal velocity uses
% |n dot (u_out-U_i)|. Traction uses the norm of the two-component
% tangential traction-jump defect, including the prescribed active stress.
[maximumPhysicalErrors,pointPhysicalErrors] = ...
    drop_validate_surface_conditions(lambda,gamma,translations, ...
    innerSources,outerSources,validationPoints,validationNormals, ...
    activeDirections,S1,interiorViscosity,validationPointsPerBatch);

% Exterior viscosity is one. Each body's natural velocity and traction
% scales are |S_i|*dropRadius and |S_i| respectively. S_i is nonzero here.
relativePointErrors = cell(1,P);
maximumRelativeErrors = zeros(P,3);
for body = 1:P
    bodyScales = [abs(S1(body))*dropRadius, ...
        abs(S1(body))*dropRadius,abs(S1(body))];
    relativePointErrors{body} = pointPhysicalErrors{body}./bodyScales;
    maximumRelativeErrors(body,:) = maximumPhysicalErrors(body,:)./bodyScales;
end
fprintf('\nValidation scales per body: velocity |S_i|*radius, traction |S_i|.\n');
fprintf('Maximum pointwise violations at independent nodes:\n');
disp(table((1:P).', ...
    maximumPhysicalErrors(:,1),maximumPhysicalErrors(:,2), ...
    maximumPhysicalErrors(:,3),maximumRelativeErrors(:,1), ...
    maximumRelativeErrors(:,2),maximumRelativeErrors(:,3), ...
    'VariableNames',{'Body','ContinuityAbs','NormalAbs', ...
    'TractionAbs','ContinuityRel','NormalRel','TractionRel'}));
fprintf('Largest relative violation across all bodies: continuity %.3e, normal velocity %.3e, tangential traction %.3e.\n', ...
    max(maximumRelativeErrors,[],1));

%% Translation speed and velocity vectors on one cluster view
if view_solution
    [sphereX,sphereY,sphereZ] = sphere(32);
    velocityFigure = figure('Name','Active drop cluster: translational velocities', ...
        'Color','w','Position',[100 100 1000 700]);
    velocityAxes = axes(velocityFigure);
    hold(velocityAxes,'on');
    for body = 1:P
        surf(velocityAxes,centres(body,1)+dropRadius*sphereX, ...
            centres(body,2)+dropRadius*sphereY, ...
            centres(body,3)+dropRadius*sphereZ, ...
            speeds(body)*ones(size(sphereX)), ...
            'EdgeColor','none','FaceColor','flat','FaceAlpha',0.7, ...
            'FaceLighting','gouraud','AmbientStrength',0.6, ...
            'SpecularStrength',0.12);
    end
    % All arrows have one display scale. Their direction is U_i and their
    % plotted length is proportional to |U_i|. The transparent spheres keep
    % arrow origins at their centres visible.
    arrowScale = maximumArrowLengthInRadii*dropRadius/ ...
        max(max(speeds),realmin);
    velocityArrows = quiver3(velocityAxes, ...
        centres(:,1),centres(:,2),centres(:,3), ...
        arrowScale*translations(:,1),arrowScale*translations(:,2), ...
        arrowScale*translations(:,3),0, ...
        'Color',[0.08 0.08 0.08],'LineWidth',2,'MaxHeadSize',0.7);
    axis(velocityAxes,'equal');
    axis(velocityAxes,'tight');
    axis(velocityAxes,'vis3d');
    grid(velocityAxes,'on');
    view(velocityAxes,35,25);
    set(velocityAxes,'FontSize',11,'LineWidth',0.8, ...
        'GridAlpha',0.15,'Projection','perspective');
    xlabel(velocityAxes,'$x$','Interpreter','latex');
    ylabel(velocityAxes,'$y$','Interpreter','latex');
    zlabel(velocityAxes,'$z$','Interpreter','latex');
    title(velocityAxes,sprintf( ...
        'Translational velocities: $P=%d$, $\\delta=%.2g$',P,gap), ...
        'Interpreter','latex');
    colormap(velocityAxes,turbo);
    speedBar = colorbar(velocityAxes);
    speedBar.Label.String = '$|\mathbf{U}_i|$';
    speedBar.Label.Interpreter = 'latex';
    speedBar.FontSize = 10;
    legend(velocityAxes,velocityArrows, ...
        sprintf('$\\mathbf{U}_i$: arrow length proportional to speed'), ...
        'Interpreter','latex','Location','southoutside','Box','off');
    camlight(velocityAxes,'headlight');

    % Independent-grid residuals, one figure for each surface condition
    conditionNames = {'Velocity continuity |u_{out}-u_{in}|', ...
        'Normal velocity |n \cdot (u_{out}-U_i)|', ...
        'Tangential traction-jump defect'};
    scaleNames = {'|S_i| a','|S_i| a','|S_i|'};
    for condition = 1:3
        figure('Name',sprintf('Drop cluster: surface condition %d',condition));
        hold on;
        for body = 1:P
            scatter3(validationPoints{body}(:,1), ...
                validationPoints{body}(:,2),validationPoints{body}(:,3), ...
                16,relativePointErrors{body}(:,condition),'filled');
        end
        axis equal;
        axis tight;
        grid on;
        view(3);
        xlabel('x'); ylabel('y'); zlabel('z');
        title(sprintf('%s / %s at independent surface nodes', ...
            conditionNames{condition},scaleNames{condition}));
        colormap(turbo);
        errorBar = colorbar;
        errorBar.Label.String = 'Pointwise relative residual';
    end
end
