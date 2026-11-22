function results = demo_random_wall_torque(P,layerThickness)
%DEMO_RANDOM_WALL_TORQUE Mobility of randomly packed torque-driven spheres above a wall.
%
%   results = DEMO_RANDOM_WALL_TORQUE(P, layerThickness)
%
%   Places P equal spheres (radius 1) at random above a no-slip wall at z = 0,
%   applies a unit torque about the z-axis to every sphere, and solves the
%   wall-bounded mobility problem (Blake fundamental solution) with GMRES.
%   The geometry is plotted before the solve; afterwards the surface speed,
%   the independent surface residual and the traction are plotted.
%
%   Usage (StokesHalfspaceFMM must be set up first):
%       setpath; blake.setup;
%       demo_random_wall_torque            % P = 40, center-height span 5
%       demo_random_wall_torque(5,2)       % small example
%       demo_random_wall_torque(40,[])     % 3D cloud instead of a layer
%
%   INPUTS:
%       P              : number of spheres (default 40).
%       layerThickness : exact span of the sphere-center heights for P >= 2
%                        (default 5); the full sphere layer extends one
%                        diameter farther. Pass [] for a 3D cloud.
%
%   OUTPUT:
%       results : struct with fields
%           P, centres (P x 3), geometry (gaps, spans, seed, mode),
%           appliedTorque (6P x 1), rigidVelocity (6P x 1, translation and
%           angular velocity per sphere), iterations (GMRES), details
%           (solver details incl. validation), surfaceTraction (cell, |sigma n|
%           per sphere), geometryFigure, velocityFigure.
%
%   Steps: 1 geometry, 2 geometry preview, 3 solver setup, 4 mobility solve,
%   5 post-processing (surface speed, traction), 6 result plots, 7 summary.

if nargin < 1
    P = 40;
end
if nargin < 2
    layerThickness = 5;
end

%% 1. Geometry: random spheres above the wall
radius = 1;            % sphere radius
mu = 1;                % viscosity
delta = 1;             % minimum gap between spheres
wallGap = 0.5;         % gap between the lowest sphere and the wall
seed = 20261027;       % fixed seed: the same configuration every run

% grow_cluster places each new sphere at the requested gap from the cluster.
rngState = rng;
rngCleanup = onCleanup(@() rng(rngState));
rng(seed,'twister');
if isempty(layerThickness)
    q = radius*grow_cluster(P,delta/radius,3);
    q(:,3) = q(:,3)-min(q(:,3))+radius+wallGap;
    geometryMode = "3D cloud";
else
    validateattributes(layerThickness,{'double'},{'scalar','nonnegative','finite'});
    lateral = radius*grow_cluster(P,delta/radius,2);
    heights = radius+wallGap+layerThickness*rand(P,1);
    heights(1) = radius+wallGap;                     % lowest sphere at wallGap
    if P > 1
        heights(2) = radius+wallGap+layerThickness;  % highest sphere sets the span
    end
    q = [real(lateral),imag(lateral),heights];
    geometryMode = sprintf('center-height span %.3g',layerThickness);
end
clear rngCleanup;

% Gap from each sphere to its nearest neighbour.
nearestBodyGaps = inf(P,1);
for p = 1:P
    others = [1:p-1,p+1:P];
    if ~isempty(others)
        nearestBodyGaps(p) = min(vecnorm(q(others,:)-q(p,:),2,2)-2*radius);
    end
end
if P == 1
    minBodyGap = NaN;
else
    minBodyGap = min(nearestBodyGaps);
end
geometry = struct('P',P,'radius',radius,'delta_body',delta, ...
    'delta_wall',wallGap,'seed',seed,'mode',geometryMode, ...
    'nearestBodyGaps',nearestBodyGaps,'minBodyGap',minBodyGap, ...
    'minWallGap',min(q(:,3)-radius), ...
    'centerHeightSpan',max(q(:,3))-min(q(:,3)), ...
    'sphereLayerExtent',max(q(:,3))+radius-(min(q(:,3))-radius));

fprintf('\nRandom wall cluster: P=%d, %s, minimum body gap %.3g, minimum wall gap %.3g.\n', ...
    P,geometryMode,geometry.minBodyGap,geometry.minWallGap);

%% 2. Geometry preview (before constructing the solver)
% Wall patch just covering the cluster; it also sets the horizontal plot
% limits, and zTop the vertical one, so all plots are tightly framed.
margin = 0.25*radius;
limits = [min(q(:,1:2),[],1)-radius-margin;max(q(:,1:2),[],1)+radius+margin];
wallX = [limits(1,1),limits(2,1);limits(1,1),limits(2,1)];
wallY = [limits(1,2),limits(1,2);limits(2,2),limits(2,2)];
zTop = max(q(:,3))+radius+margin;

geometryFigure = figure('Name',sprintf('Random wall cluster, P=%d',P));
geometryAxes = axes(geometryFigure);
setup_wall_axes(geometryAxes,wallX,wallY,zTop);
draw_spheres(geometryAxes,q,radius,32,[]);
set(geometryAxes,'Color',[0.96,0.97,0.98]);
title(geometryAxes,sprintf('Before solve: P=%d, %s, body gap >= %.2g, wall gap = %.2g', ...
    P,geometryMode,delta,wallGap));
drawnow;

%% 3. Solver setup
opt = struct( ...
    'radius',radius, ...
    'mu',mu, ...
    'wall_height',0, ...                    % wall at z = 0
    'backend',"fmm", ...                    % half-space FMM for the Blake sums
    'preconditioner',"freeSpace", ...       % set to "wall" for wall-aware one-body preconditioners
    'N',400, ...                            % proxy sources per sphere
    'proxyRadius',0.63, ...                 % proxy-sphere radius / sphere radius
    'svdCutoff',1e-12, ...
    'gmresTolerance',1e-7, ...
    'fmmTolerance',1e-9, ...
    'maxIterations',120, ...
    'computeResidual',true, ...             % independent check of the surface residual
    'validationResolution',[18,28]);
opt = wall_solver_options(opt);
opt.radius = radius;
[rin,rout] = init_spheres(q,opt.proxyRadius,opt.N,opt.collocationRatio,radius);

% Identical grids share the free-space one-body preconditioner.
fprintf('Preparing wall solver (one-body preconditioners):\n');
opt.preparedProblem = prepare_wall_problem(q,rin,rout,opt);

%% 4. Mobility solve: unit torque about z on every sphere
unitTorque = [0;0;0;0;0;1];               % [force; torque] per sphere
appliedTorque = repmat(unitTorque,P,1);

fprintf('Solving with GMRES: ');           % GMRES prints one dot per iteration
[rigidVelocity,iterations,~,~,lambda,details] = ...
    solve_mobility_wall(q,rin,rout,appliedTorque,opt);
fprintf('GMRES: %d iterations, true residual %.3e.\n', ...
    iterations,details.trueGmresResidual);

%% 5. Post-processing
% Rigid-body surface speed |U + omega x (x - c)| on every sphere.
[X,Y,Z] = sphere(24);
surfaceSpeed = cell(P,1);
relativePoints = radius*[X(:),Y(:),Z(:)];
for p = 1:P
    rows = (p-1)*6+(1:6);
    translation = rigidVelocity(rows(1:3)).';
    angularVelocity = rigidVelocity(rows(4:6)).';
    pointVelocity = repmat(translation,size(relativePoints,1),1) + ...
        cross(repmat(angularVelocity,size(relativePoints,1),1),relativePoints,2);
    surfaceSpeed{p} = reshape(vecnorm(pointVelocity,2,2),size(X));
end

% Full Blake traction on a coarse uniform grid on every sphere.
fprintf('Evaluating surface traction (Blake, %s)...\n',opt.backend);
tractionResolution = 12;
surfaceTraction = evaluate_surface_traction( ...
    rin,lambda,q,radius,mu,opt.backend,tractionResolution);

% The solver's independent check grid supplies relative velocity residuals.
validation = details.validation;

%% 6. Result plots: surface speed, surface residual, traction
figureHandle = figure('Name',sprintf('Mobility results, P=%d',P));
layout = tiledlayout(figureHandle,1,3,'TileSpacing','compact','Padding','compact');

% Surface speed
velocityAxes = nexttile(layout);
setup_wall_axes(velocityAxes,wallX,wallY,zTop);
draw_spheres(velocityAxes,q,radius,24,surfaceSpeed);
colormap(velocityAxes,parula);
colorbarHandle = colorbar(velocityAxes);
ylabel(colorbarHandle,'Surface speed |u|');
title(velocityAxes,sprintf('Rigid-body surface speed after mobility solve (%d GMRES iterations)', ...
    iterations));

% Independent surface residual (log10, at the validation points)
residualAxes = nexttile(layout);
setup_wall_axes(residualAxes,wallX,wallY,zTop);
draw_spheres(residualAxes,q,radius,20,[]);
errors = validation.pointErrors(:);
positiveErrors = errors(errors>0);
if isempty(positiveErrors)
    errorFloor = realmin;
else
    errorFloor = min(positiveErrors);
end
logErrors = log10(max(errors,errorFloor));
residualLimits = [min(logErrors),max(logErrors)];
if residualLimits(1) == residualLimits(2)
    residualLimits = residualLimits(1)+[-0.5,0.5];
end
scatter3(residualAxes,validation.points(:,1),validation.points(:,2), ...
    validation.points(:,3),12,logErrors,'filled');
clim(residualAxes,residualLimits);
colormap(residualAxes,parula);
residualColorbar = colorbar(residualAxes);
ylabel(residualColorbar,'log_{10} relative surface velocity residual');
title(residualAxes,sprintf('Independent surface residual (max %.2e)', ...
    validation.error));

% Traction magnitude
tractionAxes = nexttile(layout);
setup_wall_axes(tractionAxes,wallX,wallY,zTop);
draw_spheres(tractionAxes,q,radius,tractionResolution,surfaceTraction);
colormap(tractionAxes,hot);
tractionColorbar = colorbar(tractionAxes);
ylabel(tractionColorbar,'Blake traction magnitude |sigma n|');
title(tractionAxes,'Surface traction from Blake tensor');

%% 7. Summary
fprintf('\nSummary\n');
fprintf('  Spheres:                       P = %d (%s)\n',P,geometryMode);
fprintf('  Minimum body gap:              %.3g\n',geometry.minBodyGap);
fprintf('  Minimum wall gap:              %.3g\n',geometry.minWallGap);
fprintf('  Center-height span:            %.3g\n',geometry.centerHeightSpan);
fprintf('  Sphere-layer extent:           %.3g\n',geometry.sphereLayerExtent);
fprintf('  GMRES iterations:              %d\n',iterations);
fprintf('  GMRES true residual:           %.3e\n',details.trueGmresResidual);
fprintf('  Max relative surface residual: %.3e\n',validation.error);

results = struct('P',P,'centres',q,'geometry',geometry, ...
    'appliedTorque',appliedTorque,'rigidVelocity',rigidVelocity, ...
    'iterations',iterations,'details',details, ...
    'surfaceTraction', {surfaceTraction}, ...
    'geometryFigure',geometryFigure,'velocityFigure',figureHandle);
end

% Wall plane, labels, view, tight limits and lighting shared by all plots.
% The view is mainly from the side, turned 30 degrees towards a corner and
% raised 18 degrees, so the wall shows as a thin plane below the spheres.
% The limits are fixed to the wall patch and to [0,zTop] instead of using
% axis equal/vis3d, which zoom out; daspect([1 1 1]) keeps spheres round.
function setup_wall_axes(ax,wallX,wallY,zTop)
hold(ax,'on');
surf(ax,wallX,wallY,zeros(2),'FaceAlpha',0.18, ...
    'FaceColor',[0.35,0.42,0.48],'EdgeColor','none');
xlim(ax,[wallX(1,1),wallX(1,2)]);
ylim(ax,[wallY(1,1),wallY(2,1)]);
zlim(ax,[0,zTop]);
daspect(ax,[1 1 1]);
grid(ax,'on');
xlabel(ax,'x');
ylabel(ax,'y');
zlabel(ax,'z');
view(ax,-30,18);
camlight(ax,'headlight');
lighting(ax,'gouraud');
end

% Draw equal-radius spheres, optionally colored by one scalar field per sphere.
function draw_spheres(ax,centres,radius,resolution,surfaceData)
[X,Y,Z] = sphere(resolution);
for p = 1:size(centres,1)
    x = radius*X+centres(p,1);
    y = radius*Y+centres(p,2);
    z = radius*Z+centres(p,3);
    if isempty(surfaceData)
        surf(ax,x,y,z,'FaceColor',[0.38,0.63,0.76], ...
            'EdgeColor','none','FaceAlpha',1,'FaceLighting','gouraud', ...
            'AmbientStrength',0.28,'DiffuseStrength',0.78, ...
            'SpecularStrength',0.32,'SpecularExponent',20);
    else
        surf(ax,x,y,z,surfaceData{p},'EdgeColor','none', ...
            'FaceColor','interp','FaceLighting','gouraud');
    end
end
end

% Evaluate sigma*n from the complete Blake field on every sphere surface.
function tractionMagnitude = evaluate_surface_traction( ...
    sources,lambda,centres,radius,mu,backend,resolution)
[X,Y,Z] = sphere(resolution);
localPoints = radius*[X(:),Y(:),Z(:)].';
pointCount = size(localPoints,2);
targets = zeros(3,size(centres,1)*pointCount);
normals = zeros(size(targets));
for p = 1:size(centres,1)
    columns = (p-1)*pointCount+(1:pointCount);
    targets(:,columns) = localPoints+radius*[centres(p,1);centres(p,2);centres(p,3)];
    normals(:,columns) = localPoints/radius;
end
options = struct('method',char(backend),'mu',mu,'wall_height',0, ...
    'requested',{{'traction'}},'normals',normals);
fields = blake.fields(sources.',reshape(lambda,3,[]),targets,options);
tractionMagnitude = cell(size(centres,1),1);
for p = 1:size(centres,1)
    columns = (p-1)*pointCount+(1:pointCount);
    tractionMagnitude{p} = reshape(vecnorm(fields.traction(:,columns),2,1),size(X));
end
end
