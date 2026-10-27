function results = demo_random_wall_torque(P,layerThickness)
%DEMO_RANDOM_WALL_TORQUE Mobility of randomly packed torque-driven spheres.
% Run setpath; blake.setup; demo_random_wall_torque for P=100 and a layer
% center-height span of 5. Pass [] as the second input for a 3D cloud.
% layerThickness is the exact center-height span for P>=2; the full sphere
% layer extends one diameter farther. Example: demo_random_wall_torque(5,2).
if nargin < 1
    P = 40;
end
if nargin < 2
    layerThickness = 5;
end

% Configuration.
radius = 1;
mu = 1;
delta = 1;
wallGap = 0.5;
seed = 20261027;

% grow_cluster validates P and places each new sphere at the requested gap.
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
    heights(1) = radius+wallGap;
    if P > 1
        heights(2) = radius+wallGap+layerThickness;
    end
    q = [real(lateral),imag(lateral),heights];
    geometryMode = sprintf('center-height span %.3g',layerThickness);
end
clear rngCleanup;

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

% Preview the seeded geometry before constructing the solver or solving.
geometryFigure = figure('Name',sprintf('Random wall cluster, P=%d',P));
geometryAxes = axes(geometryFigure);
hold(geometryAxes,'on');
limits = [min(q(:,1:2),[],1)-2*radius;max(q(:,1:2),[],1)+2*radius];
wallX = [limits(1,1),limits(2,1);limits(1,1),limits(2,1)];
wallY = [limits(1,2),limits(1,2);limits(2,2),limits(2,2)];
surf(geometryAxes,wallX,wallY,zeros(2),'FaceAlpha',0.18, ...
    'FaceColor',[0.35,0.42,0.48],'EdgeColor','none');
draw_spheres(geometryAxes,q,radius,32,[]);
axis(geometryAxes,'equal');
axis(geometryAxes,'vis3d');
grid(geometryAxes,'on');
xlabel(geometryAxes,'x');
ylabel(geometryAxes,'y');
zlabel(geometryAxes,'z');
view(geometryAxes,35,25);
set(geometryAxes,'Projection','perspective','Color',[0.96,0.97,0.98]);
camlight(geometryAxes,'headlight');
lighting(geometryAxes,'gouraud');
title(geometryAxes,sprintf('Before solve: P=%d, %s, body gap >= %.2g, wall gap = %.2g', ...
    P,geometryMode,delta,wallGap));
drawnow;

% Solve mobility; identical grids share the free-space one-body preconditioner.
opt = struct('radius',radius,'mu',mu,'wall_height',0,'backend',"fmm", ...
    'preconditioner',"freeSpace",'N',400,'proxyRadius',0.63, ... %compare: preconditioner, "wall"
    'svdCutoff',1e-12,'gmresTolerance',1e-7,'fmmTolerance',1e-9, ...
    'maxIterations',120,'computeResidual',true, ...
    'validationResolution',[18,28]);
[rin,rout,opt] = init_spheres_wall(q,radius,opt);
opt.preparedProblem = prepare_wall_problem(q,rin,rout,opt);

unitTorque = [0;0;0;0;0;1];
appliedTorque = repmat(unitTorque,P,1);
[rigidVelocity,iterations,~,~,lambda,details] = ...
    solve_mobility_wall(q,rin,rout,appliedTorque,opt);

% Evaluate rigid-body surface speed from the solved translations and rotations.
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

% Evaluate full Blake traction at a coarse uniform grid on every sphere.
tractionResolution = 12;
surfaceTraction = evaluate_surface_traction( ...
    rin,lambda,q,radius,mu,opt.backend,tractionResolution);

% The solver's independent check grid supplies relative velocity residuals.
validation = details.validation;
figureHandle = figure('Name',sprintf('Mobility results, P=%d',P));
layout = tiledlayout(figureHandle,1,3,'TileSpacing','compact','Padding','compact');

velocityAxes = nexttile(layout);
hold(velocityAxes,'on');
surf(velocityAxes,wallX,wallY,zeros(2),'FaceAlpha',0.18, ...
    'FaceColor',[0.35,0.42,0.48],'EdgeColor','none');
draw_spheres(velocityAxes,q,radius,24,surfaceSpeed);
axis(velocityAxes,'equal');
axis(velocityAxes,'vis3d');
grid(velocityAxes,'on');
xlabel(velocityAxes,'x');
ylabel(velocityAxes,'y');
zlabel(velocityAxes,'z');
view(velocityAxes,35,25);
camlight(velocityAxes,'headlight');
lighting(velocityAxes,'gouraud');
colormap(velocityAxes,parula);
colorbarHandle = colorbar(velocityAxes);
ylabel(colorbarHandle,'Surface speed |u|');
title(velocityAxes,sprintf('Rigid-body surface speed after mobility solve (%d GMRES iterations)', ...
    iterations));

residualAxes = nexttile(layout);
hold(residualAxes,'on');
surf(residualAxes,wallX,wallY,zeros(2),'FaceAlpha',0.18, ...
    'FaceColor',[0.35,0.42,0.48],'EdgeColor','none');
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
axis(residualAxes,'equal');
axis(residualAxes,'vis3d');
grid(residualAxes,'on');
xlabel(residualAxes,'x');
ylabel(residualAxes,'y');
zlabel(residualAxes,'z');
view(residualAxes,35,25);
camlight(residualAxes,'headlight');
lighting(residualAxes,'gouraud');
colormap(residualAxes,parula);
residualColorbar = colorbar(residualAxes);
ylabel(residualColorbar,'log_{10} relative surface velocity residual');
title(residualAxes,sprintf('Independent surface residual (max %.2e)', ...
    validation.error));

tractionAxes = nexttile(layout);
hold(tractionAxes,'on');
surf(tractionAxes,wallX,wallY,zeros(2),'FaceAlpha',0.18, ...
    'FaceColor',[0.35,0.42,0.48],'EdgeColor','none');
draw_spheres(tractionAxes,q,radius,tractionResolution,surfaceTraction);
axis(tractionAxes,'equal');
axis(tractionAxes,'vis3d');
grid(tractionAxes,'on');
xlabel(tractionAxes,'x');
ylabel(tractionAxes,'y');
zlabel(tractionAxes,'z');
view(tractionAxes,35,25);
camlight(tractionAxes,'headlight');
lighting(tractionAxes,'gouraud');
colormap(tractionAxes,hot);
tractionColorbar = colorbar(tractionAxes);
ylabel(tractionColorbar,'Blake traction magnitude |sigma n|');
title(tractionAxes,'Surface traction from Blake tensor');

fprintf('Random wall cluster: P=%d, minimum body gap %.3g, minimum wall gap %.3g.\n', ...
    P,geometry.minBodyGap,geometry.minWallGap);
fprintf('Center-height span %.3g, sphere-layer extent %.3g.\n', ...
    geometry.centerHeightSpan,geometry.sphereLayerExtent);
fprintf('Maximum relative surface residual: %.3e.\n',validation.error);
fprintf('Mobility solve: %d GMRES iterations, true residual %.3e.\n', ...
    iterations,details.trueGmresResidual);
results = struct('P',P,'centres',q,'geometry',geometry, ...
    'appliedTorque',appliedTorque,'rigidVelocity',rigidVelocity, ...
    'iterations',iterations,'details',details, ...
    'surfaceTraction', {surfaceTraction}, ...
    'geometryFigure',geometryFigure,'velocityFigure',figureHandle);
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