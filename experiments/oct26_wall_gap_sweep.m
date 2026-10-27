%OKT26_WALL_GAP_SWEEP Solve resistance and mobility as the wall gap varies.

% Fixed experiment parameters.
P = 4;
radius = 1;
mu = 1;
seed = 20261026;
bodyGap = radius;
wallGaps = radius*logspace(-1,1,11);
N = 700;
proxyRadius = 0.63;
svdCutoff = 1e-12;
gmresTolerance = 1e-9;
fmmTolerance = 1e-12;
maxIterations = 400;
wallHeight = 0;

% Keep the same particle arrangement and imposed data at every wall gap.
q = grow_cluster_wall(P,radius,bodyGap,wallGaps(1),seed,"equalHeight");
rngState = rng;
rng(seed+1,'twister');
velocityUnits = repmat([1;1;1;1/radius;1/radius;1/radius],P,1);
forceUnits = repmat([mu*radius*ones(3,1);mu*radius^2*ones(3,1)],P,1);
U = randn(6*P,1);
U = U/(norm(U)/sqrt(P)).*velocityUnits;
FT = randn(6*P,1);
FT = FT/(norm(FT)/sqrt(P)).*forceUnits;
rng(rngState);

opt = struct('radius',radius,'mu',mu,'wall_height',wallHeight, ...
    'backend',"blake",'preconditioner',"freeSpace",'N',N, ... %'backend', "fmm", 'preconditioner',"wall"
    'proxyRadius',proxyRadius,'svdCutoff',svdCutoff, ...
    'gmresTolerance',gmresTolerance,'fmmTolerance',fmmTolerance, ...
    'maxIterations',maxIterations);
gapCount = numel(wallGaps);
resistance = zeros(6*P,gapCount);
mobility = zeros(6*P,gapCount);
resistanceIterations = zeros(1,gapCount);
mobilityIterations = zeros(1,gapCount);
resistanceResidual = zeros(1,gapCount);
mobilityResidual = zeros(1,gapCount);

for j = 1:gapCount
    j
    q(:,3) = wallHeight+radius+wallGaps(j);
    opt.preparedProblem = [];
    [sources,targets,opt] = init_spheres_wall(q,radius,opt);
    opt.preparedProblem = prepare_wall_problem(q,sources,targets,opt);

    [resistance(:,j),resistanceIterations(j),~,resistanceResidual(j)] = ...
        solve_resistance_wall(q,sources,targets,U,opt);
    [mobility(:,j),mobilityIterations(j),~,mobilityResidual(j)] = ...
        solve_mobility_wall(q,sources,targets,FT,opt);

    fprintf('wall gap/a %.3g: resistance %d iterations, mobility %d iterations\n', ...
        wallGaps(j)/radius,resistanceIterations(j),mobilityIterations(j));
end

figure;
layout = tiledlayout(1,2);

ax = nexttile(layout);
semilogx(ax,wallGaps/radius,resistanceIterations,'-o', ...
    'DisplayName','Resistance');
hold(ax,'on');
semilogx(ax,wallGaps/radius,mobilityIterations,'-s', ...
    'DisplayName','Mobility');
xlabel(ax,'Wall gap, \delta/a');
ylabel(ax,'GMRES iterations');
title(ax,'GMRES iterations vs wall gap');
grid(ax,'on');
legend(ax,'Location','best');

ax = nexttile(layout);
loglog(ax,wallGaps/radius,resistanceResidual,'-o', ...
    'DisplayName','Resistance');
hold(ax,'on');
loglog(ax,wallGaps/radius,mobilityResidual,'-s', ...
    'DisplayName','Mobility');
xlabel(ax,'Wall gap, \delta/a');
ylabel(ax,'Residual');
title(ax,'GMRES residuals vs wall gap');
grid(ax,'on');
legend(ax,'Location','best');
