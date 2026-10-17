function result = sep10_single_ellipsoid_laplace_focal_proxy_sweep
%SEP10_SINGLE_ELLIPSOID_LAPLACE_FOCAL_PROXY_SWEEP Sweep the proxy offset.
%
% Solve the exterior Laplace problem u=1 on one triaxial ellipsoid.  The
% source nodes are moved inward along physical-surface normals.  Errors are
% evaluated on one independent fine surface grid for several resolutions.

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'laplace'));

semiAxes = [1.2,1.1,0.7];
sourceResolutions = [18,12; 24,15; 30,19; 40, 28; 50, 35];
testResolution = [56,36];
sourceOffsets = unique([0.04:0.02:0.4,0.125]);
comparisonOffsets = [0.125,0.30];
svdCutoff = 1e-12;

centre = [0,0,0];
rotation = {eye(3)};

% The focal ellipse belongs to the x-y plane.  The proxy ceases to enclose
% it when the normal offset reaches its minimum physical-surface clearance.
focalEllipseSemiAxes = sqrt(semiAxes(1:2).^2-semiAxes(3)^2);
focalClearance = min(semiAxes(1:2)-focalEllipseSemiAxes);

testSurface = ellipsoid_param(semiAxes(1),semiAxes(2),semiAxes(3));
testSurface = setupsurfquad(testSurface,testResolution);
testPoints = testSurface.x.';

numberOfOffsets = numel(sourceOffsets);
numberOfResolutions = size(sourceResolutions,1);
testMaxResidual = zeros(numberOfOffsets,numberOfResolutions);
collocationMaxResidual = zeros(numberOfOffsets,numberOfResolutions);
charges = zeros(numberOfOffsets,numberOfResolutions);
coefficientNorms = zeros(numberOfOffsets,numberOfResolutions);
keptConditionNumbers = zeros(numberOfOffsets,numberOfResolutions);
retainedRanks = zeros(numberOfOffsets,numberOfResolutions);
solveTimes = zeros(numberOfOffsets,numberOfResolutions);
sourceCounts = zeros(1,numberOfResolutions);
collocationCounts = zeros(1,numberOfResolutions);

fprintf('\nSingle-ellipsoid Laplace proxy-offset sweep\n');
fprintf('  semiaxes                 [%.2f, %.2f, %.2f]\n',semiAxes);
fprintf('  focal-ellipse semiaxes   [%.8f, %.8f]\n',focalEllipseSemiAxes);
fprintf('  enclosure threshold      %.8f\n',focalClearance);
fprintf('  boundary data            u = 1\n');
fprintf('  SVD cutoff               %.1e\n\n',svdCutoff);

for resolutionIndex = 1:numberOfResolutions
    resolution = sourceResolutions(resolutionIndex,:);

    for offsetIndex = 1:numberOfOffsets
        sourceOffset = sourceOffsets(offsetIndex);
        [sourcePoints,collocationPoints] = getEllipsoidGrids( ...
            semiAxes,1,1,resolution(1),resolution(2),sourceOffset, ...
            rotation,centre);

        sourceLevel = sum((sourcePoints./semiAxes).^2,2);
        collocationLevel = sum((collocationPoints./semiAxes).^2,2);
        assert(all(sourceLevel < 1));
        assert(max(abs(collocationLevel-1)) < 1e-12);
        assert(size(collocationPoints,1) > size(sourcePoints,1));

        mfsGrid.sources = sourcePoints;
        mfsGrid.collocation = collocationPoints;
        mfsGrid.sourceBody = ones(size(sourcePoints,1),1);
        mfsGrid.collocationBody = ones(size(collocationPoints,1),1);
        mfsGrid.testPoints = testPoints;
        mfsGrid.testBody = ones(size(testPoints,1),1);

        solution = solve_laplace_capacitance_mfs(mfsGrid,1,svdCutoff);
        testMaxResidual(offsetIndex,resolutionIndex) = ...
            solution.testMaxResidual;
        collocationMaxResidual(offsetIndex,resolutionIndex) = ...
            solution.collocationMaxResidual;
        charges(offsetIndex,resolutionIndex) = solution.charges;
        coefficientNorms(offsetIndex,resolutionIndex) = solution.coefficientNorm;
        keptConditionNumbers(offsetIndex,resolutionIndex) = ...
            solution.keptConditionNumber;
        retainedRanks(offsetIndex,resolutionIndex) = solution.retainedRank;
        solveTimes(offsetIndex,resolutionIndex) = ...
            solution.assemblyTime+solution.solveTime+solution.testEvaluationTime;
    end

    sourceCounts(resolutionIndex) = size(sourcePoints,1);
    collocationCounts(resolutionIndex) = size(collocationPoints,1);
    [bestResidual,bestIndex] = min(testMaxResidual(:,resolutionIndex));
    fprintf('resolution [%d,%d], sources/collocation %d/%d\n', ...
        resolution,sourceCounts(resolutionIndex), ...
        collocationCounts(resolutionIndex));
    fprintf('  best offset %.3f, independent max relative residual %.3e\n', ...
        sourceOffsets(bestIndex),bestResidual);
end

% Plot the actual source nodes against the focal ellipse for two offsets.
focalAngle = linspace(0,2*pi,400);
focalEllipse = [focalEllipseSemiAxes(1)*cos(focalAngle); ...
    focalEllipseSemiAxes(2)*sin(focalAngle); zeros(size(focalAngle))];
[surfaceX,surfaceY,surfaceZ] = ellipsoid( ...
    0,0,0,semiAxes(1),semiAxes(2),semiAxes(3),48);

figure('Name','Laplace proxy nodes and focal ellipse');
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for comparisonIndex = 1:2
    sourceOffset = comparisonOffsets(comparisonIndex);
    resolution = sourceResolutions(end,:);
    sourcePoints = getEllipsoidGrids(semiAxes,1,1, ...
        resolution(1),resolution(2),sourceOffset,rotation,centre);

    ax = nexttile(layout);
    bodyPlot = surf(ax,surfaceX,surfaceY,surfaceZ, ...
        'FaceColor',[0.55,0.65,0.85],'FaceAlpha',0.14,'EdgeColor','none');
    hold(ax,'on');
    sourcePlot = scatter3(ax,sourcePoints(:,1),sourcePoints(:,2), ...
        sourcePoints(:,3),12,[0.80,0.15,0.12],'filled');
    focalPlot = plot3(ax,focalEllipse(1,:),focalEllipse(2,:), ...
        focalEllipse(3,:),'k-','LineWidth',2.5);
    axis(ax,'equal');
    xlim(ax,1.08*[-semiAxes(1),semiAxes(1)]);
    ylim(ax,1.08*[-semiAxes(2),semiAxes(2)]);
    zlim(ax,1.08*[-semiAxes(3),semiAxes(3)]);
    view(ax,34,24);
    grid(ax,'on');
    xlabel(ax,'x');
    ylabel(ax,'y');
    zlabel(ax,'z');
    axis(ax,'tight');

    if sourceOffset < focalClearance
        enclosureText = 'encloses focal ellipse';
    else
        enclosureText = 'does not enclose focal ellipse';
    end
    title(ax,sprintf('offset %.3f: %s',sourceOffset,enclosureText));
    if comparisonIndex == 1
        legend(ax,[bodyPlot,sourcePlot,focalPlot], ...
            {'physical surface','proxy sources','focal ellipse'}, ...
            'Location','northeast');
    end
end
sgtitle(layout,'Laplace proxy nodes and the focal ellipse');

% Plot independent residual versus offset and mark loss of enclosure.
figure('Name','Laplace residual versus proxy offset');
ax = axes;
hold(ax,'on');
set(ax,'YScale','log');
colours = lines(numberOfResolutions);
curveHandles = gobjects(numberOfResolutions,1);
curveLabels = strings(numberOfResolutions,1);
for resolutionIndex = 1:numberOfResolutions
    curveHandles(resolutionIndex) = semilogy(ax,sourceOffsets, ...
        testMaxResidual(:,resolutionIndex),'-o','Color', ...
        colours(resolutionIndex,:),'LineWidth',1.5,'MarkerSize',5);
    resolution = sourceResolutions(resolutionIndex,:);
    curveLabels(resolutionIndex) = sprintf('[%d,%d], N=%d', ...
        resolution,sourceCounts(resolutionIndex));
end
xline(ax,focalClearance,'k--','focal ellipse no longer enclosed', ...
    'LineWidth',1.5,'LabelHorizontalAlignment','left');
xlabel(ax,'normal offset from physical surface');
ylabel(ax,'maximum relative residual on independent surface nodes');
title(ax,'Laplace MFS accuracy versus proxy offset');
legend(ax,curveHandles,curveLabels,'Location','best');
grid(ax,'on');
axis(ax,'tight');

% Plot coefficient growth using the same resolutions and threshold marker.
figure('Name','Laplace coefficient norm versus proxy offset');
ax = axes;
hold(ax,'on');
set(ax,'YScale','log');
coefficientHandles = gobjects(numberOfResolutions,1);
for resolutionIndex = 1:numberOfResolutions
    coefficientHandles(resolutionIndex) = semilogy(ax,sourceOffsets, ...
        coefficientNorms(:,resolutionIndex),'-o','Color', ...
        colours(resolutionIndex,:),'LineWidth',1.5,'MarkerSize',5);
end
xline(ax,focalClearance,'k--','focal ellipse no longer enclosed', ...
    'LineWidth',1.5,'LabelHorizontalAlignment','left');
xlabel(ax,'normal offset from physical surface');
ylabel(ax,'coefficient infinity norm');
title(ax,'Laplace MFS coefficient norm versus proxy offset');
legend(ax,coefficientHandles,curveLabels,'Location','best');
grid(ax,'on');
axis(ax,'tight');

result.semiAxes = semiAxes;
result.sourceResolutions = sourceResolutions;
result.testResolution = testResolution;
result.sourceOffsets = sourceOffsets;
result.svdCutoff = svdCutoff;
result.focalEllipseSemiAxes = focalEllipseSemiAxes;
result.focalClearance = focalClearance;
result.enclosesFocalEllipse = sourceOffsets < focalClearance;
result.sourceCounts = sourceCounts;
result.collocationCounts = collocationCounts;
result.testMaxResidual = testMaxResidual;
result.collocationMaxResidual = collocationMaxResidual;
result.charges = charges;
result.coefficientNorms = coefficientNorms;
result.keptConditionNumbers = keptConditionNumbers;
result.retainedRanks = retainedRanks;
result.solveTimes = solveTimes;
end
