function result = sep10_single_ellipsoid_stokes_focal_proxy_sweep
%SEP10_SINGLE_ELLIPSOID_STOKES_FOCAL_PROXY_SWEEP Sweep the proxy offset.
%
% Solve one resistance problem for one translating triaxial ellipsoid.
% Residuals are evaluated both at the collocation nodes and on an
% independent surface grid.

% ----------------------- USER SETTING -------------------------------
% Change this 3-vector to prescribe the ellipsoid translation [Ux;Uy;Uz].
translationVelocity = [1;0;0];
% --------------------------------------------------------------------

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'common','dense'));

semiAxes = [1.2,1.1,0.7];
sourceResolutions = [18,12; 24,15; 30,19; 40, 28; 50, 35];
testResolution = [56,36];
sourceOffsets = unique([0.04:0.02:0.4,0.125]);
comparisonOffsets = [0.125,0.30];
svdCutoff = 1e-12;

centre = [0,0,0];
rotation = {eye(3)};
if ~isnumeric(translationVelocity) || numel(translationVelocity) ~= 3 || ...
        any(~isfinite(translationVelocity)) || norm(translationVelocity) == 0
    error('translationVelocity must be a finite, nonzero 3-vector.');
end
translationVelocity = translationVelocity(:);
rigidMotion = [translationVelocity;zeros(3,1)];

% The focal ellipse lies in the x-y plane.  The normal-offset proxy stops
% enclosing it at the minimum clearance between the ellipse and the body.
focalEllipseSemiAxes = sqrt(semiAxes(1:2).^2-semiAxes(3)^2);
focalClearance = min(semiAxes(1:2)-focalEllipseSemiAxes);

testSurface = ellipsoid_param(semiAxes(1),semiAxes(2),semiAxes(3));
testSurface = setupsurfquad(testSurface,testResolution);
testPoints = testSurface.x.';
testVelocity = getKmat(testPoints,centre)*rigidMotion;

numberOfOffsets = numel(sourceOffsets);
numberOfResolutions = size(sourceResolutions,1);

collocationMaxResidual = zeros(numberOfOffsets,numberOfResolutions);
testMaxResidual = zeros(numberOfOffsets,numberOfResolutions);
coefficientNorms = zeros(numberOfOffsets,numberOfResolutions);
forceAndTorque = zeros(6,numberOfOffsets,numberOfResolutions);
keptConditionNumbers = zeros(numberOfOffsets,numberOfResolutions);
retainedRanks = zeros(numberOfOffsets,numberOfResolutions);
solveTimes = zeros(numberOfOffsets,numberOfResolutions);
sourceCounts = zeros(1,numberOfResolutions);
collocationCounts = zeros(1,numberOfResolutions);

fprintf('\nSingle-ellipsoid Stokes proxy-offset sweep\n');
fprintf('  semiaxes                 [%.2f, %.2f, %.2f]\n',semiAxes);
fprintf('  focal-ellipse semiaxes   [%.8f, %.8f]\n',focalEllipseSemiAxes);
fprintf('  enclosure threshold      %.8f\n',focalClearance);
fprintf('  translation velocity     [%g, %g, %g]^T\n',translationVelocity);
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

        boundaryVelocity = getKmat(collocationPoints,centre)*rigidMotion;

        solveStart = tic;
        stokesMatrix = stokes_SLP_mat(sourcePoints,collocationPoints);
        [leftVectors,singularMatrix,rightVectors] = svd(stokesMatrix,'econ');
        singularValues = diag(singularMatrix);
        retained = singularValues >= svdCutoff*singularValues(1);
        projectedData = leftVectors(:,retained).'*boundaryVelocity;
        coefficients = rightVectors(:,retained)*bsxfun( ...
            @rdivide,projectedData,singularValues(retained));
        solveTimes(offsetIndex,resolutionIndex) = toc(solveStart);

        collocationError = stokesMatrix*coefficients-boundaryVelocity;
        testError = stokes_SLP_mat(sourcePoints,testPoints)* ...
            coefficients-testVelocity;
        sourceRigidMap = getKmat(sourcePoints,centre);

        collocationNodeError = vecnorm( ...
            reshape(collocationError,3,[]),2,1);
        testNodeError = vecnorm(reshape(testError,3,[]),2,1);
        collocationSpeed = vecnorm( ...
            reshape(boundaryVelocity,3,[]),2,1);
        testSpeed = vecnorm(reshape(testVelocity,3,[]),2,1);

        collocationMaxResidual(offsetIndex,resolutionIndex) = ...
            max(collocationNodeError)/max(collocationSpeed);
        testMaxResidual(offsetIndex,resolutionIndex) = ...
            max(testNodeError)/max(testSpeed);
        coefficientNorms(offsetIndex,resolutionIndex) = norm(coefficients);
        forceAndTorque(:,offsetIndex,resolutionIndex) = ...
            sourceRigidMap.'*coefficients;

        retainedRanks(offsetIndex,resolutionIndex) = nnz(retained);
        keptConditionNumbers(offsetIndex,resolutionIndex) = ...
            singularValues(1)/singularValues(find(retained,1,'last'));
    end

    sourceCounts(resolutionIndex) = size(sourcePoints,1);
    collocationCounts(resolutionIndex) = size(collocationPoints,1);
    fprintf('resolution [%d,%d], sources/collocation %d/%d\n', ...
        resolution,sourceCounts(resolutionIndex), ...
        collocationCounts(resolutionIndex));
    [bestResidual,bestIndex] = min(testMaxResidual(:,resolutionIndex));
    bestCollocationResidual = collocationMaxResidual( ...
        bestIndex,resolutionIndex);
    fprintf('  best offset %.3f\n',sourceOffsets(bestIndex));
    fprintf('    collocation max relative residual %.3e\n', ...
        bestCollocationResidual);
    fprintf('    test-node max relative residual   %.3e\n',bestResidual);
end

% Show the physical surface, actual source nodes, and focal ellipse.
focalAngle = linspace(0,2*pi,400);
focalEllipse = [focalEllipseSemiAxes(1)*cos(focalAngle); ...
    focalEllipseSemiAxes(2)*sin(focalAngle); zeros(size(focalAngle))];
[surfaceX,surfaceY,surfaceZ] = ellipsoid( ...
    0,0,0,semiAxes(1),semiAxes(2),semiAxes(3),48);

figure('Name','Stokes proxy nodes and focal ellipse');
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
    axis(ax,'tight');
    view(ax,34,24);
    grid(ax,'on');
    xlabel(ax,'x');
    ylabel(ax,'y');
    zlabel(ax,'z');

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
title(layout,'One ellipsoid: normal-offset proxies and focal ellipse');

% Plot the only two reported errors: collocation and independent-test max.
figure('Name','Stokes residual versus proxy offset');
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
colours = lines(numberOfResolutions);
curveLabels = strings(numberOfResolutions,1);
for errorIndex = 1:2
    ax = nexttile(layout);
    hold(ax,'on');
    set(ax,'YScale','log');
    curveHandles = gobjects(numberOfResolutions,1);
    for resolutionIndex = 1:numberOfResolutions
        if errorIndex == 1
            residuals = collocationMaxResidual(:,resolutionIndex);
        else
            residuals = testMaxResidual(:,resolutionIndex);
        end
        curveHandles(resolutionIndex) = semilogy(ax,sourceOffsets, ...
            residuals,'-o','Color',colours(resolutionIndex,:), ...
            'LineWidth',1.4,'MarkerSize',4);
        resolution = sourceResolutions(resolutionIndex,:);
        curveLabels(resolutionIndex) = sprintf('[%d,%d], N=%d', ...
            resolution,sourceCounts(resolutionIndex));
    end
    xline(ax,focalClearance,'k--','focal ellipse not enclosed', ...
        'LineWidth',1.3,'LabelHorizontalAlignment','left');
    xlabel(ax,'normal offset from physical surface');
    if errorIndex == 1
        ylabel(ax,'collocation max relative residual');
    else
        ylabel(ax,'test-node max relative residual');
        legend(ax,curveHandles,curveLabels,'Location','best');
    end
    grid(ax,'on');
    axis(ax,'tight');
end
title(layout,sprintf('Stokes MFS accuracy for U = [%g, %g, %g]^T', ...
    translationVelocity));

% Plot coefficient growth for the prescribed translation.
figure('Name','Stokes coefficient norm versus proxy offset');
ax = axes;
hold(ax,'on');
set(ax,'YScale','log');
coefficientHandles = gobjects(numberOfResolutions,1);
for resolutionIndex = 1:numberOfResolutions
    coefficientHandles(resolutionIndex) = semilogy(ax,sourceOffsets, ...
        coefficientNorms(:,resolutionIndex),'-o', ...
        'Color',colours(resolutionIndex,:), ...
        'LineWidth',1.4,'MarkerSize',4);
end
xline(ax,focalClearance,'k--','focal ellipse not enclosed', ...
    'LineWidth',1.3,'LabelHorizontalAlignment','left');
xlabel(ax,'normal offset from physical surface');
ylabel(ax,'coefficient 2-norm');
title(ax,sprintf('Stokes MFS coefficient norm for U = [%g, %g, %g]^T', ...
    translationVelocity));
legend(ax,coefficientHandles,curveLabels,'Location','best');
grid(ax,'on');
axis(ax,'tight');

result.semiAxes = semiAxes;
result.sourceResolutions = sourceResolutions;
result.testResolution = testResolution;
result.sourceOffsets = sourceOffsets;
result.svdCutoff = svdCutoff;
result.translationVelocity = translationVelocity;
result.rigidMotion = rigidMotion;
result.focalEllipseSemiAxes = focalEllipseSemiAxes;
result.focalClearance = focalClearance;
result.enclosesFocalEllipse = sourceOffsets < focalClearance;
result.sourceCounts = sourceCounts;
result.collocationCounts = collocationCounts;
result.collocationMaxResidual = collocationMaxResidual;
result.testMaxResidual = testMaxResidual;
result.coefficientNorms = coefficientNorms;
result.forceAndTorque = forceAndTorque;
result.keptConditionNumbers = keptConditionNumbers;
result.retainedRanks = retainedRanks;
result.solveTimes = solveTimes;
end
