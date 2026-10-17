function result = sep10_single_ellipsoid_stokes_focal_proxy_experiment
%SEP10_SINGLE_ELLIPSOID_STOKES_FOCAL_PROXY_EXPERIMENT Compare two proxies.
%
% A triaxial ellipsoid translates with unit velocity along its major axis.
% The same Stokes resistance problem is solved with two normal-offset proxy
% surfaces: one encloses the focal set and one does not.  Accuracy is
% measured at a new, finer set of physical-surface points.

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'common','dense'));

% semiAxes = [1.2,1.1,0.7];
semiAxes = [1.1,1.1,0.9];
% sourceResolution = [28,18];
sourceResolution = [35,28];
% testResolution = [47,31];
testResolution = [57,41];
%testResolution = [35,41];
% sourceOffsets = [0.125,0.20];
sourceOffsets = [0.3,0.60];
svdCutoff = 1e-13;

centre = [0,0,0];
rotation = {eye(3)};
rigidMotion = [1;0;0;0;0;0];  % unit translation along the major axis

% The focal ellipse lies in the x-y plane.  Its semiaxes follow from the
% standard confocal family.  For these axes, its minimum clearance from the
% physical surface occurs at an axial endpoint.
focalEllipseSemiAxes = sqrt(semiAxes(1:2).^2-semiAxes(3)^2);
focalClearance = min(semiAxes(1:2)-focalEllipseSemiAxes);
enclosesFocalSet = sourceOffsets < focalClearance;
assert(isequal(enclosesFocalSet,[true,false]));

testSurface = ellipsoid_param(semiAxes(1),semiAxes(2),semiAxes(3));
testSurface = setupsurfquad(testSurface,testResolution);
testPoints = testSurface.x.';
testVelocity = getKmat(testPoints,centre)*rigidMotion;

fprintf('\nSingle-ellipsoid Stokes resistance experiment\n');
fprintf('  semiaxes                 [%.2f, %.2f, %.2f]\n',semiAxes);
fprintf('  focal-ellipse semiaxes   [%.8f, %.8f]\n',focalEllipseSemiAxes);
fprintf('  focal-set tip clearance  %.8f\n',focalClearance);
fprintf('  prescribed motion        unit major-axis translation\n\n');

for proxyCase = 1:2
    sourceOffset = sourceOffsets(proxyCase);

    [sourcePoints,collocationPoints] = getEllipsoidGrids( ...
        semiAxes,1,1,sourceResolution(1),sourceResolution(2), ...
        sourceOffset,rotation,centre);

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
    coefficients = rightVectors(:,retained)*( ...
        (leftVectors(:,retained).'*boundaryVelocity)./singularValues(retained));
    solveTime = toc(solveStart);

    collocationError = stokesMatrix*coefficients-boundaryVelocity;
    testError = stokes_SLP_mat(sourcePoints,testPoints)*coefficients-testVelocity;
    collocationNodeError = vecnorm(reshape(collocationError,3,[]),2,1);
    testNodeError = vecnorm(reshape(testError,3,[]),2,1);
    testSpeed = vecnorm(reshape(testVelocity,3,[]),2,1);

    maxTestResidual = max(testNodeError)/max(testSpeed);
    maxCollocationResidual = max(collocationNodeError)/max(testSpeed);
    forceAndTorque = getKmat(sourcePoints,centre).'*coefficients;

    caseResult.sourceOffset = sourceOffset;
    caseResult.enclosesFocalSet = enclosesFocalSet(proxyCase);
    caseResult.sourcePoints = sourcePoints;
    caseResult.collocationPoints = collocationPoints;
    caseResult.coefficients = coefficients;
    caseResult.forceAndTorque = forceAndTorque;
    caseResult.singularValues = singularValues;
    caseResult.retainedRank = nnz(retained);
    caseResult.keptConditionNumber = singularValues(1)/ ...
        singularValues(find(retained,1,'last'));
    caseResult.coefficientNorm = norm(coefficients);
    caseResult.collocationMaxResidual = maxCollocationResidual;
    caseResult.testMaxResidual = maxTestResidual;
    caseResult.solveTime = solveTime;
    cases(proxyCase) = caseResult; %#ok<AGROW>

    if enclosesFocalSet(proxyCase)
        description = 'encloses focal set';
    else
        description = 'does not enclose focal set';
    end
    fprintf('offset %.3f (%s)\n',sourceOffset,description);
    fprintf('  sources / collocation    %d / %d\n', ...
        size(sourcePoints,1),size(collocationPoints,1));
    fprintf('  retained rank            %d of %d\n', ...
        nnz(retained),numel(singularValues));
    fprintf('  kept condition number    %.3e\n',caseResult.keptConditionNumber);
    fprintf('  coefficient norm         %.3e\n',norm(coefficients));
    fprintf('  collocation max residual %.3e\n',maxCollocationResidual);
    fprintf('  test-node max residual   %.3e\n\n',maxTestResidual);
end

result.semiAxes = semiAxes;
result.sourceResolution = sourceResolution;
result.testResolution = testResolution;
result.svdCutoff = svdCutoff;
result.rigidMotion = rigidMotion;
result.focalEllipseSemiAxes = focalEllipseSemiAxes;
result.focalClearance = focalClearance;
result.testPoints = testPoints;
result.cases = cases;
result.maxResidualRatio = cases(2).testMaxResidual/cases(1).testMaxResidual;
result.relativeForceDifference = norm( ...
    cases(2).forceAndTorque-cases(1).forceAndTorque)/ ...
    norm(cases(1).forceAndTorque);

fprintf('not-enclosed / enclosed new-point residual ratios\n');
fprintf('  maximum                  %.3e\n',result.maxResidualRatio);
fprintf('  relative force change    %.3e\n',result.relativeForceDifference);

% Show the actual source nodes and focal ellipse for both offsets.
focalAngle = linspace(0,2*pi,400);
focalEllipse = [focalEllipseSemiAxes(1)*cos(focalAngle); ...
    focalEllipseSemiAxes(2)*sin(focalAngle); zeros(size(focalAngle))];
[surfaceX,surfaceY,surfaceZ] = ellipsoid( ...
    0,0,0,semiAxes(1),semiAxes(2),semiAxes(3),48);

figure('Name','Focal-ellipse enclosure by normal-offset proxies');
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for proxyCase = 1:2
    ax = nexttile(layout);
    bodyPlot = surf(ax,surfaceX,surfaceY,surfaceZ, ...
        'FaceColor',[0.55,0.65,0.85],'FaceAlpha',0.14,'EdgeColor','none');
    hold(ax,'on');
    sourcePoints = cases(proxyCase).sourcePoints;
    sourcePlot = scatter3(ax,sourcePoints(:,1),sourcePoints(:,2), ...
        sourcePoints(:,3), ...
        12,[0.80,0.15,0.12],'filled');
    focalPlot = plot3(ax,focalEllipse(1,:),focalEllipse(2,:), ...
        focalEllipse(3,:), ...
        'k-','LineWidth',2.5);
    axis(ax,'equal');
    xlim(ax,1.08*[-semiAxes(1),semiAxes(1)]);
    ylim(ax,1.08*[-semiAxes(2),semiAxes(2)]);
    zlim(ax,1.08*[-semiAxes(3),semiAxes(3)]);
    view(ax,34,24);
    grid(ax,'on');
    xlabel(ax,'x');
    ylabel(ax,'y');
    zlabel(ax,'z');

    if cases(proxyCase).enclosesFocalSet
        enclosureText = 'encloses focal ellipse';
    else
        enclosureText = 'does not enclose focal ellipse';
    end
    title(ax,sprintf('offset %.3f: %s', ...
        cases(proxyCase).sourceOffset,enclosureText));
    if proxyCase == 1
        legend(ax,[bodyPlot,sourcePlot,focalPlot], ...
            {'physical surface','proxy sources','focal ellipse'}, ...
            'Location','northeast');
    end
end
title(layout,'Normal-offset proxy nodes and the focal ellipse');
end
