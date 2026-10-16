function result = sep10_single_ellipsoid_laplace_confocal_multipole_experiment( ...
    show_linear)
%SEP10_SINGLE_ELLIPSOID_LAPLACE_CONFOCAL_MULTIPOLE_EXPERIMENT Compare bases.
%
% result = SEP10_SINGLE_ELLIPSOID_LAPLACE_CONFOCAL_MULTIPOLE_EXPERIMENT
% result = SEP10_SINGLE_ELLIPSOID_LAPLACE_CONFOCAL_MULTIPOLE_EXPERIMENT(show_linear)
%
% Solve u=1 on one triaxial ellipsoid using sources on its filled focal
% ellipse.  The same source nodes are used for two MFS bases:
%
%   1. one Laplace monopole per focal-ellipse node;
%   2. one monopole and three Cartesian dipoles per node.
%
% Both systems use the same fixed collocation nodes on the true ellipsoid
% and are validated on a separate, finer true-surface grid.  The dipole
% matrix is multiplied by dipoleLengthScale, so all fitted coefficients
% have charge units and their Euclidean norms can be compared.  Physical
% dipole moments equal dipoleLengthScale times the fitted dipole vectors.
% By default, source magnitudes and surface residuals are shown only on a
% log10 scale.  Set show_linear to 1 to add their linear-scale panels.

if nargin < 1 || isempty(show_linear)
    show_linear = 0;
end
validateattributes(show_linear,{'numeric','logical'}, ...
    {'real','finite','scalar'},mfilename,'show_linear',1);
if show_linear ~= 0 && show_linear ~= 1
    error('sep10_single_ellipsoid_laplace_confocal_multipole_experiment:BadShowLinear', ...
        'show_linear must be 0 or 1.');
end
show_linear = logical(show_linear);

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory);
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'laplace'));

% User-adjustable experiment parameters.
semiAxes = [1.2,1.1,0.7];
ellipseResolutions = [8,3; 12,4; 16,5; 20,6; 24,7; 28,8; 32 9; 36 10; 40 14];
sourceDistributionType = 4;
collocationResolution = [80,60];
testResolution = [83,55];
svdCutoff = 1e-13;  % Change this value to alter every SVD truncation.
dipoleLengthScale = semiAxes(3);
methodNames = ["monopoles","monopoles + dipoles"];

surface = ellipsoid_param(semiAxes(1),semiAxes(2),semiAxes(3));
collocationSurface = setupsurfquad(surface,collocationResolution);
testSurface = setupsurfquad(surface,testResolution);
collocationPoints = collocationSurface.x.';
testPoints = testSurface.x.';
collocationData = ones(size(collocationPoints,1),1);
testData = ones(size(testPoints,1),1);

numberOfResolutions = size(ellipseResolutions,1);
numberOfMethods = numel(methodNames);
sourceNodeCounts = prod(ellipseResolutions,2);
maximumCombinedUnknowns = 4*max(sourceNodeCounts);
assert(size(collocationPoints,1) > maximumCombinedUnknowns, ...
    'The fixed collocation system must be overdetermined at every resolution.');

% With G=1/(4*pi*r), the exact isolated-ellipsoid capacity is 8*pi/I.
capacityIntegral = integral(@(s) 1./sqrt( ...
    (semiAxes(1)^2+s).*(semiAxes(2)^2+s).*(semiAxes(3)^2+s)), ...
    0,Inf,'RelTol',1e-13,'AbsTol',1e-14);
exactCapacity = 8*pi/capacityIntegral;

collocationMaxResidual = zeros(numberOfResolutions,numberOfMethods);
testMaxResidual = zeros(numberOfResolutions,numberOfMethods);
capacities = zeros(numberOfResolutions,numberOfMethods);
capacityRelativeError = zeros(numberOfResolutions,numberOfMethods);
coefficientNorms = zeros(numberOfResolutions,numberOfMethods);
monopoleCoefficientNorms = zeros(numberOfResolutions,numberOfMethods);
scaledDipoleCoefficientNorms = zeros(numberOfResolutions,numberOfMethods);
retainedRanks = zeros(numberOfResolutions,numberOfMethods);
keptConditionNumbers = zeros(numberOfResolutions,numberOfMethods);
assemblyTimes = zeros(numberOfResolutions,1);
solveTimes = zeros(numberOfResolutions,numberOfMethods);

finestSourcePoints = [];
finestMonopoleOnlyCoefficients = [];
finestCombinedMonopoleCoefficients = [];
finestCombinedScaledDipoleCoefficients = [];
finestTestResiduals = zeros(size(testPoints,1),numberOfMethods);

fprintf('\nSingle-ellipsoid Laplace confocal multipole experiment\n');
fprintf('  semiaxes                    [%.2f, %.2f, %.2f]\n',semiAxes);
fprintf('  focal-node distribution     type %d\n',sourceDistributionType);
fprintf('  collocation/test nodes      %d / %d\n', ...
    size(collocationPoints,1),size(testPoints,1));
fprintf('  dipole column length scale  %.3f\n',dipoleLengthScale);
fprintf('  SVD cutoff                  %.1e\n',svdCutoff);
fprintf('  exact capacity              %.12f\n\n',exactCapacity);

for resolutionIndex = 1:numberOfResolutions
    resolution = ellipseResolutions(resolutionIndex,:);
    sourcePoints = discretize_confocal_ellipse( ...
        semiAxes(1),semiAxes(2),semiAxes(3), ...
        resolution(1),resolution(2),sourceDistributionType);
    numberOfSources = size(sourcePoints,1);

    sourceLevel = sourcePoints(:,1).^2/(semiAxes(1)^2-semiAxes(3)^2) + ...
        sourcePoints(:,2).^2/(semiAxes(2)^2-semiAxes(3)^2);
    assert(all(sourceLevel > 0));
    if sourceDistributionType ~= 4
        assert(all(sourceLevel < 1));
    end
    assert(all(sourcePoints(:,3) == 0));

    assemblyStart = tic;
    collocationMonopoles = lap3d_charge_potential_matrix( ...
        collocationPoints,sourcePoints);
    collocationDipoles = laplace_dipole_potential_matrix( ...
        collocationPoints,sourcePoints);
    testMonopoles = lap3d_charge_potential_matrix(testPoints,sourcePoints);
    testDipoles = laplace_dipole_potential_matrix(testPoints,sourcePoints);
    assemblyTimes(resolutionIndex) = toc(assemblyStart);

    collocationMatrices = {collocationMonopoles, ...
        [collocationMonopoles,dipoleLengthScale*collocationDipoles]};
    testMatrices = {testMonopoles, ...
        [testMonopoles,dipoleLengthScale*testDipoles]};

    fprintf('resolution [%d,%d], focal nodes %d\n', ...
        resolution,numberOfSources);
    for methodIndex = 1:numberOfMethods
        systemMatrix = collocationMatrices{methodIndex};
        assert(size(systemMatrix,1) > size(systemMatrix,2));

        solveStart = tic;
        [coefficients,retainedRank,keptConditionNumber] = ...
            truncated_svd_solve(systemMatrix,collocationData,svdCutoff);
        solveTimes(resolutionIndex,methodIndex) = toc(solveStart);

        monopoleCoefficients = coefficients(1:numberOfSources);
        if methodIndex == 1
            scaledDipoleCoefficients = zeros(numberOfSources,3);
        else
            % The combined ordering is [q; d_x; d_y; d_z].
            scaledDipoleCoefficients = reshape( ...
                coefficients(numberOfSources+1:end),numberOfSources,3);
        end

        collocationError = systemMatrix*coefficients-collocationData;
        testError = testMatrices{methodIndex}*coefficients-testData;
        capacity = sum(monopoleCoefficients);

        collocationMaxResidual(resolutionIndex,methodIndex) = ...
            norm(collocationError,inf);
        testMaxResidual(resolutionIndex,methodIndex) = norm(testError,inf);
        capacities(resolutionIndex,methodIndex) = capacity;
        capacityRelativeError(resolutionIndex,methodIndex) = ...
            abs(capacity-exactCapacity)/exactCapacity;
        coefficientNorms(resolutionIndex,methodIndex) = norm(coefficients);
        monopoleCoefficientNorms(resolutionIndex,methodIndex) = ...
            norm(monopoleCoefficients);
        scaledDipoleCoefficientNorms(resolutionIndex,methodIndex) = ...
            norm(scaledDipoleCoefficients,'fro');
        retainedRanks(resolutionIndex,methodIndex) = retainedRank;
        keptConditionNumbers(resolutionIndex,methodIndex) = ...
            keptConditionNumber;

        fprintf('  %-20s unknowns %4d, rank %4d\n', ...
            methodNames(methodIndex),size(systemMatrix,2),retainedRank);
        fprintf('    collocation/test max residual  %.3e / %.3e\n', ...
            collocationMaxResidual(resolutionIndex,methodIndex), ...
            testMaxResidual(resolutionIndex,methodIndex));
        fprintf('    capacity relative error        %.3e\n', ...
            capacityRelativeError(resolutionIndex,methodIndex));
        fprintf('    scaled coefficient 2-norm      %.3e\n', ...
            coefficientNorms(resolutionIndex,methodIndex));

        if resolutionIndex == numberOfResolutions
            finestSourcePoints = sourcePoints;
            finestTestResiduals(:,methodIndex) = testError;
            if methodIndex == 1
                finestMonopoleOnlyCoefficients = monopoleCoefficients;
            else
                finestCombinedMonopoleCoefficients = monopoleCoefficients;
                finestCombinedScaledDipoleCoefficients = ...
                    scaledDipoleCoefficients;
            end
        end
    end
end

colours = lines(numberOfMethods);
figure('Name','Laplace confocal-basis convergence');
layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
for methodIndex = 1:numberOfMethods
    semilogy(ax,sourceNodeCounts,testMaxResidual(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-ellipse nodes');
ylabel(ax,'independent max residual');
title(ax,'Boundary accuracy');
grid(ax,'on');
axis(ax,'tight');

ax = nexttile(layout);
for methodIndex = 1:numberOfMethods
    semilogy(ax,sourceNodeCounts,capacityRelativeError(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-ellipse nodes');
ylabel(ax,'relative capacity error');
title(ax,'Capacity accuracy');
grid(ax,'on');
axis(ax,'tight');

ax = nexttile(layout);
curveHandles = gobjects(numberOfMethods,1);
for methodIndex = 1:numberOfMethods
    curveHandles(methodIndex) = semilogy(ax,sourceNodeCounts, ...
        coefficientNorms(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-ellipse nodes');
ylabel(ax,'scaled coefficient 2-norm');
title(ax,'Coefficient growth');
grid(ax,'on');
axis(ax,'tight');
legend(ax,curveHandles,methodNames,'Location','best');
title(layout,'Single ellipsoid: monopoles versus monopoles and dipoles');

figure('Name','Laplace ellipsoid collocation grid');
scatter3(collocationPoints(:,1),collocationPoints(:,2), ...
    collocationPoints(:,3),14,[0.10,0.35,0.80],'filled');
axis equal;
axis tight;
view(34,24);
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('True-surface collocation grid: [%d,%d], M=%d', ...
    collocationResolution,size(collocationPoints,1)));

% Show logarithmic fields alone unless linear panels are explicitly enabled.
if show_linear
    displayedScaleIndices = 1:2;
else
    displayedScaleIndices = 2;
end
numberOfDisplayedScales = numel(displayedScaleIndices);

surfaceTriangles = convhull(testPoints(:,1),testPoints(:,2),testPoints(:,3));
figure('Name','Laplace residual on the ellipsoid surface');
layout = tiledlayout(numberOfDisplayedScales,numberOfMethods, ...
    'TileSpacing','compact','Padding','compact');
for methodIndex = 1:numberOfMethods
    surfaceResidual = finestTestResiduals(:,methodIndex);
    maximumAbsoluteResidual = max(abs(surfaceResidual));
    residualFloor = max(maximumAbsoluteResidual*1e-12,realmin('double'));
    surfaceFields = {surfaceResidual, ...
        log10(max(abs(surfaceResidual),residualFloor))};
    colourBarLabels = {'signed residual u_{MFS}-1', ...
        'log_{10}|u_{MFS}-1|'};
    scaleNames = {'linear','log_{10} magnitude'};

    for displayedScaleIndex = 1:numberOfDisplayedScales
        scaleIndex = displayedScaleIndices(displayedScaleIndex);
        tileIndex = (displayedScaleIndex-1)*numberOfMethods+methodIndex;
        ax = nexttile(layout,tileIndex);
        surfaceField = surfaceFields{scaleIndex};
        trisurf(surfaceTriangles,testPoints(:,1),testPoints(:,2), ...
            testPoints(:,3),surfaceField,'Parent',ax, ...
            'EdgeColor','none','FaceColor','interp');
        fieldLimits = [min(surfaceField),max(surfaceField)];
        if fieldLimits(1) < fieldLimits(2)
            % Every panel uses its own observed range; limits are not shared.
            clim(ax,fieldLimits);
        end
        colormap(ax,turbo(256));
        colourBar = colorbar(ax);
        colourBar.Label.String = colourBarLabels{scaleIndex};
        if fieldLimits(1) < fieldLimits(2)
            colourBar.Ticks = linspace(fieldLimits(1),fieldLimits(2),5);
        end
        axis(ax,'equal');
        axis(ax,'tight');
        view(ax,34,24);
        grid(ax,'on');
        xlabel(ax,'x');
        ylabel(ax,'y');
        zlabel(ax,'z');
        title(ax,sprintf('%s: %s', ...
            methodNames(methodIndex),scaleNames{scaleIndex}));
    end
end
maximumSurfaceResiduals = max(abs(finestTestResiduals),[],1);
if show_linear
    residualFigureTitle = sprintf( ...
        'Independent surface residual, finest focal grid N=%d; max |r| %.2e / %.2e', ...
        sourceNodeCounts(end),maximumSurfaceResiduals);
else
    logMaximumSurfaceResiduals = log10(max( ...
        maximumSurfaceResiduals,realmin('double')));
    residualFigureTitle = sprintf( ...
        ['Independent surface residual, finest focal grid N=%d; ' ...
         'log_{10}(max |r|) %.2f / %.2f'], ...
        sourceNodeCounts(end),logMaximumSurfaceResiduals);
end
title(layout,residualFigureTitle);

focalSemiAxes = sqrt(semiAxes(1:2).^2-semiAxes(3)^2);
boundaryAngle = linspace(0,2*pi,500);
focalBoundary = [focalSemiAxes(1)*cos(boundaryAngle); ...
    focalSemiAxes(2)*sin(boundaryAngle)];
ellipsoidCrossSection = [semiAxes(1)*cos(boundaryAngle); ...
    semiAxes(2)*sin(boundaryAngle)];
physicalDipoleMoments = dipoleLengthScale* ...
    finestCombinedScaledDipoleCoefficients;
dipoleMagnitudes = vecnorm(physicalDipoleMoments,2,2);

figure('Name','Laplace confocal source strengths');
layout = tiledlayout(numberOfDisplayedScales,3, ...
    'TileSpacing','compact','Padding','compact');
sourceMagnitudeSets = {abs(finestMonopoleOnlyCoefficients), ...
    abs(finestCombinedMonopoleCoefficients),dipoleMagnitudes};
sourceTitles = {'monopoles-only |q|', ...
    'combined-basis |q|','combined-basis |p|'};
sourceColourBarLabels = {'|q|','|q|','|p|'};
maximumDipoleMagnitude = max(dipoleMagnitudes);

for displayedScaleIndex = 1:numberOfDisplayedScales
    scaleIndex = displayedScaleIndices(displayedScaleIndex);
    for panelIndex = 1:3
        tileIndex = (displayedScaleIndex-1)*3+panelIndex;
        ax = nexttile(layout,tileIndex);
        sourceMagnitudes = sourceMagnitudeSets{panelIndex};
        maximumSourceMagnitude = max(sourceMagnitudes);
        magnitudeFloor = max(maximumSourceMagnitude*1e-12,realmin('double'));
        if scaleIndex == 1
            sourceField = sourceMagnitudes;
            scaleName = 'linear';
            colourBarLabel = sourceColourBarLabels{panelIndex};
        else
            sourceField = log10(max(sourceMagnitudes,magnitudeFloor));
            scaleName = 'log_{10}';
            colourBarLabel = ['log_{10}',sourceColourBarLabels{panelIndex}];
        end

        nodePlot = scatter(ax,finestSourcePoints(:,1), ...
            finestSourcePoints(:,2),36,sourceField,'filled');
        hold(ax,'on');
        focalPlot = plot(ax,focalBoundary(1,:),focalBoundary(2,:),'k-', ...
            'LineWidth',1.3);
        bodyPlot = plot(ax,ellipsoidCrossSection(1,:), ...
            ellipsoidCrossSection(2,:),'k--','LineWidth',1.0);
        fieldLimits = [min(sourceField),max(sourceField)];
        if fieldLimits(1) < fieldLimits(2)
            % Every source panel uses its own observed range.
            clim(ax,fieldLimits);
        end
        colormap(ax,turbo(256));

        if panelIndex == 3 && maximumDipoleMagnitude > 0
            displayScale = 0.12*max(focalSemiAxes)/maximumDipoleMagnitude;
            quiver(ax,finestSourcePoints(:,1),finestSourcePoints(:,2), ...
                displayScale*physicalDipoleMoments(:,1), ...
                displayScale*physicalDipoleMoments(:,2),0,'k', ...
                'LineWidth',0.8,'MaxHeadSize',0.8);
        end

        axis(ax,'equal');
        grid(ax,'on');
        xlabel(ax,'x');
        ylabel(ax,'y');
        title(ax,sprintf('%s: %s',sourceTitles{panelIndex},scaleName));
        colourBar = colorbar(ax);
        colourBar.Label.String = colourBarLabel;
        if fieldLimits(1) < fieldLimits(2)
            colourBar.Ticks = linspace(fieldLimits(1),fieldLimits(2),5);
        end
        if tileIndex == 1
            legend(ax,[nodePlot,focalPlot,bodyPlot], ...
                {'source node','focal boundary','ellipsoid cross-section'}, ...
                'Location','best');
        end
    end
end
if show_linear
    sourceFigureTitle = ...
        'Source magnitudes, top: linear, bottom: log_{10}';
else
    sourceFigureTitle = 'Source magnitudes, log_{10} scale';
end
title(layout,sprintf('%s; grid [%d,%d], N=%d',sourceFigureTitle, ...
    ellipseResolutions(end,:),sourceNodeCounts(end)));

result.semiAxes = semiAxes;
result.ellipseResolutions = ellipseResolutions;
result.sourceNodeCounts = sourceNodeCounts;
result.sourceDistributionType = sourceDistributionType;
result.showLinear = show_linear;
result.collocationResolution = collocationResolution;
result.collocationCount = size(collocationPoints,1);
result.testResolution = testResolution;
result.testCount = size(testPoints,1);
result.svdCutoff = svdCutoff;
result.dipoleLengthScale = dipoleLengthScale;
result.methodNames = methodNames;
result.exactCapacity = exactCapacity;
result.collocationMaxResidual = collocationMaxResidual;
result.testMaxResidual = testMaxResidual;
result.capacities = capacities;
result.capacityRelativeError = capacityRelativeError;
result.coefficientNorms = coefficientNorms;
result.monopoleCoefficientNorms = monopoleCoefficientNorms;
result.scaledDipoleCoefficientNorms = scaledDipoleCoefficientNorms;
result.retainedRanks = retainedRanks;
result.keptConditionNumbers = keptConditionNumbers;
result.assemblyTimes = assemblyTimes;
result.solveTimes = solveTimes;
result.finestSourcePoints = finestSourcePoints;
result.finestMonopoleOnlyCoefficients = finestMonopoleOnlyCoefficients;
result.finestCombinedMonopoleCoefficients = ...
    finestCombinedMonopoleCoefficients;
result.finestCombinedDipoleMoments = physicalDipoleMoments;
result.finestTestResiduals = finestTestResiduals;
end

function D = laplace_dipole_potential_matrix(targets,sources)
%LAPLACE_DIPOLE_POTENTIAL_MATRIX Build Cartesian source-dipole columns.
dx = targets(:,1)-sources(:,1).';
dy = targets(:,2)-sources(:,2).';
dz = targets(:,3)-sources(:,3).';
distanceSquared = dx.^2+dy.^2+dz.^2;
if any(distanceSquared == 0,'all')
    error('laplace_dipole_potential_matrix:CoincidentPoints', ...
        'A target coincides with a source.');
end
inverseDistanceCubed = distanceSquared.^(-3/2);
D = [dx.*inverseDistanceCubed,dy.*inverseDistanceCubed, ...
    dz.*inverseDistanceCubed]/(4*pi);
end

function [coefficients,retainedRank,keptConditionNumber] = ...
    truncated_svd_solve(A,rhs,svdCutoff)
%TRUNCATED_SVD_SOLVE Apply the experiment's common relative SVD cutoff.
[leftVectors,singularMatrix,rightVectors] = svd(A,'econ');
singularValues = diag(singularMatrix);
retained = singularValues >= svdCutoff*singularValues(1);
if ~any(retained)
    error('truncated_svd_solve:EmptyTruncation', ...
        'The SVD cutoff discarded every singular value.');
end
projectedData = leftVectors(:,retained).'*rhs;
coefficients = rightVectors(:,retained)*( ...
    projectedData./singularValues(retained));
retainedRank = nnz(retained);
keptConditionNumber = singularValues(1)/ ...
    singularValues(find(retained,1,'last'));
end
