function result = sep11_two_spheroid_laplace_line_gap_sweep(deltas,visualise)
%SEP11_TWO_SPHEROID_LAPLACE_LINE_GAP_SWEEP Sweep the gap in three pairs.
%
% Two equal prolate spheroids are held at distinct voltages.  Each body has
% 40 source nodes on its focal line.  Two bases are compared:
%
%   1. SLP: one scalar monopole per source node;
%   2. SLP+DLP: a monopole and one transverse point dipole per node.
%
% The DLP columns are multiplied by the minor semiaxis, so their fitted
% coefficients have the same nominal charge units as the SLP coefficients.
% The DLP direction is the second local body axis, perpendicular to the
% focal line, and is rotated with each spheroid.
% For tip-to-tip, side-by-side, and non-symmetric oblique configurations,
% the experiment plots the independent maximum relative residual and the
% maximum scaled nodal coefficient magnitude against the closest gap delta.
%
%   result = sep11_two_spheroid_laplace_line_gap_sweep
%   result = sep11_two_spheroid_laplace_line_gap_sweep(deltas,false)

if nargin < 1 || isempty(deltas)
    deltas = [1e-3,3e-3,1e-2,3e-2,0.1,0.3,1, 10, 100,1000];
end
if nargin < 2 || isempty(visualise)
    visualise = true;
end

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory);
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'laplace'));

semiAxes = [1.6;0.7;0.7];
numberOfLineSourcesPerBody = 80;
bodyVoltages = [1;-0.35];
configurationNames = ["tip","side","oblique"];
methodNames = ["SLP only","SLP + DLP"];
collocationResolution = [48,36];
testResolution = [55,37];
svdCutoff = 1e-12;
dipoleLengthScale = semiAxes(2);
dlpLocalDirection = [0;1;0];

deltas = validate_gap_vector(deltas);
if ~isscalar(visualise)
    error('visualise must be a scalar logical value.');
end
visualise = logical(visualise);

surface = ellipsoid_param(semiAxes(1),semiAxes(2),semiAxes(3));
collocationSurface = setupsurfquad(surface,collocationResolution);
testSurface = setupsurfquad(surface,testResolution);
localCollocationPoints = collocationSurface.x.';
localTestPoints = testSurface.x.';
numberOfCollocationPointsPerBody = size(localCollocationPoints,1);
numberOfTestPointsPerBody = size(localTestPoints,1);

focalHalfLength = sqrt(semiAxes(1)^2-semiAxes(2)^2);
lineParameters = linspace(0,1,numberOfLineSourcesPerBody).';
normalizedLineCoordinates = -cos(pi*lineParameters);
localSourcePoints = [focalHalfLength*normalizedLineCoordinates, ...
    zeros(numberOfLineSourcesPerBody,2)];

numberOfDeltas = numel(deltas);
numberOfConfigurations = numel(configurationNames);
numberOfMethods = numel(methodNames);
metricSize = [numberOfDeltas,numberOfConfigurations,numberOfMethods];
collocationMaxRelativeResidual = zeros(metricSize);
testMaxRelativeResidual = zeros(metricSize);
maximumCoefficientMagnitude = zeros(metricSize);
maximumCoefficientMagnitudeByBody = zeros( ...
    numberOfDeltas,2,numberOfConfigurations,numberOfMethods);
coefficientNorms = zeros(metricSize);
retainedRanks = zeros(metricSize);
keptConditionNumbers = zeros(metricSize);
assemblyTimes = zeros(numberOfDeltas,numberOfConfigurations);
solveTimes = zeros(metricSize);
testEvaluationTimes = zeros(metricSize);
coefficients = cell(metricSize);
geometries = cell(numberOfDeltas,numberOfConfigurations);
potentialScale = max(norm(bodyVoltages,inf),eps);

numberOfSources = 2*numberOfLineSourcesPerBody;
numberOfCollocationPoints = 2*numberOfCollocationPointsPerBody;
unknownCounts = [numberOfSources,2*numberOfSources];
assert(numberOfCollocationPoints > max(unknownCounts), ...
    'Both Laplace systems must be overdetermined.');

fprintf('\nTwo-spheroid Laplace focal-line gap sweep\n');
fprintf('  semiaxes                         [%.2f, %.2f, %.2f]\n', ...
    semiAxes);
fprintf('  body voltages                    [%g, %g]^T\n',bodyVoltages);
fprintf('  focal-line sources per body      %d (Chebyshev--Lobatto)\n', ...
    numberOfLineSourcesPerBody);
fprintf('  methods                          SLP / SLP+DLP at every node\n');
fprintf('  DLP direction per node           local e_2 (perpendicular to line)\n');
fprintf('  dipole column length scale       %.3g\n',dipoleLengthScale);
fprintf('  collocation/test nodes per body  %d / %d\n', ...
    numberOfCollocationPointsPerBody,numberOfTestPointsPerBody);
fprintf('  SVD cutoff                       %.1e\n',svdCutoff);
fprintf('  relative residual scale          %.3g\n\n',potentialScale);

for configurationIndex = 1:numberOfConfigurations
    caseName = configurationNames(configurationIndex);
    fprintf('%s configuration\n',caseName);
    for deltaIndex = 1:numberOfDeltas
        delta = deltas(deltaIndex);
        pair = two_spheroid_gap_configuration(caseName,semiAxes,delta);
        geometries{deltaIndex,configurationIndex} = pair;

        [sourcePoints,collocationPoints,testPoints,dlpDirections] = ...
            transform_pair_points(pair,localSourcePoints, ...
            localCollocationPoints,localTestPoints,dlpLocalDirection);
        collocationData = repelem( ...
            bodyVoltages,numberOfCollocationPointsPerBody);
        testData = repelem(bodyVoltages,numberOfTestPointsPerBody);

        assemblyStart = tic;
        collocationSlp = lap3d_charge_potential_matrix( ...
            collocationPoints,sourcePoints);
        collocationDlp = laplace_point_dlp_matrix( ...
            collocationPoints,sourcePoints,dlpDirections);
        testSlp = lap3d_charge_potential_matrix(testPoints,sourcePoints);
        testDlp = laplace_point_dlp_matrix( ...
            testPoints,sourcePoints,dlpDirections);
        assemblyTimes(deltaIndex,configurationIndex) = toc(assemblyStart);

        % Columns are [q] for SLP and [q;d_perp] for SLP+DLP.
        collocationMatrices = {collocationSlp, ...
            [collocationSlp,dipoleLengthScale*collocationDlp]};
        testMatrices = {testSlp, ...
            [testSlp,dipoleLengthScale*testDlp]};

        fprintf('  delta %.3e\n',delta);
        for methodIndex = 1:numberOfMethods
            systemMatrix = collocationMatrices{methodIndex};
            assert(isequal(size(systemMatrix), ...
                [numberOfCollocationPoints,unknownCounts(methodIndex)]));

            solveStart = tic;
            [currentCoefficients,retainedRank,keptConditionNumber] = ...
                truncated_svd_solve(systemMatrix,collocationData,svdCutoff);
            solveTimes(deltaIndex,configurationIndex,methodIndex) = ...
                toc(solveStart);

            collocationError = systemMatrix*currentCoefficients- ...
                collocationData;
            testStart = tic;
            testError = testMatrices{methodIndex}*currentCoefficients- ...
                testData;
            testEvaluationTimes( ...
                deltaIndex,configurationIndex,methodIndex) = toc(testStart);

            [nodalMagnitudes,nodalMagnitudesByBody] = ...
                laplace_nodal_coefficient_magnitudes( ...
                currentCoefficients,numberOfSources, ...
                numberOfLineSourcesPerBody,methodIndex == 2);
            collocationMaxRelativeResidual( ...
                deltaIndex,configurationIndex,methodIndex) = ...
                norm(collocationError,inf)/potentialScale;
            testMaxRelativeResidual( ...
                deltaIndex,configurationIndex,methodIndex) = ...
                norm(testError,inf)/potentialScale;
            maximumCoefficientMagnitude( ...
                deltaIndex,configurationIndex,methodIndex) = ...
                max(nodalMagnitudes);
            maximumCoefficientMagnitudeByBody( ...
                deltaIndex,:,configurationIndex,methodIndex) = ...
                nodalMagnitudesByBody;
            coefficientNorms(deltaIndex,configurationIndex,methodIndex) = ...
                norm(currentCoefficients);
            retainedRanks(deltaIndex,configurationIndex,methodIndex) = ...
                retainedRank;
            keptConditionNumbers( ...
                deltaIndex,configurationIndex,methodIndex) = ...
                keptConditionNumber;
            coefficients{deltaIndex,configurationIndex,methodIndex} = ...
                currentCoefficients;

            fprintf(['    %-10s: test max rel residual %.3e, ', ...
                'max |c_j| %.3e, rank %d\n'],methodNames(methodIndex), ...
                testMaxRelativeResidual( ...
                deltaIndex,configurationIndex,methodIndex), ...
                maximumCoefficientMagnitude( ...
                deltaIndex,configurationIndex,methodIndex),retainedRank);
        end
    end
    fprintf('\n');
end

configurationLabels = strings(1,numberOfConfigurations);
for configurationIndex = 1:numberOfConfigurations
    configurationLabels(configurationIndex) = ...
        geometries{1,configurationIndex}.label;
end

if visualise
    plot_gap_sweep(deltas,configurationLabels,methodNames, ...
        testMaxRelativeResidual,maximumCoefficientMagnitude);
end

result.problem = "exterior Laplace Dirichlet";
result.semiAxes = semiAxes;
result.deltas = deltas;
result.configurationNames = configurationNames;
result.configurationLabels = configurationLabels;
result.methodNames = methodNames;
result.bodyVoltages = bodyVoltages;
result.numberOfLineSourcesPerBody = numberOfLineSourcesPerBody;
result.sourceDistribution = "Chebyshev-Lobatto on each full focal line";
result.focalHalfLength = focalHalfLength;
result.normalizedLineCoordinates = normalizedLineCoordinates;
result.collocationResolution = collocationResolution;
result.collocationCountPerBody = numberOfCollocationPointsPerBody;
result.testResolution = testResolution;
result.testCountPerBody = numberOfTestPointsPerBody;
result.residualNormalization = potentialScale;
result.svdCutoff = svdCutoff;
result.dipoleLengthScale = dipoleLengthScale;
result.dlpLocalDirection = dlpLocalDirection;
result.dlpDirectionDescription = ...
    "second local body axis, perpendicular to each focal line";
result.unknownCounts = unknownCounts;
result.collocationMaxRelativeResidual = ...
    collocationMaxRelativeResidual;
result.testMaxRelativeResidual = testMaxRelativeResidual;
result.maximumCoefficientMagnitude = maximumCoefficientMagnitude;
result.maximumCoefficientMagnitudeByBody = ...
    maximumCoefficientMagnitudeByBody;
result.coefficientNorms = coefficientNorms;
result.retainedRanks = retainedRanks;
result.keptConditionNumbers = keptConditionNumbers;
result.assemblyTimes = assemblyTimes;
result.solveTimes = solveTimes;
result.testEvaluationTimes = testEvaluationTimes;
result.coefficients = coefficients;
result.geometries = geometries;
end


function deltas = validate_gap_vector(deltas)
%VALIDATE_GAP_VECTOR Normalize and check the requested physical gaps.
if ~isnumeric(deltas) || ~isreal(deltas) || isempty(deltas) || ...
        ~isvector(deltas) || any(~isfinite(deltas)) || any(deltas <= 0)
    error('deltas must be a nonempty vector of positive finite numbers.');
end
deltas = sort(deltas(:));
end


function [sources,collocation,testPoints,dlpDirections] = ...
        transform_pair_points(pair,localSources,localCollocation, ...
        localTestPoints,dlpLocalDirection)
%TRANSFORM_PAIR_POINTS Place nodes and the transverse DLP directions.
numberOfSourcesPerBody = size(localSources,1);
numberOfCollocationPerBody = size(localCollocation,1);
numberOfTestPointsPerBody = size(localTestPoints,1);
sources = zeros(2*numberOfSourcesPerBody,3);
collocation = zeros(2*numberOfCollocationPerBody,3);
testPoints = zeros(2*numberOfTestPointsPerBody,3);
dlpDirections = zeros(2*numberOfSourcesPerBody,3);

for bodyIndex = 1:2
    sourceRows = (bodyIndex-1)*numberOfSourcesPerBody+ ...
        (1:numberOfSourcesPerBody);
    collocationRows = (bodyIndex-1)*numberOfCollocationPerBody+ ...
        (1:numberOfCollocationPerBody);
    testRows = (bodyIndex-1)*numberOfTestPointsPerBody+ ...
        (1:numberOfTestPointsPerBody);
    rotation = pair.rotations{bodyIndex};
    centre = pair.centres(:,bodyIndex).';
    sources(sourceRows,:) = localSources*rotation.'+centre;
    collocation(collocationRows,:) = localCollocation*rotation.'+centre;
    testPoints(testRows,:) = localTestPoints*rotation.'+centre;
    globalDlpDirection = rotation*dlpLocalDirection;
    assert(abs(globalDlpDirection.'*rotation(:,1)) < 100*eps);
    dlpDirections(sourceRows,:) = repmat( ...
        globalDlpDirection.',numberOfSourcesPerBody,1);
end
end


function matrix = laplace_point_dlp_matrix(targets,sources,directions)
%LAPLACE_POINT_DLP_MATRIX Build one directed point-DLP per source node.
assert(isequal(size(directions),size(sources)));
dx = targets(:,1)-sources(:,1).';
dy = targets(:,2)-sources(:,2).';
dz = targets(:,3)-sources(:,3).';
distanceSquared = dx.^2+dy.^2+dz.^2;
if any(distanceSquared == 0,'all')
    error('laplace_point_dlp_matrix:CoincidentPoints', ...
        'A target coincides with a source.');
end
inverseDistanceCubed = distanceSquared.^(-3/2);
directionProjection = dx.*directions(:,1).'+ ...
    dy.*directions(:,2).'+dz.*directions(:,3).';
matrix = directionProjection.*inverseDistanceCubed/(4*pi);
end


function [magnitudes,maximumByBody] = ...
        laplace_nodal_coefficient_magnitudes( ...
        coefficients,numberOfSources,numberOfSourcesPerBody,hasDlp)
%LAPLACE_NODAL_COEFFICIENT_MAGNITUDES Combine scaled blocks at each node.
monopoles = coefficients(1:numberOfSources);
squaredMagnitudes = abs(monopoles).^2;
if hasDlp
    dipoles = coefficients(numberOfSources+1:end);
    assert(numel(dipoles) == numberOfSources);
    squaredMagnitudes = squaredMagnitudes+abs(dipoles).^2;
end
magnitudes = sqrt(squaredMagnitudes);
maximumByBody = zeros(1,2);
for bodyIndex = 1:2
    bodyIndices = (bodyIndex-1)*numberOfSourcesPerBody+ ...
        (1:numberOfSourcesPerBody);
    maximumByBody(bodyIndex) = max(magnitudes(bodyIndices));
end
end


function [coefficients,retainedRank,keptConditionNumber] = ...
        truncated_svd_solve(matrix,rhs,svdCutoff)
%TRUNCATED_SVD_SOLVE Apply one relative truncation rule to every gap.
[leftVectors,singularMatrix,rightVectors] = svd(matrix,'econ');
singularValues = diag(singularMatrix);
retained = singularValues >= svdCutoff*singularValues(1);
if ~any(retained)
    error('truncated_svd_solve:EmptyTruncation', ...
        'The SVD cutoff discarded every singular value.');
end
coefficients = rightVectors(:,retained)*( ...
    (leftVectors(:,retained).'*rhs)./singularValues(retained));
retainedRank = nnz(retained);
keptConditionNumber = singularValues(1)/ ...
    singularValues(find(retained,1,'last'));
end


function plot_gap_sweep(deltas,configurationLabels,methodNames, ...
        residuals,maxCoefficients)
%PLOT_GAP_SWEEP Visualize both bases separately for every configuration.
numberOfConfigurations = numel(configurationLabels);
numberOfMethods = numel(methodNames);
methodStyles = {'-o','--s'};
figure('Name','Two-spheroid Laplace focal-line gap sweep');
layout = tiledlayout(2,numberOfConfigurations, ...
    'TileSpacing','compact','Padding','compact');

for configurationIndex = 1:numberOfConfigurations
    ax = nexttile(layout,configurationIndex);
    for methodIndex = 1:numberOfMethods
        loglog(ax,deltas,residuals(:,configurationIndex,methodIndex), ...
            methodStyles{methodIndex},'LineWidth',1.5,'MarkerSize',5);
        hold(ax,'on');
    end
    xlabel(ax,'surface gap \delta');
    ylabel(ax,'independent max relative residual');
    title(ax,configurationLabels(configurationIndex));
    grid(ax,'on');

    ax = nexttile(layout,numberOfConfigurations+configurationIndex);
    for methodIndex = 1:numberOfMethods
        loglog(ax,deltas,maxCoefficients( ...
            :,configurationIndex,methodIndex),methodStyles{methodIndex}, ...
            'LineWidth',1.5,'MarkerSize',5);
        hold(ax,'on');
    end
    xlabel(ax,'surface gap \delta');
    ylabel(ax,'max scaled nodal coefficient');
    grid(ax,'on');
    if configurationIndex == numberOfConfigurations
        legend(ax,methodNames,'Location','southwest');
    end
end
title(layout,'Two spheroids: Laplace SLP versus SLP+DLP, 40 nodes per body');
end
