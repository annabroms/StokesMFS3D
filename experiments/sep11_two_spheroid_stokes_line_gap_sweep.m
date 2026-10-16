function result = sep11_two_spheroid_stokes_line_gap_sweep(deltas,visualise)
%SEP11_TWO_SPHEROID_STOKES_LINE_GAP_SWEEP Sweep the gap in three pairs.
%
% Two equal prolate spheroids have distinct prescribed translations and no
% rotation.  Each body has 40 source nodes on its focal line.  An S-only
% baseline is compared with a user-selected combination of:
%
%   S: three-component Stokeslets;
%   R: three-component rotlets;
%   T: two stresslet/DLP blocks along the body's transverse local axes;
%   D: three-component potential dipoles.
%
% R and T columns are multiplied by the minor semiaxis, and D columns by
% its square, so all fitted vector coefficients have nominal force units.
% For tip-to-tip, side-by-side, and non-symmetric oblique pairs, the
% experiment plots the independent maximum relative velocity residual and
% maximum scaled nodal coefficient magnitude against the closest gap.
%
%   result = sep11_two_spheroid_stokes_line_gap_sweep
%   result = sep11_two_spheroid_stokes_line_gap_sweep(deltas,false)

if nargin < 1 || isempty(deltas)
    deltas = [1e-3,3e-3,1e-2,3e-2,0.1,0.3,1,10,100,1000];
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
addpath(fullfile(projectDirectory,'common','dense'));

% ----------------------- USER SETTINGS -------------------------------
translationVelocities = [1,0,0; -0.35,0.55,0.2];
selectedSingularityTypes = ["S","T","R","D"];
% ---------------------------------------------------------------------

semiAxes = [1.6;0.7;0.7];
numberOfLineSourcesPerBody = 80;
configurationNames = ["tip","side","oblique"];
collocationResolution = [48,36];
testResolution = [55,37];
svdCutoff = 1e-12;
singularityLengthScale = semiAxes(2);
familyNames = ["S","R","T","D"];

selectedSingularityTypes = upper(string(selectedSingularityTypes(:).'));
if isempty(selectedSingularityTypes) || ...
        any(ismissing(selectedSingularityTypes)) || ...
        any(strlength(selectedSingularityTypes) == 0) || ...
        any(~ismember(selectedSingularityTypes,familyNames))
    error(['selectedSingularityTypes must be a nonempty subset of ', ...
        '"S", "R", "T", and "D".']);
end
if numel(unique(selectedSingularityTypes)) ~= ...
        numel(selectedSingularityTypes)
    error('selectedSingularityTypes must not contain duplicates.');
end
selectedFamilies = ismember(familyNames,selectedSingularityTypes);
selectedCombinationName = strjoin(familyNames(selectedFamilies),'+');
methodNames = ["S only","selected: "+selectedCombinationName];
activeFamiliesByMethod = [true,false,false,false; selectedFamilies];

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
maximumFamilyCoefficientMagnitude = zeros( ...
    numberOfDeltas,numberOfConfigurations,numberOfMethods, ...
    numel(familyNames));
coefficientNorms = zeros(metricSize);
retainedRanks = zeros(metricSize);
keptConditionNumbers = zeros(metricSize);
assemblyTimes = zeros(numberOfDeltas,numberOfConfigurations);
solveTimes = zeros(metricSize);
testEvaluationTimes = zeros(metricSize);
coefficients = cell(metricSize);
geometries = cell(numberOfDeltas,numberOfConfigurations);
velocityScale = max(vecnorm(translationVelocities,2,2));

numberOfSources = 2*numberOfLineSourcesPerBody;
numberOfCollocationPoints = 2*numberOfCollocationPointsPerBody;
familyColumnsPerNode = [3,3,6,3];
familyColumnCounts = numberOfSources*familyColumnsPerNode;
unknownCounts = activeFamiliesByMethod*familyColumnCounts.';
assert(3*numberOfCollocationPoints > max(unknownCounts), ...
    'Both Stokes systems must be overdetermined.');

fprintf('\nTwo-spheroid Stokes focal-line gap sweep\n');
fprintf('  semiaxes                         [%.2f, %.2f, %.2f]\n', ...
    semiAxes);
fprintf('  body 1 translation               [%g, %g, %g]^T\n', ...
    translationVelocities(1,:));
fprintf('  body 2 translation               [%g, %g, %g]^T\n', ...
    translationVelocities(2,:));
fprintf('  focal-line sources per body      %d (Chebyshev--Lobatto)\n', ...
    numberOfLineSourcesPerBody);
fprintf('  baseline                         S only\n');
fprintf('  selected combination             %s\n',selectedCombinationName);
fprintf('  T directions per node            2 transverse directions\n');
fprintf('  singularity length scale         %.3g\n', ...
    singularityLengthScale);
fprintf('  collocation/test nodes per body  %d / %d\n', ...
    numberOfCollocationPointsPerBody,numberOfTestPointsPerBody);
fprintf('  SVD cutoff                       %.1e\n',svdCutoff);
fprintf('  relative residual speed scale    %.3g\n\n',velocityScale);

for configurationIndex = 1:numberOfConfigurations
    caseName = configurationNames(configurationIndex);
    fprintf('%s configuration\n',caseName);
    for deltaIndex = 1:numberOfDeltas
        delta = deltas(deltaIndex);
        pair = two_spheroid_gap_configuration(caseName,semiAxes,delta);
        geometries{deltaIndex,configurationIndex} = pair;

        [sourcePoints,collocationPoints,testPoints,transverseDirections] = ...
            transform_pair_points(pair,localSourcePoints, ...
            localCollocationPoints,localTestPoints);
        collocationVelocity = prescribed_translation_data( ...
            numberOfCollocationPointsPerBody,translationVelocities);
        testVelocity = prescribed_translation_data( ...
            numberOfTestPointsPerBody,translationVelocities);

        assemblyStart = tic;
        collocationSlp = stokes_SLP_mat(sourcePoints,collocationPoints);
        testSlp = stokes_SLP_mat(sourcePoints,testPoints);
        [collocationRotlets,~,collocationPotentialDipoles] = ...
            stokes_singularity_matrices(sourcePoints,collocationPoints, ...
            [0,1,0],logical([selectedFamilies(2),false, ...
            selectedFamilies(4)]));
        [testRotlets,~,testPotentialDipoles] = ...
            stokes_singularity_matrices(sourcePoints,testPoints, ...
            [0,1,0],logical([selectedFamilies(2),false, ...
            selectedFamilies(4)]));

        collocationStresslets = zeros(3*numberOfCollocationPoints,0);
        testStresslets = zeros(3*size(testPoints,1),0);
        if selectedFamilies(3)
            collocationStresslets = [stokes_DLP_mat( ...
                sourcePoints,collocationPoints,transverseDirections{1}), ...
                stokes_DLP_mat(sourcePoints,collocationPoints, ...
                transverseDirections{2})];
            testStresslets = [stokes_DLP_mat( ...
                sourcePoints,testPoints,transverseDirections{1}), ...
                stokes_DLP_mat(sourcePoints,testPoints, ...
                transverseDirections{2})];
        end
        assemblyTimes(deltaIndex,configurationIndex) = toc(assemblyStart);

        % Canonical family order is [S,R,T,D].  Every vector block uses
        % interleaved [x1;y1;z1;x2;...] component ordering.
        collocationFamilyMatrices = {collocationSlp, ...
            singularityLengthScale*collocationRotlets, ...
            singularityLengthScale*collocationStresslets, ...
            singularityLengthScale^2*collocationPotentialDipoles};
        testFamilyMatrices = {testSlp, ...
            singularityLengthScale*testRotlets, ...
            singularityLengthScale*testStresslets, ...
            singularityLengthScale^2*testPotentialDipoles};
        collocationMatrices = cell(1,numberOfMethods);
        testMatrices = cell(1,numberOfMethods);
        for methodIndex = 1:numberOfMethods
            activeFamilies = activeFamiliesByMethod(methodIndex,:);
            collocationMatrices{methodIndex} = horzcat( ...
                collocationFamilyMatrices{activeFamilies});
            testMatrices{methodIndex} = horzcat( ...
                testFamilyMatrices{activeFamilies});
        end

        fprintf('  delta %.3e\n',delta);
        for methodIndex = 1:numberOfMethods
            systemMatrix = collocationMatrices{methodIndex};
            assert(isequal(size(systemMatrix), ...
                [3*numberOfCollocationPoints,unknownCounts(methodIndex)]));

            solveStart = tic;
            [currentCoefficients,retainedRank,keptConditionNumber] = ...
                truncated_svd_solve( ...
                systemMatrix,collocationVelocity,svdCutoff);
            solveTimes(deltaIndex,configurationIndex,methodIndex) = ...
                toc(solveStart);

            collocationError = systemMatrix*currentCoefficients- ...
                collocationVelocity;
            testStart = tic;
            testError = testMatrices{methodIndex}*currentCoefficients- ...
                testVelocity;
            testEvaluationTimes( ...
                deltaIndex,configurationIndex,methodIndex) = toc(testStart);
            collocationErrorVectors = reshape( ...
                collocationError,3,[]).';
            testErrorVectors = reshape(testError,3,[]).';

            activeFamilies = activeFamiliesByMethod(methodIndex,:);
            [nodalMagnitudes,nodalMagnitudesByBody,familyMaximum] = ...
                stokes_nodal_coefficient_magnitudes( ...
                currentCoefficients,numberOfSources, ...
                numberOfLineSourcesPerBody,activeFamilies, ...
                familyColumnsPerNode);
            collocationMaxRelativeResidual( ...
                deltaIndex,configurationIndex,methodIndex) = ...
                max(vecnorm(collocationErrorVectors,2,2))/velocityScale;
            testMaxRelativeResidual( ...
                deltaIndex,configurationIndex,methodIndex) = ...
                max(vecnorm(testErrorVectors,2,2))/velocityScale;
            maximumCoefficientMagnitude( ...
                deltaIndex,configurationIndex,methodIndex) = ...
                max(nodalMagnitudes);
            maximumCoefficientMagnitudeByBody( ...
                deltaIndex,:,configurationIndex,methodIndex) = ...
                nodalMagnitudesByBody;
            maximumFamilyCoefficientMagnitude( ...
                deltaIndex,configurationIndex,methodIndex,:) = ...
                familyMaximum;
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
        testMaxRelativeResidual,maximumCoefficientMagnitude, ...
        selectedCombinationName);
end

result.problem = "Stokes resistance with prescribed translations";
result.semiAxes = semiAxes;
result.deltas = deltas;
result.configurationNames = configurationNames;
result.configurationLabels = configurationLabels;
result.methodNames = methodNames;
result.familyNames = familyNames;
result.selectedSingularityTypes = selectedSingularityTypes;
result.selectedCombinationName = selectedCombinationName;
result.selectedFamilies = selectedFamilies;
result.activeFamiliesByMethod = activeFamiliesByMethod;
result.familyColumnsPerNode = familyColumnsPerNode;
result.translationVelocities = translationVelocities;
result.angularVelocities = zeros(2,3);
result.numberOfLineSourcesPerBody = numberOfLineSourcesPerBody;
result.sourceDistribution = "Chebyshev-Lobatto on each full focal line";
result.focalHalfLength = focalHalfLength;
result.normalizedLineCoordinates = normalizedLineCoordinates;
result.stressletDirectionDescription = ...
    "both local directions transverse to each focal line";
result.singularityLengthScale = singularityLengthScale;
result.collocationResolution = collocationResolution;
result.collocationCountPerBody = numberOfCollocationPointsPerBody;
result.testResolution = testResolution;
result.testCountPerBody = numberOfTestPointsPerBody;
result.residualNormalization = velocityScale;
result.svdCutoff = svdCutoff;
result.unknownCounts = unknownCounts;
result.collocationMaxRelativeResidual = ...
    collocationMaxRelativeResidual;
result.testMaxRelativeResidual = testMaxRelativeResidual;
result.maximumCoefficientMagnitude = maximumCoefficientMagnitude;
result.maximumCoefficientMagnitudeByBody = ...
    maximumCoefficientMagnitudeByBody;
result.maximumFamilyCoefficientMagnitude = ...
    maximumFamilyCoefficientMagnitude;
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


function [sources,collocation,testPoints,transverseDirections] = ...
        transform_pair_points( ...
        pair,localSources,localCollocation,localTestPoints)
%TRANSFORM_PAIR_POINTS Place nodes and both transverse stresslet frames.
numberOfSourcesPerBody = size(localSources,1);
numberOfCollocationPerBody = size(localCollocation,1);
numberOfTestPointsPerBody = size(localTestPoints,1);
sources = zeros(2*numberOfSourcesPerBody,3);
collocation = zeros(2*numberOfCollocationPerBody,3);
testPoints = zeros(2*numberOfTestPointsPerBody,3);
transverseDirections = {zeros(2*numberOfSourcesPerBody,3), ...
    zeros(2*numberOfSourcesPerBody,3)};

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
    transverseDirections{1}(sourceRows,:) = repmat( ...
        rotation(:,2).',numberOfSourcesPerBody,1);
    transverseDirections{2}(sourceRows,:) = repmat( ...
        rotation(:,3).',numberOfSourcesPerBody,1);
end
end


function data = prescribed_translation_data( ...
        numberOfPointsPerBody,translationVelocities)
%PRESCRIBED_TRANSLATION_DATA Stack interleaved velocities body by body.
pointVelocities = [ ...
    repmat(translationVelocities(1,:),numberOfPointsPerBody,1); ...
    repmat(translationVelocities(2,:),numberOfPointsPerBody,1)];
data = reshape(pointVelocities.',[],1);
end


function [magnitudes,maximumByBody,familyMaximum] = ...
        stokes_nodal_coefficient_magnitudes( ...
        coefficients,numberOfSources,numberOfSourcesPerBody, ...
        activeFamilies,familyColumnsPerNode)
%STOKES_NODAL_COEFFICIENT_MAGNITUDES Combine every active family by node.
squaredMagnitudes = zeros(numberOfSources,1);
familyMaximum = zeros(1,numel(activeFamilies));
coefficientOffset = 0;
for familyIndex = 1:numel(activeFamilies)
    if ~activeFamilies(familyIndex)
        continue
    end
    directionsPerNode = familyColumnsPerNode(familyIndex)/3;
    familyBlockLength = 3*numberOfSources*directionsPerNode;
    familyIndices = coefficientOffset+(1:familyBlockLength);
    familyBlock = coefficients(familyIndices);
    familyVectorMagnitudes = zeros(numberOfSources,directionsPerNode);
    for directionIndex = 1:directionsPerNode
        directionOffset = (directionIndex-1)*3*numberOfSources;
        directionIndices = directionOffset+(1:3*numberOfSources);
        nodalVectors = reshape( ...
            familyBlock(directionIndices),3,numberOfSources).';
        familyVectorMagnitudes(:,directionIndex) = ...
            vecnorm(nodalVectors,2,2);
    end
    squaredMagnitudes = squaredMagnitudes+ ...
        sum(familyVectorMagnitudes.^2,2);
    familyMaximum(familyIndex) = max(familyVectorMagnitudes,[],'all');
    coefficientOffset = coefficientOffset+familyBlockLength;
end
assert(coefficientOffset == numel(coefficients));
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
        residuals,maxCoefficients,selectedCombinationName)
%PLOT_GAP_SWEEP Visualize both bases separately for every configuration.
numberOfConfigurations = numel(configurationLabels);
numberOfMethods = numel(methodNames);
methodStyles = {'-o','--s'};
figure('Name','Two-spheroid Stokes focal-line gap sweep');
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
title(layout,sprintf( ...
    'Two spheroids: Stokes S only versus %s, 40 nodes per body', ...
    selectedCombinationName));
end
