function result = sep11_single_spheroid_stokes_line_singularity_combination_sweep
%SEP11_SINGLE_SPHEROID_STOKES_LINE_SINGULARITY_COMBINATION_SWEEP Compare bases.
%
% Solve one prolate-spheroid resistance problem with singularities placed
% at every node of a Chebyshev--Lobatto focal-line grid.  Every tested basis
% contains Stokeslets.  The sweep selects all requested combinations of:
%
%   R       three-component rotlets;
%   T_perp  stresslets with both independent normals transverse to the line;
%   D       three-component potential dipoles.
%
% The unit-viscosity normalizations are 1/(8*pi) for rotlets and
% stresslets and 1/(4*pi) for potential dipoles.  Rotlet and stresslet
% columns are multiplied by singularityLengthScale, and potential-dipole
% columns by singularityLengthScale^2.  All fitted coefficient blocks then
% have nominal force units, making their Euclidean norms comparable.

% ----------------------- USER SETTINGS -------------------------------
translationVelocity = [1;0;0];
angularVelocity = [0;0;0];
useTooShortSourceLine = false;
tooShortSourceLineFraction = 0.8;
numberOfLineSources = 32;

% Rows specify combinations; columns are [R, T_perp, D].  Stokeslets are
% always present.  Remove rows to run only selected combinations.
combinationMasks = logical([ ...
    0,0,0; ...
    1,0,0; ...
    0,1,0; ...
    0,0,1; ...
    1,1,0; ...
    1,0,1; ...
    0,1,1; ...
    1,1,1]);

% The focal line is the x-axis.  These two unit directions span its normal
% plane and are included together whenever T_perp is selected.
stressletNormalDirections = [0,1,0; 0,0,1];
% ---------------------------------------------------------------------

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory);
addpath(fullfile(projectDirectory,'images'));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'common','dense'));

% Other user-adjustable experiment parameters.
semiAxes = [1.2,0.7,0.7];
collocationResolution = [80,60];
testResolution = [83,55];
svdCutoff = 1e-12;
singularityLengthScale = semiAxes(2);

% Edit this map to change the line-node distribution.  It maps uniformly
% spaced parameters in [0,1] to ordered coordinates in [-1,1].
sourceDistributionName = "Chebyshev-Lobatto";
sourceCoordinateMap = @(parameter) -cos(pi*parameter);

if ~isnumeric(translationVelocity) || numel(translationVelocity) ~= 3 || ...
        any(~isfinite(translationVelocity))
    error('translationVelocity must be a finite numeric 3-vector.');
end
if ~isnumeric(angularVelocity) || numel(angularVelocity) ~= 3 || ...
        any(~isfinite(angularVelocity))
    error('angularVelocity must be a finite numeric 3-vector.');
end
translationVelocity = translationVelocity(:);
angularVelocity = angularVelocity(:);
rigidMotion = [translationVelocity;angularVelocity];
if norm(rigidMotion) == 0
    error('At least one prescribed velocity component must be nonzero.');
end

validateattributes(numberOfLineSources,{'numeric'}, ...
    {'real','finite','integer','scalar','>=',2},mfilename, ...
    'numberOfLineSources');
if size(combinationMasks,2) ~= 3 || isempty(combinationMasks) || ...
        any(~ismember(combinationMasks(:),[false,true]))
    error('combinationMasks must be a nonempty logical array with 3 columns.');
end
combinationMasks = logical(combinationMasks);
if size(unique(combinationMasks,'rows'),1) ~= size(combinationMasks,1)
    error('combinationMasks must not contain duplicate rows.');
end

assert(semiAxes(1) > semiAxes(2) && semiAxes(2) == semiAxes(3), ...
    'The semiaxes must define an x-aligned prolate spheroid [a,b,b].');
focalHalfLength = sqrt(semiAxes(1)^2-semiAxes(2)^2);
assert(islogical(useTooShortSourceLine) && isscalar(useTooShortSourceLine), ...
    'useTooShortSourceLine must be a logical scalar.');
assert(tooShortSourceLineFraction > 0 && tooShortSourceLineFraction < 1, ...
    'tooShortSourceLineFraction must lie strictly between zero and one.');
sourceLineFractions = [1,tooShortSourceLineFraction];
sourceLineModes = ["full focal line","deliberately too short"];
sourceLineChoice = 1+double(useTooShortSourceLine);
sourceLineHalfLength = ...
    sourceLineFractions(sourceLineChoice)*focalHalfLength;
sourceLineMode = sourceLineModes(sourceLineChoice);

if size(stressletNormalDirections,2) ~= 3 || ...
        isempty(stressletNormalDirections) || ...
        any(~isfinite(stressletNormalDirections),'all')
    error('stressletNormalDirections must be a finite K-by-3 array.');
end
normalLengths = vecnorm(stressletNormalDirections,2,2);
if any(normalLengths == 0)
    error('Every stresslet normal direction must be nonzero.');
end
stressletNormalDirections = stressletNormalDirections./normalLengths;
lineDirection = [1,0,0];
if any(abs(stressletNormalDirections*lineDirection.') > 100*eps)
    error('Every stresslet direction must be normal to the focal line.');
end

uniformParameter = linspace(0,1,numberOfLineSources).';
normalizedLineCoordinates = sourceCoordinateMap(uniformParameter);
normalizedLineCoordinates = normalizedLineCoordinates(:);
assert(numel(normalizedLineCoordinates) == numberOfLineSources && ...
    all(isfinite(normalizedLineCoordinates)), ...
    'sourceCoordinateMap must return one finite coordinate per parameter.');
assert(all(diff(normalizedLineCoordinates) > 0), ...
    'sourceCoordinateMap must return strictly increasing coordinates.');
assert(all(abs(normalizedLineCoordinates) <= 1+100*eps), ...
    'Mapped source coordinates must lie in [-1,1].');
assert(max(abs(normalizedLineCoordinates([1,end])-[-1;1])) < 100*eps, ...
    'Both source-line endpoints must be included.');
sourcePoints = [sourceLineHalfLength*normalizedLineCoordinates, ...
    zeros(numberOfLineSources,2)];
assert(~useTooShortSourceLine || ...
    max(abs(sourcePoints(:,1))) < focalHalfLength, ...
    'The shortened source line must not reach either focal point.');

centre = [0,0,0];
surface = ellipsoid_param(semiAxes(1),semiAxes(2),semiAxes(3));
collocationSurface = setupsurfquad(surface,collocationResolution);
testSurface = setupsurfquad(surface,testResolution);
collocationPoints = collocationSurface.x.';
testPoints = testSurface.x.';
collocationVelocity = getKmat(collocationPoints,centre)*rigidMotion;
testVelocity = getKmat(testPoints,centre)*rigidMotion;
testVelocityVectors = reshape(testVelocity,3,[]).';
velocityScale = max(vecnorm(testVelocityVectors,2,2));
assert(velocityScale > 0 && isfinite(velocityScale));

assemblyStart = tic;
collocationStokeslets = stokes_SLP_mat(sourcePoints,collocationPoints);
testStokeslets = stokes_SLP_mat(sourcePoints,testPoints);
[collocationRotlets,collocationNormalStresslets, ...
    collocationPotentialDipoles] = stokes_singularity_matrices( ...
    sourcePoints,collocationPoints,stressletNormalDirections, ...
    logical([1,1,1]));
[testRotlets,testNormalStresslets,testPotentialDipoles] = ...
    stokes_singularity_matrices( ...
    sourcePoints,testPoints,stressletNormalDirections,logical([1,1,1]));
assemblyTime = toc(assemblyStart);

% Scaling gives every fitted block nominal force units: R and T strengths
% have force-times-length units, while D strengths have force-times-length^2.
familyNames = ["S","R","T_perp","D"];
collocationFamilyMatrices = {collocationStokeslets, ...
    singularityLengthScale*collocationRotlets, ...
    singularityLengthScale*collocationNormalStresslets, ...
    singularityLengthScale^2*collocationPotentialDipoles};
testFamilyMatrices = {testStokeslets, ...
    singularityLengthScale*testRotlets, ...
    singularityLengthScale*testNormalStresslets, ...
    singularityLengthScale^2*testPotentialDipoles};
familyColumnCounts = cellfun(@(matrix) size(matrix,2), ...
    collocationFamilyMatrices);
assert(isequal(familyColumnCounts,3*numberOfLineSources* ...
    [1,1,size(stressletNormalDirections,1),1]));

numberOfCombinations = size(combinationMasks,1);
combinationNames = strings(numberOfCombinations,1);
coefficientVectors = cell(numberOfCombinations,1);
coefficientBlocks = cell(numberOfCombinations,numel(familyNames));
unknownCounts = zeros(numberOfCombinations,1);
retainedRanks = zeros(numberOfCombinations,1);
keptConditionNumbers = zeros(numberOfCombinations,1);
coefficientNorms = zeros(numberOfCombinations,1);
familyCoefficientNorms = zeros(numberOfCombinations,numel(familyNames));
collocationMaxRelativeResidual = zeros(numberOfCombinations,1);
testMaxRelativeResidual = zeros(numberOfCombinations,1);
forceAndTorque = zeros(6,numberOfCombinations);
solveTimes = zeros(numberOfCombinations,1);

fprintf('\nSingle-spheroid Stokes line-singularity combination sweep\n');
fprintf('  semiaxes                    [%.2f, %.2f, %.2f]\n',semiAxes);
fprintf('  translation velocity        [%g, %g, %g]^T\n', ...
    translationVelocity);
fprintf('  angular velocity            [%g, %g, %g]^T\n',angularVelocity);
fprintf('  source nodes                %d (%s)\n', ...
    numberOfLineSources,sourceDistributionName);
fprintf('  source-line mode            %s\n',sourceLineMode);
fprintf('  source line                 [%.6f, %.6f]\n', ...
    -sourceLineHalfLength,sourceLineHalfLength);
fprintf('  true focal line             [%.6f, %.6f]\n', ...
    -focalHalfLength,focalHalfLength);
fprintf('  transverse stresslet dirs   %d\n', ...
    size(stressletNormalDirections,1));
fprintf('  singularity length scale    %.3f\n',singularityLengthScale);
fprintf('  collocation/test nodes      %d / %d\n', ...
    size(collocationPoints,1),size(testPoints,1));
fprintf('  SVD cutoff                  %.1e\n',svdCutoff);
fprintf('  combinations                %d\n\n',numberOfCombinations);

for combinationIndex = 1:numberOfCombinations
    activeFamilies = [true,combinationMasks(combinationIndex,:)];
    combinationName = "S";
    for familyIndex = 2:numel(familyNames)
        if activeFamilies(familyIndex)
            combinationName = combinationName+"+"+familyNames(familyIndex);
        end
    end
    combinationNames(combinationIndex) = combinationName;

    systemMatrix = horzcat(collocationFamilyMatrices{activeFamilies});
    testMatrix = horzcat(testFamilyMatrices{activeFamilies});
    unknownCounts(combinationIndex) = size(systemMatrix,2);
    assert(isequal(size(testMatrix), ...
        [3*size(testPoints,1),unknownCounts(combinationIndex)]));
    assert(size(systemMatrix,1) > size(systemMatrix,2), ...
        'Every combination must give an overdetermined system.');

    solveStart = tic;
    [coefficients,retainedRank,keptConditionNumber] = ...
        truncated_svd_solve(systemMatrix,collocationVelocity,svdCutoff);
    solveTimes(combinationIndex) = toc(solveStart);
    coefficientVectors{combinationIndex} = coefficients;
    retainedRanks(combinationIndex) = retainedRank;
    keptConditionNumbers(combinationIndex) = keptConditionNumber;
    coefficientNorms(combinationIndex) = norm(coefficients);

    coefficientOffset = 0;
    for familyIndex = 1:numel(familyNames)
        if activeFamilies(familyIndex)
            blockLength = familyColumnCounts(familyIndex);
            blockIndices = coefficientOffset+(1:blockLength);
            coefficientBlock = coefficients(blockIndices);
            coefficientBlocks{combinationIndex,familyIndex} = coefficientBlock;
            familyCoefficientNorms(combinationIndex,familyIndex) = ...
                norm(coefficientBlock);
            coefficientOffset = coefficientOffset+blockLength;
        else
            coefficientBlocks{combinationIndex,familyIndex} = [];
        end
    end
    assert(coefficientOffset == numel(coefficients));

    collocationError = systemMatrix*coefficients-collocationVelocity;
    testError = testMatrix*coefficients-testVelocity;
    collocationErrorVectors = reshape(collocationError,3,[]).';
    testErrorVectors = reshape(testError,3,[]).';
    collocationMaxRelativeResidual(combinationIndex) = ...
        max(vecnorm(collocationErrorVectors,2,2))/velocityScale;
    testMaxRelativeResidual(combinationIndex) = ...
        max(vecnorm(testErrorVectors,2,2))/velocityScale;

    stokesletCoefficients = coefficientBlocks{combinationIndex,1};
    sourceForceAndTorque = getKmat(sourcePoints,centre).'* ...
        stokesletCoefficients;
    rotletTorque = zeros(3,1);
    if activeFamilies(2)
        scaledRotletCoefficients = coefficientBlocks{combinationIndex,2};
        physicalRotletStrengths = singularityLengthScale*reshape( ...
            scaledRotletCoefficients,3,numberOfLineSources).';
        rotletTorque = sum(physicalRotletStrengths,1).';
    end
    forceAndTorque(:,combinationIndex) = ...
        sourceForceAndTorque+[zeros(3,1);rotletTorque];

    fprintf('%-20s unknowns %4d, rank %4d\n', ...
        combinationName,unknownCounts(combinationIndex),retainedRank);
    fprintf('  coefficient norm                    %.6e\n', ...
        coefficientNorms(combinationIndex));
    fprintf('  collocation/test max rel residual   %.3e / %.3e\n', ...
        collocationMaxRelativeResidual(combinationIndex), ...
        testMaxRelativeResidual(combinationIndex));
end

baselineIndex = find(~any(combinationMasks,2),1);
coefficientNormRatio = nan(numberOfCombinations,1);
if ~isempty(baselineIndex)
    coefficientNormRatio = coefficientNorms/coefficientNorms(baselineIndex);
end

summaryTable = table(combinationNames,combinationMasks(:,1), ...
    combinationMasks(:,2),combinationMasks(:,3),unknownCounts, ...
    retainedRanks,keptConditionNumbers,coefficientNorms, ...
    coefficientNormRatio,collocationMaxRelativeResidual, ...
    testMaxRelativeResidual,'VariableNames', ...
    {'combination','rotlets','normalStresslets','potentialDipoles', ...
    'unknowns','rank','keptConditionNumber','coefficientNorm', ...
    'normRelativeToStokeslets','collocationResidual','testResidual'});

disp(summaryTable);
[~,coefficientNormRanking] = sort(coefficientNorms);
fprintf('\nRanking by scaled coefficient 2-norm\n');
for rankingIndex = 1:numberOfCombinations
    combinationIndex = coefficientNormRanking(rankingIndex);
    fprintf('  %2d. %-20s %.6e  (test residual %.3e)\n', ...
        rankingIndex,combinationNames(combinationIndex), ...
        coefficientNorms(combinationIndex), ...
        testMaxRelativeResidual(combinationIndex));
end

x = 1:numberOfCombinations;
figure('Name','Stokes line-singularity combination sweep');
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
bar(ax,x,coefficientNorms,'FaceColor',[0.20,0.50,0.80]);
set(ax,'YScale','log');
ylabel(ax,'scaled coefficient 2-norm');
title(ax,'Total coefficient norm');
grid(ax,'on');

ax = nexttile(layout);
if isempty(baselineIndex)
    bar(ax,x,nan(size(x)));
    text(ax,0.5,0.5,'Include the S-only row for baseline ratios', ...
        'Units','normalized','HorizontalAlignment','center');
else
    bar(ax,x,coefficientNormRatio,'FaceColor',[0.35,0.70,0.35]);
    yline(ax,1,'k--','S-only baseline');
end
ylabel(ax,'coefficient norm / S-only norm');
title(ax,'Change from Stokeslets');
grid(ax,'on');

ax = nexttile(layout);
semilogy(ax,x,collocationMaxRelativeResidual,'-o','LineWidth',1.4);
hold(ax,'on');
semilogy(ax,x,testMaxRelativeResidual,'-s','LineWidth',1.4);
ylabel(ax,'maximum relative residual');
title(ax,'Boundary accuracy');
legend(ax,{'collocation','independent test'},'Location','best');
grid(ax,'on');

ax = nexttile(layout);
plottedFamilyNorms = familyCoefficientNorms;
plottedFamilyNorms(plottedFamilyNorms == 0) = NaN;
familyHandles = gobjects(numel(familyNames),1);
for familyIndex = 1:numel(familyNames)
    familyHandles(familyIndex) = semilogy(ax,x, ...
        plottedFamilyNorms(:,familyIndex),'-o','LineWidth',1.3);
    hold(ax,'on');
end
ylabel(ax,'coefficient-block 2-norm');
title(ax,'Norm by singularity family');
legend(ax,familyHandles,familyNames,'Location','best');
grid(ax,'on');

axesInLayout = findall(layout,'Type','axes');
for axisIndex = 1:numel(axesInLayout)
    ax = axesInLayout(axisIndex);
    ax.XTick = x;
    ax.XTickLabel = combinationNames;
    ax.XTickLabelRotation = 35;
end
title(layout,sprintf( ...
    'N=%d line nodes, U=[%g,%g,%g]^T, \\Omega=[%g,%g,%g]^T', ...
    numberOfLineSources,translationVelocity,angularVelocity));

crossSectionAngle = linspace(0,2*pi,500);
spheroidCrossSection = [semiAxes(1)*cos(crossSectionAngle); ...
    semiAxes(2)*sin(crossSectionAngle)];
figure('Name','Stokes combination-sweep source line');
nodePlot = scatter(sourcePoints(:,1),sourcePoints(:,2),40, ...
    normalizedLineCoordinates,'filled');
hold on;
sourceLinePlot = plot( ...
    [-sourceLineHalfLength,sourceLineHalfLength],[0,0],'-', ...
    'Color',[0.80,0.20,0.10],'LineWidth',1.5);
focalPlot = plot([-focalHalfLength,focalHalfLength],[0,0], ...
    'kx--','LineWidth',1.0,'MarkerSize',7);
bodyPlot = plot(spheroidCrossSection(1,:),spheroidCrossSection(2,:), ...
    'k--','LineWidth',1.0);
axis equal;
grid on;
xlabel('x');
ylabel('y');
colorbar;
legend([nodePlot,sourceLinePlot,focalPlot,bodyPlot], ...
    {'source node','source line','true focal line and foci', ...
    'spheroid cross-section'},'Location','best');
title(sprintf('%s nodes on the %s', ...
    sourceDistributionName,sourceLineMode));

result.semiAxes = semiAxes;
result.centre = centre;
result.translationVelocity = translationVelocity;
result.angularVelocity = angularVelocity;
result.rigidMotion = rigidMotion;
result.focalHalfLength = focalHalfLength;
result.useTooShortSourceLine = useTooShortSourceLine;
result.tooShortSourceLineFraction = tooShortSourceLineFraction;
result.sourceLineHalfLength = sourceLineHalfLength;
result.sourceLineMode = sourceLineMode;
result.numberOfLineSources = numberOfLineSources;
result.sourcePoints = sourcePoints;
result.normalizedLineCoordinates = normalizedLineCoordinates;
result.sourceDistributionName = sourceDistributionName;
result.sourceCoordinateMap = sourceCoordinateMap;
result.stressletNormalDirections = stressletNormalDirections;
result.singularityLengthScale = singularityLengthScale;
result.familyNames = familyNames;
result.familyColumnCounts = familyColumnCounts;
result.combinationMasks = combinationMasks;
result.combinationNames = combinationNames;
result.summaryTable = summaryTable;
result.coefficientNormRanking = coefficientNormRanking;
result.coefficientVectors = coefficientVectors;
result.coefficientBlocks = coefficientBlocks;
result.coefficientNorms = coefficientNorms;
result.familyCoefficientNorms = familyCoefficientNorms;
result.coefficientNormRatio = coefficientNormRatio;
result.unknownCounts = unknownCounts;
result.retainedRanks = retainedRanks;
result.keptConditionNumbers = keptConditionNumbers;
result.collocationMaxRelativeResidual = collocationMaxRelativeResidual;
result.testMaxRelativeResidual = testMaxRelativeResidual;
result.forceAndTorque = forceAndTorque;
result.collocationResolution = collocationResolution;
result.testResolution = testResolution;
result.velocityScale = velocityScale;
result.svdCutoff = svdCutoff;
result.assemblyTime = assemblyTime;
result.solveTimes = solveTimes;
end

function [coefficients,retainedRank,keptConditionNumber] = ...
    truncated_svd_solve(A,rhs,svdCutoff)
%TRUNCATED_SVD_SOLVE Apply the sweep's common relative SVD cutoff.
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
