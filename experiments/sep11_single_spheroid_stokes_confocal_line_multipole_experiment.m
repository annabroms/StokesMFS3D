function result = sep11_single_spheroid_stokes_confocal_line_multipole_experiment
%SEP11_SINGLE_SPHEROID_STOKES_CONFOCAL_LINE_MULTIPOLE_EXPERIMENT Compare bases.
%
% Solve a resistance problem for one prolate spheroid using singularities
% only on its focal line.  The same source nodes are used for two bases:
%
%   1. one three-component Stokeslet (S) per focal-line node;
%   2. the user-selected combination of Stokeslets (S), rotlets (R),
%      transverse-normal stresslets (T), and potential dipoles (D).
%
% Nodal vector components are interleaved as
% [x1;y1;z1;x2;y2;z2;...].  The selected matrix follows the canonical
% block order [S,R,T,D], independent of the order typed below.  Both
% transverse stresslet normals are included together in T.  R and T
% columns are scaled by singularityLengthScale, and D columns by its
% square, so every fitted block has nominal force units.

% ----------------------- USER SETTINGS -------------------------------
translationVelocity = [1;0;0];
angularVelocity = [0;1;0];
useTooShortSourceLine = false;
tooShortSourceLineFraction = 0.8;

% Choose any nonempty subset of "S", "R", "T", and "D".  This selected
% basis is always compared with a separate Stokeslets-only solve.
selectedSingularityTypes = ["S","R","D"];

% The focal line is the x-axis.  T includes one stresslet family for each
% transverse direction below; keep both rows to span the normal plane.
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
lineSourceCounts = [3,8,12,16,20,24,32,40,48];
collocationResolution = [80,60];
testResolution = [83,55];
svdCutoff = 1e-12;  % Change this value to alter every SVD truncation.
singularityLengthScale = semiAxes(2);
familyNames = ["S","R","T","D"];

selectedSingularityTypes = upper(string(selectedSingularityTypes(:).'));
if isempty(selectedSingularityTypes) || ...
        any(ismissing(selectedSingularityTypes)) || ...
        any(strlength(selectedSingularityTypes) == 0)
    error('selectedSingularityTypes must be a nonempty string list.');
end
if any(~ismember(selectedSingularityTypes,familyNames))
    error('Choose selectedSingularityTypes only from "S", "R", "T", and "D".');
end
if numel(unique(selectedSingularityTypes)) ~= numel(selectedSingularityTypes)
    error('selectedSingularityTypes must not contain duplicates.');
end
selectedFamilies = ismember(familyNames,selectedSingularityTypes);
selectedCombinationName = strjoin(familyNames(selectedFamilies),'+');
methodNames = ["S only","selected: "+selectedCombinationName];
activeFamiliesByMethod = [true,false,false,false; selectedFamilies];

% Edit this map to change the line-node distribution.  It maps uniformly
% spaced parameters in [0,1] to ordered coordinates on the selected source
% line in [-1,1].  The Chebyshev--Lobatto map includes both source-line
% endpoints and clusters near both.
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

numberOfResolutions = numel(lineSourceCounts);
numberOfMethods = numel(methodNames);
sourceNodeCounts = lineSourceCounts(:);
familyColumnsPerNode = [3,3,3*size(stressletNormalDirections,1),3];
unknownsPerNode = activeFamiliesByMethod*familyColumnsPerNode.';
unknownCounts = sourceNodeCounts*unknownsPerNode.';
assert(3*size(collocationPoints,1) > max(unknownCounts,[],'all'), ...
    'The fixed collocation system must be overdetermined at every resolution.');

collocationMaxRelativeResidual = zeros( ...
    numberOfResolutions,numberOfMethods);
testMaxRelativeResidual = zeros(numberOfResolutions,numberOfMethods);
coefficientNorms = zeros(numberOfResolutions,numberOfMethods);
familyCoefficientNorms = zeros( ...
    numberOfResolutions,numberOfMethods,numel(familyNames));
coefficientBlocks = cell( ...
    numberOfResolutions,numberOfMethods,numel(familyNames));
forceAndTorque = zeros(6,numberOfResolutions,numberOfMethods);
forceMagnitudes = zeros(numberOfResolutions,numberOfMethods);
torqueMagnitudes = zeros(numberOfResolutions,numberOfMethods);
retainedRanks = zeros(numberOfResolutions,numberOfMethods);
keptConditionNumbers = zeros(numberOfResolutions,numberOfMethods);
assemblyTimes = zeros(numberOfResolutions,1);
solveTimes = zeros(numberOfResolutions,numberOfMethods);

finestSourcePoints = [];
finestNormalizedLineCoordinates = [];
finestCoefficientBlocks = cell(numberOfMethods,numel(familyNames));
finestTestResidualVectors = zeros( ...
    size(testPoints,1),3,numberOfMethods);

fprintf('\nSingle-spheroid Stokes confocal-line multipole experiment\n');
fprintf('  semiaxes                    [%.2f, %.2f, %.2f]\n',semiAxes);
fprintf('  translation velocity        [%g, %g, %g]^T\n', ...
    translationVelocity);
fprintf('  angular velocity            [%g, %g, %g]^T\n',angularVelocity);
fprintf('  focal line                  x in [%.6f, %.6f]\n', ...
    -focalHalfLength,focalHalfLength);
fprintf('  source line                 x in [%.6f, %.6f]\n', ...
    -sourceLineHalfLength,sourceLineHalfLength);
fprintf('  source-line mode            %s\n',sourceLineMode);
fprintf('  line-node distribution      %s\n',sourceDistributionName);
fprintf('  selected singularities      %s\n',selectedCombinationName);
fprintf('  comparison baseline         S only\n');
fprintf('  transverse stresslet dirs   %d\n', ...
    size(stressletNormalDirections,1));
fprintf('  collocation/test nodes      %d / %d\n', ...
    size(collocationPoints,1),size(testPoints,1));
fprintf('  velocity normalization      %.6e\n',velocityScale);
fprintf('  singularity length scale    %.3f\n',singularityLengthScale);
fprintf('  SVD cutoff                  %.1e\n\n',svdCutoff);

for resolutionIndex = 1:numberOfResolutions
    numberOfSources = sourceNodeCounts(resolutionIndex);
    uniformParameter = linspace(0,1,numberOfSources).';
    normalizedLineCoordinates = sourceCoordinateMap(uniformParameter);
    normalizedLineCoordinates = normalizedLineCoordinates(:);

    assert(numel(normalizedLineCoordinates) == numberOfSources && ...
        all(isfinite(normalizedLineCoordinates)), ...
        'sourceCoordinateMap must return one finite coordinate per parameter.');
    assert(all(diff(normalizedLineCoordinates) > 0), ...
        'sourceCoordinateMap must return strictly increasing coordinates.');
    assert(all(abs(normalizedLineCoordinates) <= 1+100*eps), ...
        'Mapped source coordinates must lie in [-1,1].');
    assert(max(abs(normalizedLineCoordinates([1,end])-[-1;1])) < 100*eps, ...
        'This experiment requires both source-line endpoints to be included.');

    sourcePoints = [sourceLineHalfLength*normalizedLineCoordinates, ...
        zeros(numberOfSources,2)];
    assert(all(sourcePoints(:,2) == 0) && all(sourcePoints(:,3) == 0));
    assert(~useTooShortSourceLine || ...
        max(abs(sourcePoints(:,1))) < focalHalfLength, ...
        'The shortened source line must not reach either focal point.');

    assemblyStart = tic;
    collocationStokeslets = stokes_SLP_mat( ...
        sourcePoints,collocationPoints);
    testStokeslets = stokes_SLP_mat(sourcePoints,testPoints);
    [collocationRotlets,collocationNormalStresslets, ...
        collocationPotentialDipoles] = stokes_singularity_matrices( ...
        sourcePoints,collocationPoints,stressletNormalDirections, ...
        selectedFamilies(2:4));
    [testRotlets,testNormalStresslets,testPotentialDipoles] = ...
        stokes_singularity_matrices( ...
        sourcePoints,testPoints,stressletNormalDirections, ...
        selectedFamilies(2:4));
    assemblyTimes(resolutionIndex) = toc(assemblyStart);

    % This is the canonical block order [S,R,T,D].  Scaling gives each
    % fitted coefficient block nominal force units.
    collocationFamilyMatrices = {collocationStokeslets, ...
        singularityLengthScale*collocationRotlets, ...
        singularityLengthScale*collocationNormalStresslets, ...
        singularityLengthScale^2*collocationPotentialDipoles};
    testFamilyMatrices = {testStokeslets, ...
        singularityLengthScale*testRotlets, ...
        singularityLengthScale*testNormalStresslets, ...
        singularityLengthScale^2*testPotentialDipoles};
    familyColumnCounts = numberOfSources*familyColumnsPerNode;

    collocationMatrices = cell(1,numberOfMethods);
    testMatrices = cell(1,numberOfMethods);
    for methodIndex = 1:numberOfMethods
        activeFamilies = activeFamiliesByMethod(methodIndex,:);
        collocationMatrices{methodIndex} = horzcat( ...
            collocationFamilyMatrices{activeFamilies});
        testMatrices{methodIndex} = horzcat( ...
            testFamilyMatrices{activeFamilies});
    end

    fprintf('line sources %d\n',numberOfSources);
    for methodIndex = 1:numberOfMethods
        systemMatrix = collocationMatrices{methodIndex};
        assert(isequal(size(systemMatrix), ...
            [3*size(collocationPoints,1),unknownCounts( ...
            resolutionIndex,methodIndex)]));
        assert(size(systemMatrix,1) > size(systemMatrix,2));

        solveStart = tic;
        [coefficients,retainedRank,keptConditionNumber] = ...
            truncated_svd_solve(systemMatrix,collocationVelocity,svdCutoff);
        solveTimes(resolutionIndex,methodIndex) = toc(solveStart);

        activeFamilies = activeFamiliesByMethod(methodIndex,:);
        coefficientOffset = 0;
        for familyIndex = 1:numel(familyNames)
            if activeFamilies(familyIndex)
                blockLength = familyColumnCounts(familyIndex);
                blockIndices = coefficientOffset+(1:blockLength);
                coefficientBlock = coefficients(blockIndices);
                coefficientBlocks{resolutionIndex,methodIndex,familyIndex} = ...
                    coefficientBlock;
                familyCoefficientNorms( ...
                    resolutionIndex,methodIndex,familyIndex) = ...
                    norm(coefficientBlock);
                coefficientOffset = coefficientOffset+blockLength;
            else
                coefficientBlocks{resolutionIndex,methodIndex,familyIndex} = [];
            end
        end
        assert(coefficientOffset == numel(coefficients));

        collocationError = systemMatrix*coefficients-collocationVelocity;
        testError = testMatrices{methodIndex}*coefficients-testVelocity;
        collocationErrorVectors = reshape(collocationError,3,[]).';
        testErrorVectors = reshape(testError,3,[]).';
        collocationMaxRelativeResidual(resolutionIndex,methodIndex) = ...
            max(vecnorm(collocationErrorVectors,2,2))/velocityScale;
        testMaxRelativeResidual(resolutionIndex,methodIndex) = ...
            max(vecnorm(testErrorVectors,2,2))/velocityScale;

        sourceForceAndTorque = zeros(6,1);
        if activeFamilies(1)
            stokesletCoefficientVector = ...
                coefficientBlocks{resolutionIndex,methodIndex,1};
            sourceForceAndTorque = getKmat(sourcePoints,centre).'* ...
                stokesletCoefficientVector;
        end
        % R has physical strength L*c_R and contributes its summed strength
        % as torque.  The stresslets and potential dipoles have zero net
        % force and torque, so S and R are the only resultant contributors.
        rotletTorque = zeros(3,1);
        if activeFamilies(2)
            scaledRotletCoefficients = ...
                coefficientBlocks{resolutionIndex,methodIndex,2};
            physicalRotletStrengths = singularityLengthScale*reshape( ...
                scaledRotletCoefficients,3,numberOfSources).';
            rotletTorque = sum(physicalRotletStrengths,1).';
        end
        forceAndTorque(:,resolutionIndex,methodIndex) = ...
            sourceForceAndTorque+[zeros(3,1);rotletTorque];

        coefficientNorms(resolutionIndex,methodIndex) = norm(coefficients);
        forceMagnitudes(resolutionIndex,methodIndex) = norm( ...
            forceAndTorque(1:3,resolutionIndex,methodIndex));
        torqueMagnitudes(resolutionIndex,methodIndex) = norm( ...
            forceAndTorque(4:6,resolutionIndex,methodIndex));
        retainedRanks(resolutionIndex,methodIndex) = retainedRank;
        keptConditionNumbers(resolutionIndex,methodIndex) = ...
            keptConditionNumber;

        currentForceAndTorque = forceAndTorque( ...
            :,resolutionIndex,methodIndex);
        fprintf('  %-22s unknowns %4d, rank %4d\n', ...
            methodNames(methodIndex),size(systemMatrix,2),retainedRank);
        fprintf('    collocation/test max relative residual  %.3e / %.3e\n', ...
            collocationMaxRelativeResidual(resolutionIndex,methodIndex), ...
            testMaxRelativeResidual(resolutionIndex,methodIndex));
        fprintf('    force   [% .6e, % .6e, % .6e]^T\n', ...
            currentForceAndTorque(1:3));
        fprintf('    torque  [% .6e, % .6e, % .6e]^T\n', ...
            currentForceAndTorque(4:6));
        fprintf('    scaled coefficient 2-norm              %.3e\n', ...
            coefficientNorms(resolutionIndex,methodIndex));

        if resolutionIndex == numberOfResolutions
            finestSourcePoints = sourcePoints;
            finestNormalizedLineCoordinates = normalizedLineCoordinates;
            finestTestResidualVectors(:,:,methodIndex) = testErrorVectors;
            for familyIndex = 1:numel(familyNames)
                finestCoefficientBlocks{methodIndex,familyIndex} = ...
                    coefficientBlocks{resolutionIndex,methodIndex,familyIndex};
            end
        end
    end
end
coefficientNormRatioToStokeslets = ...
    coefficientNorms(:,2)./coefficientNorms(:,1);

colours = lines(numberOfMethods);
figure('Name','Stokes confocal-line basis convergence');
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
for methodIndex = 1:numberOfMethods
    semilogy(ax,sourceNodeCounts, ...
        testMaxRelativeResidual(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-line nodes');
ylabel(ax,'independent max relative residual');
title(ax,'Boundary accuracy');
grid(ax,'on');

ax = nexttile(layout);
for methodIndex = 1:numberOfMethods
    semilogy(ax,sourceNodeCounts,coefficientNorms(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-line nodes');
ylabel(ax,'scaled coefficient 2-norm');
title(ax,'Coefficient growth');
grid(ax,'on');

ax = nexttile(layout);
for methodIndex = 1:numberOfMethods
    plot(ax,sourceNodeCounts,forceMagnitudes(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-line nodes');
ylabel(ax,'|F|');
title(ax,'Resultant force');
grid(ax,'on');

ax = nexttile(layout);
curveHandles = gobjects(numberOfMethods,1);
for methodIndex = 1:numberOfMethods
    curveHandles(methodIndex) = plot(ax,sourceNodeCounts, ...
        torqueMagnitudes(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-line nodes');
ylabel(ax,'|T|');
title(ax,'Resultant torque');
grid(ax,'on');
legend(ax,curveHandles,methodNames,'Location','best');
title(layout,sprintf( ...
    'Single spheroid, U=[%g,%g,%g]^T, \\Omega=[%g,%g,%g]^T', ...
    translationVelocity,angularVelocity));

figure('Name','Stokes spheroid collocation grid');
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

% Show linear and logarithmic magnitudes of the vector residual.
surfaceTriangles = convhull(testPoints(:,1),testPoints(:,2),testPoints(:,3));
figure('Name','Stokes residual on the spheroid surface');
layout = tiledlayout(2,numberOfMethods, ...
    'TileSpacing','compact','Padding','compact');
for methodIndex = 1:numberOfMethods
    relativeResidualMagnitude = vecnorm( ...
        finestTestResidualVectors(:,:,methodIndex),2,2)/velocityScale;
    maximumResidual = max(relativeResidualMagnitude);
    residualFloor = max(maximumResidual*1e-12,realmin('double'));
    surfaceFields = {relativeResidualMagnitude, ...
        log10(max(relativeResidualMagnitude,residualFloor))};
    colourBarLabels = {'|u_{MFS}-u_{rigid}|/U_{scale}', ...
        'log_{10}(|u_{MFS}-u_{rigid}|/U_{scale})'};
    scaleNames = {'linear','log_{10} magnitude'};

    for scaleIndex = 1:2
        tileIndex = (scaleIndex-1)*numberOfMethods+methodIndex;
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
title(layout,sprintf( ...
    'Independent velocity residual, finest focal line N=%d; max %.2e / %.2e', ...
    sourceNodeCounts(end), ...
    testMaxRelativeResidual(end,1),testMaxRelativeResidual(end,2)));

crossSectionAngle = linspace(0,2*pi,500);
spheroidCrossSection = [semiAxes(1)*cos(crossSectionAngle); ...
    semiAxes(2)*sin(crossSectionAngle)];

numberOfSourcePanels = 1+nnz(selectedFamilies);
sourceMagnitudeSets = cell(1,numberOfSourcePanels);
sourceTitles = cell(1,numberOfSourcePanels);
sourceColourBarLabels = cell(1,numberOfSourcePanels);
sourceVectorSets = cell(1,numberOfSourcePanels);

baselineStokesletCoefficients = reshape( ...
    finestCoefficientBlocks{1,1},3,sourceNodeCounts(end)).';
sourceMagnitudeSets{1} = vecnorm(baselineStokesletCoefficients,2,2);
sourceTitles{1} = 'S-only coefficients';
sourceColourBarLabels{1} = '|c_S|';
sourceVectorSets{1} = baselineStokesletCoefficients;

panelIndex = 1;
for familyIndex = find(selectedFamilies)
    panelIndex = panelIndex+1;
    coefficientBlock = finestCoefficientBlocks{2,familyIndex};
    directionsPerNode = familyColumnsPerNode(familyIndex)/3;
    sourceMagnitudeSets{panelIndex} = nodal_coefficient_magnitudes( ...
        coefficientBlock,sourceNodeCounts(end),directionsPerNode);
    sourceTitles{panelIndex} = sprintf( ...
        'selected %s coefficients',familyNames(familyIndex));
    sourceColourBarLabels{panelIndex} = sprintf( ...
        '|c_%s|',familyNames(familyIndex));
    if directionsPerNode == 1
        sourceVectorSets{panelIndex} = reshape( ...
            coefficientBlock,3,sourceNodeCounts(end)).';
    else
        sourceVectorSets{panelIndex} = [];
    end
end

figure('Name','Stokes confocal-line source strengths');
layout = tiledlayout(2,numberOfSourcePanels, ...
    'TileSpacing','compact','Padding','compact');

for scaleIndex = 1:2
    for panelIndex = 1:numberOfSourcePanels
        tileIndex = (scaleIndex-1)*numberOfSourcePanels+panelIndex;
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
            finestSourcePoints(:,2),40,sourceField,'filled');
        hold(ax,'on');
        focalPlot = plot(ax,[-focalHalfLength,focalHalfLength],[0,0], ...
            'kx--','LineWidth',1.0,'MarkerSize',7);
        sourceLinePlot = plot(ax, ...
            [-sourceLineHalfLength,sourceLineHalfLength],[0,0],'-', ...
            'Color',[0.80,0.20,0.10],'LineWidth',1.5);
        bodyPlot = plot(ax,spheroidCrossSection(1,:), ...
            spheroidCrossSection(2,:),'k--','LineWidth',1.0);
        fieldLimits = [min(sourceField),max(sourceField)];
        if fieldLimits(1) < fieldLimits(2)
            % Every source panel uses its own observed range.
            clim(ax,fieldLimits);
        end
        colormap(ax,turbo(256));

        sourceVectors = sourceVectorSets{panelIndex};
        if ~isempty(sourceVectors) && maximumSourceMagnitude > 0
            displayScale = 0.12*semiAxes(1)/maximumSourceMagnitude;
            quiver(ax,finestSourcePoints(:,1),finestSourcePoints(:,2), ...
                displayScale*sourceVectors(:,1), ...
                displayScale*sourceVectors(:,2),0,'k', ...
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
            legend(ax,[nodePlot,sourceLinePlot,focalPlot,bodyPlot], ...
                {'source node','source line','true focal line and foci', ...
                'spheroid cross-section'}, ...
                'Location','best');
        end
    end
end
title(layout,sprintf( ...
    ['Scaled coefficient magnitudes, top: linear, bottom: log_{10}; ', ...
    '%s nodes, N=%d; %s'], ...
    sourceDistributionName,sourceNodeCounts(end),sourceLineMode));

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
result.lineSourceCounts = lineSourceCounts;
result.sourceNodeCounts = sourceNodeCounts;
result.unknownCounts = unknownCounts;
result.sourceDistributionName = sourceDistributionName;
result.sourceCoordinateMap = sourceCoordinateMap;
result.collocationResolution = collocationResolution;
result.collocationCount = size(collocationPoints,1);
result.testResolution = testResolution;
result.testCount = size(testPoints,1);
result.velocityScale = velocityScale;
result.svdCutoff = svdCutoff;
result.selectedSingularityTypes = selectedSingularityTypes;
result.selectedCombinationName = selectedCombinationName;
result.familyNames = familyNames;
result.selectedFamilies = selectedFamilies;
result.activeFamiliesByMethod = activeFamiliesByMethod;
result.familyColumnsPerNode = familyColumnsPerNode;
result.stressletNormalDirections = stressletNormalDirections;
result.singularityLengthScale = singularityLengthScale;
result.methodNames = methodNames;
result.collocationMaxRelativeResidual = collocationMaxRelativeResidual;
result.testMaxRelativeResidual = testMaxRelativeResidual;
result.coefficientNorms = coefficientNorms;
result.coefficientNormRatioToStokeslets = ...
    coefficientNormRatioToStokeslets;
result.familyCoefficientNorms = familyCoefficientNorms;
result.coefficientBlocks = coefficientBlocks;
result.stokesletCoefficientNorms = familyCoefficientNorms(:,:,1);
result.forceAndTorque = forceAndTorque;
result.forceMagnitudes = forceMagnitudes;
result.torqueMagnitudes = torqueMagnitudes;
result.retainedRanks = retainedRanks;
result.keptConditionNumbers = keptConditionNumbers;
result.assemblyTimes = assemblyTimes;
result.solveTimes = solveTimes;
result.finestSourcePoints = finestSourcePoints;
result.finestNormalizedLineCoordinates = finestNormalizedLineCoordinates;
result.finestCoefficientBlocks = finestCoefficientBlocks;
result.finestTestResidualVectors = finestTestResidualVectors;
end

function magnitudes = nodal_coefficient_magnitudes( ...
    coefficientBlock,numberOfSources,directionsPerNode)
%NODAL_COEFFICIENT_MAGNITUDES Collapse vector directions at each node.
assert(numel(coefficientBlock) == ...
    3*numberOfSources*directionsPerNode);
squaredMagnitudes = zeros(numberOfSources,1);
for directionIndex = 1:directionsPerNode
    directionOffset = (directionIndex-1)*3*numberOfSources;
    directionIndices = directionOffset+(1:3*numberOfSources);
    nodalVectors = reshape( ...
        coefficientBlock(directionIndices),3,numberOfSources).';
    squaredMagnitudes = squaredMagnitudes+sum(nodalVectors.^2,2);
end
magnitudes = sqrt(squaredMagnitudes);
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
