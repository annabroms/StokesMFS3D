function result = sep11_single_spheroid_laplace_confocal_line_multipole_experiment
%SEP11_SINGLE_SPHEROID_LAPLACE_CONFOCAL_LINE_MULTIPOLE_EXPERIMENT Compare bases.
%
% Solve a Dirichlet problem on one prolate spheroid using sources only on
% its focal line.  The same source nodes are used for two MFS bases:
%
%   1. one Laplace monopole per focal-line node;
%   2. one monopole, three Cartesian dipoles, and optionally five
%      independent trace-free Cartesian quadrupoles per node.
%
% Both systems use the same fixed collocation nodes on the true spheroid
% and are validated on a separate, finer true-surface grid.  The dipole and
% quadrupole matrices are multiplied by length and length-squared
% scales, respectively, so all fitted coefficients have charge units and
% their Euclidean norms can be compared.  Physical dipole moments and
% quadrupole coefficients are obtained by applying those scales afterward.

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory);
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
addpath(fullfile(projectDirectory,'laplace'));

% User-adjustable experiment parameters.
useTooShortSourceLine = false;
tooShortSourceLineFraction = 0.7;
% This adds five independent trace-free modes per node.  Sources confined
% to the x-axis then span azimuthal orders only through |m|=2.
useQuadrupoles = true;

% Choose "constant", "linear", "quadratic", "mixed-multipoles", or
% "custom".  The first four are configured below by
% make_boundary_data_function.  Edit this handle when using "custom".
boundaryDataChoice = "quadratic";
customBoundaryDataName = "custom harmonic polynomial";
customBoundaryDataFunction = @(points) 1+0.25*points(:,1) ...
    -0.20*points(:,2)+0.15*points(:,1).*points(:,3);
% Set this only when the custom exterior solution's total charge is known.
customExactTotalCharge = NaN;

semiAxes = [1.2,0.7,0.7];
lineSourceCounts = [8,12,16,20,24,32,40,48,64,80];
collocationResolution = [80,60];
testResolution = [83,55];
svdCutoff = 1e-12;  % Change this value to alter every SVD truncation.
dipoleLengthScale = semiAxes(2);
quadrupoleLengthScale = semiAxes(2);
combinedMethodNames = ["monopoles + dipoles", ...
    "monopoles + dipoles + quadrupoles"];
methodNames = ["monopoles",combinedMethodNames(1+double(useQuadrupoles))];

% Edit this map to change the line-node distribution.  It maps uniformly
% spaced parameters in [0,1] to ordered coordinates on the selected source
% line in [-1,1].  The Chebyshev--Lobatto map includes both source-line
% endpoints and clusters near both.
sourceDistributionName = "Chebyshev-Lobatto";
sourceCoordinateMap = @(parameter) -cos(pi*parameter);

assert(semiAxes(1) > semiAxes(2) && semiAxes(2) == semiAxes(3), ...
    'The semiaxes must define an x-aligned prolate spheroid [a,b,b].');
focalHalfLength = sqrt(semiAxes(1)^2-semiAxes(2)^2);
assert(islogical(useTooShortSourceLine) && isscalar(useTooShortSourceLine), ...
    'useTooShortSourceLine must be a logical scalar.');
assert(islogical(useQuadrupoles) && isscalar(useQuadrupoles), ...
    'useQuadrupoles must be a logical scalar.');
assert(tooShortSourceLineFraction > 0 && tooShortSourceLineFraction < 1, ...
    'tooShortSourceLineFraction must lie strictly between zero and one.');
sourceLineFractions = [1,tooShortSourceLineFraction];
sourceLineModes = ["full focal line","deliberately too short"];
sourceLineChoice = 1+double(useTooShortSourceLine);
sourceLineHalfLength = ...
    sourceLineFractions(sourceLineChoice)*focalHalfLength;
sourceLineMode = sourceLineModes(sourceLineChoice);

% With G=1/(4*pi*r), the exact isolated-spheroid capacity is 8*pi/I.
capacityIntegral = integral(@(s) 1./sqrt( ...
    (semiAxes(1)^2+s).*(semiAxes(2)^2+s).*(semiAxes(3)^2+s)), ...
    0,Inf,'RelTol',1e-13,'AbsTol',1e-14);
exactCapacity = 8*pi/capacityIntegral;
[boundaryDataFunction,boundaryDataName,exactTotalCharge, ...
    boundaryDataParameters] = make_boundary_data_function( ...
    boundaryDataChoice,semiAxes,exactCapacity,customBoundaryDataName, ...
    customBoundaryDataFunction,customExactTotalCharge);
if isfinite(exactTotalCharge) && exactTotalCharge ~= 0
    totalChargeErrorScale = abs(exactTotalCharge);
    totalChargeDiagnosticLabel = 'relative total-charge error';
elseif isfinite(exactTotalCharge)
    totalChargeErrorScale = 1;
    totalChargeDiagnosticLabel = 'absolute total-charge error';
else
    totalChargeErrorScale = NaN;
    totalChargeDiagnosticLabel = '|total monopole strength|';
end

surface = ellipsoid_param(semiAxes(1),semiAxes(2),semiAxes(3));
collocationSurface = setupsurfquad(surface,collocationResolution);
testSurface = setupsurfquad(surface,testResolution);
collocationPoints = collocationSurface.x.';
testPoints = testSurface.x.';
collocationData = boundaryDataFunction(collocationPoints);
testData = boundaryDataFunction(testPoints);
collocationData = collocationData(:);
testData = testData(:);
if numel(collocationData) ~= size(collocationPoints,1) || ...
        numel(testData) ~= size(testPoints,1) || ...
        ~isreal(collocationData) || ~isreal(testData) || ...
        any(~isfinite(collocationData)) || any(~isfinite(testData))
    error('The boundary-data function must return one finite real value per point.');
end
boundaryDataScale = max(abs([collocationData;testData]));
if boundaryDataScale == 0
    error('The selected boundary data is identically zero on both grids.');
end

numberOfResolutions = numel(lineSourceCounts);
numberOfMethods = numel(methodNames);
sourceNodeCounts = lineSourceCounts(:);
maximumCombinedUnknowns = (4+5*double(useQuadrupoles))* ...
    max(sourceNodeCounts);
assert(size(collocationPoints,1) > maximumCombinedUnknowns, ...
    'The fixed collocation system must be overdetermined at every resolution.');

collocationMaxResidual = zeros(numberOfResolutions,numberOfMethods);
testMaxResidual = zeros(numberOfResolutions,numberOfMethods);
collocationMaxRelativeResidual = zeros( ...
    numberOfResolutions,numberOfMethods);
testMaxRelativeResidual = zeros(numberOfResolutions,numberOfMethods);
totalMonopoleStrengths = zeros(numberOfResolutions,numberOfMethods);
totalMonopoleStrengthError = nan(numberOfResolutions,numberOfMethods);
coefficientNorms = zeros(numberOfResolutions,numberOfMethods);
monopoleCoefficientNorms = zeros(numberOfResolutions,numberOfMethods);
scaledDipoleCoefficientNorms = zeros(numberOfResolutions,numberOfMethods);
scaledQuadrupoleCoefficientNorms = zeros( ...
    numberOfResolutions,numberOfMethods);
retainedRanks = zeros(numberOfResolutions,numberOfMethods);
keptConditionNumbers = zeros(numberOfResolutions,numberOfMethods);
assemblyTimes = zeros(numberOfResolutions,1);
solveTimes = zeros(numberOfResolutions,numberOfMethods);

finestSourcePoints = [];
finestNormalizedLineCoordinates = [];
finestMonopoleOnlyCoefficients = [];
finestCombinedMonopoleCoefficients = [];
finestCombinedScaledDipoleCoefficients = [];
finestCombinedScaledQuadrupoleCoefficients = [];
finestTestResiduals = zeros(size(testPoints,1),numberOfMethods);

fprintf('\nSingle-spheroid Laplace confocal-line multipole experiment\n');
fprintf('  semiaxes                    [%.2f, %.2f, %.2f]\n',semiAxes);
fprintf('  focal line                  x in [%.6f, %.6f]\n', ...
    -focalHalfLength,focalHalfLength);
fprintf('  source line                 x in [%.6f, %.6f]\n', ...
    -sourceLineHalfLength,sourceLineHalfLength);
fprintf('  source-line mode            %s\n',sourceLineMode);
fprintf('  line-node distribution      %s\n',sourceDistributionName);
fprintf('  boundary data               %s\n',boundaryDataName);
fprintf('  boundary-data scale         %.6e\n',boundaryDataScale);
fprintf('  collocation/test nodes      %d / %d\n', ...
    size(collocationPoints,1),size(testPoints,1));
fprintf('  dipole column length scale  %.3f\n',dipoleLengthScale);
fprintf('  use quadrupoles             %s\n',mat2str(useQuadrupoles));
fprintf('  quadrupole length scale     %.3f\n',quadrupoleLengthScale);
fprintf('  SVD cutoff                  %.1e\n',svdCutoff);
fprintf('  spheroid capacity           %.12f\n',exactCapacity);
if isfinite(exactTotalCharge)
    fprintf('  exact total charge          %.12f\n\n',exactTotalCharge);
else
    fprintf('  exact total charge          not supplied\n\n');
end

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
    collocationMonopoles = lap3d_charge_potential_matrix( ...
        collocationPoints,sourcePoints);
    collocationDipoles = laplace_dipole_potential_matrix( ...
        collocationPoints,sourcePoints);
    collocationQuadrupoles = laplace_quadrupole_potential_matrix( ...
        collocationPoints,sourcePoints,useQuadrupoles);
    testMonopoles = lap3d_charge_potential_matrix(testPoints,sourcePoints);
    testDipoles = laplace_dipole_potential_matrix(testPoints,sourcePoints);
    testQuadrupoles = laplace_quadrupole_potential_matrix( ...
        testPoints,sourcePoints,useQuadrupoles);
    assemblyTimes(resolutionIndex) = toc(assemblyStart);

    collocationMatrices = {collocationMonopoles, ...
        [collocationMonopoles,dipoleLengthScale*collocationDipoles, ...
        quadrupoleLengthScale^2*collocationQuadrupoles]};
    testMatrices = {testMonopoles, ...
        [testMonopoles,dipoleLengthScale*testDipoles, ...
        quadrupoleLengthScale^2*testQuadrupoles]};

    fprintf('line sources %d\n',numberOfSources);
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
            scaledQuadrupoleCoefficients = zeros(numberOfSources,5);
        else
            % Ordering is [q; d_x; d_y; d_z; Q_1; ...; Q_5].
            scaledDipoleCoefficients = reshape( ...
                coefficients(numberOfSources+1:4*numberOfSources), ...
                numberOfSources,3);
            activeQuadrupoleCoefficients = reshape( ...
                coefficients(4*numberOfSources+1:end),numberOfSources, ...
                5*double(useQuadrupoles));
            scaledQuadrupoleCoefficients = zeros(numberOfSources,5);
            scaledQuadrupoleCoefficients( ...
                :,1:size(activeQuadrupoleCoefficients,2)) = ...
                activeQuadrupoleCoefficients;
        end

        collocationError = systemMatrix*coefficients-collocationData;
        testError = testMatrices{methodIndex}*coefficients-testData;
        totalMonopoleStrength = sum(monopoleCoefficients);

        collocationMaxResidual(resolutionIndex,methodIndex) = ...
            norm(collocationError,inf);
        testMaxResidual(resolutionIndex,methodIndex) = norm(testError,inf);
        collocationMaxRelativeResidual(resolutionIndex,methodIndex) = ...
            norm(collocationError,inf)/boundaryDataScale;
        testMaxRelativeResidual(resolutionIndex,methodIndex) = ...
            norm(testError,inf)/boundaryDataScale;
        totalMonopoleStrengths(resolutionIndex,methodIndex) = ...
            totalMonopoleStrength;
        if isfinite(exactTotalCharge)
            totalMonopoleStrengthError(resolutionIndex,methodIndex) = ...
                abs(totalMonopoleStrength-exactTotalCharge)/ ...
                totalChargeErrorScale;
        end
        coefficientNorms(resolutionIndex,methodIndex) = norm(coefficients);
        monopoleCoefficientNorms(resolutionIndex,methodIndex) = ...
            norm(monopoleCoefficients);
        scaledDipoleCoefficientNorms(resolutionIndex,methodIndex) = ...
            norm(scaledDipoleCoefficients,'fro');
        scaledQuadrupoleCoefficientNorms(resolutionIndex,methodIndex) = ...
            norm(scaledQuadrupoleCoefficients,'fro');
        retainedRanks(resolutionIndex,methodIndex) = retainedRank;
        keptConditionNumbers(resolutionIndex,methodIndex) = ...
            keptConditionNumber;

        fprintf('  %-20s unknowns %4d, rank %4d\n', ...
            methodNames(methodIndex),size(systemMatrix,2),retainedRank);
        fprintf('    collocation/test max relative residual  %.3e / %.3e\n', ...
            collocationMaxRelativeResidual(resolutionIndex,methodIndex), ...
            testMaxRelativeResidual(resolutionIndex,methodIndex));
        fprintf('    total monopole strength                 %.6e\n', ...
            totalMonopoleStrength);
        fprintf('    scaled coefficient 2-norm      %.3e\n', ...
            coefficientNorms(resolutionIndex,methodIndex));

        if resolutionIndex == numberOfResolutions
            finestSourcePoints = sourcePoints;
            finestNormalizedLineCoordinates = normalizedLineCoordinates;
            finestTestResiduals(:,methodIndex) = testError;
            if methodIndex == 1
                finestMonopoleOnlyCoefficients = monopoleCoefficients;
            else
                finestCombinedMonopoleCoefficients = monopoleCoefficients;
                finestCombinedScaledDipoleCoefficients = ...
                    scaledDipoleCoefficients;
                finestCombinedScaledQuadrupoleCoefficients = ...
                    scaledQuadrupoleCoefficients;
            end
        end
    end
end

colours = lines(numberOfMethods);
figure('Name','Laplace confocal-line basis convergence');
layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

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
if isfinite(exactTotalCharge)
    chargeDiagnostic = totalMonopoleStrengthError;
else
    chargeDiagnostic = abs(totalMonopoleStrengths);
end
for methodIndex = 1:numberOfMethods
    semilogy(ax,sourceNodeCounts,chargeDiagnostic(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-line nodes');
ylabel(ax,totalChargeDiagnosticLabel);
title(ax,'Total charge diagnostic');
grid(ax,'on');

ax = nexttile(layout);
curveHandles = gobjects(numberOfMethods,1);
for methodIndex = 1:numberOfMethods
    curveHandles(methodIndex) = semilogy(ax,sourceNodeCounts, ...
        coefficientNorms(:,methodIndex),'-o', ...
        'Color',colours(methodIndex,:),'LineWidth',1.5,'MarkerSize',5);
    hold(ax,'on');
end
xlabel(ax,'number of focal-line nodes');
ylabel(ax,'scaled coefficient 2-norm');
title(ax,'Coefficient growth');
grid(ax,'on');
legend(ax,curveHandles,methodNames,'Location','best');
title(layout,sprintf('%s: focal-line basis comparison',boundaryDataName));

figure('Name','Laplace spheroid boundary data');
scatter3(collocationPoints(:,1),collocationPoints(:,2), ...
    collocationPoints(:,3),14,collocationData,'filled');
colormap(turbo(256));
colourBar = colorbar;
colourBar.Label.String = 'Dirichlet data g';
axis equal;
axis tight;
view(34,24);
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('%s: collocation grid [%d,%d], M=%d', ...
    boundaryDataName,collocationResolution,size(collocationPoints,1)));

% Show linear signed residuals and logarithmic residual magnitudes.
surfaceTriangles = convhull(testPoints(:,1),testPoints(:,2),testPoints(:,3));
figure('Name','Laplace residual on the spheroid surface');
layout = tiledlayout(2,numberOfMethods, ...
    'TileSpacing','compact','Padding','compact');
for methodIndex = 1:numberOfMethods
    surfaceResidual = finestTestResiduals(:,methodIndex)/boundaryDataScale;
    maximumAbsoluteResidual = max(abs(surfaceResidual));
    residualFloor = max(maximumAbsoluteResidual*1e-12,realmin('double'));
    surfaceFields = {surfaceResidual, ...
        log10(max(abs(surfaceResidual),residualFloor))};
    colourBarLabels = {'signed relative residual (u_{MFS}-g)/g_{scale}', ...
        'log_{10}|(u_{MFS}-g)/g_{scale}|'};
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
    ['%s, finest focal line N=%d; ', ...
    'max relative residual %.2e / %.2e'], ...
    boundaryDataName, ...
    sourceNodeCounts(end), ...
    max(abs(finestTestResiduals(:,1)))/boundaryDataScale, ...
    max(abs(finestTestResiduals(:,2)))/boundaryDataScale));

crossSectionAngle = linspace(0,2*pi,500);
spheroidCrossSection = [semiAxes(1)*cos(crossSectionAngle); ...
    semiAxes(2)*sin(crossSectionAngle)];
physicalDipoleMoments = dipoleLengthScale* ...
    finestCombinedScaledDipoleCoefficients;
dipoleMagnitudes = vecnorm(physicalDipoleMoments,2,2);
physicalQuadrupoleCoefficients = quadrupoleLengthScale^2* ...
    finestCombinedScaledQuadrupoleCoefficients;
quadrupoleMagnitudes = vecnorm(physicalQuadrupoleCoefficients,2,2);

figure('Name','Laplace confocal-line source strengths');
numberOfSourcePanels = 3+double(useQuadrupoles);
layout = tiledlayout(2,numberOfSourcePanels, ...
    'TileSpacing','compact','Padding','compact');
sourceMagnitudeSets = {abs(finestMonopoleOnlyCoefficients), ...
    abs(finestCombinedMonopoleCoefficients),dipoleMagnitudes, ...
    quadrupoleMagnitudes};
sourceTitles = {'monopoles-only |q|', ...
    'combined-basis |q|','combined-basis |p|', ...
    'combined-basis |Q|_F'};
sourceColourBarLabels = {'|q|','|q|','|p|','|Q|_F'};
maximumDipoleMagnitude = max(dipoleMagnitudes);

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

        if panelIndex == 3 && maximumDipoleMagnitude > 0
            displayScale = 0.12*semiAxes(1)/maximumDipoleMagnitude;
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
            legend(ax,[nodePlot,sourceLinePlot,focalPlot,bodyPlot], ...
                {'source node','source line','true focal line and foci', ...
                'spheroid cross-section'}, ...
                'Location','best');
        end
    end
end
title(layout,sprintf( ...
    ['Source magnitudes, top: linear, bottom: log_{10}; ', ...
    '%s nodes, N=%d; %s'], ...
    sourceDistributionName,sourceNodeCounts(end),sourceLineMode));

result.semiAxes = semiAxes;
result.focalHalfLength = focalHalfLength;
result.useTooShortSourceLine = useTooShortSourceLine;
result.tooShortSourceLineFraction = tooShortSourceLineFraction;
result.sourceLineHalfLength = sourceLineHalfLength;
result.sourceLineMode = sourceLineMode;
result.lineSourceCounts = lineSourceCounts;
result.sourceNodeCounts = sourceNodeCounts;
result.sourceDistributionName = sourceDistributionName;
result.sourceCoordinateMap = sourceCoordinateMap;
result.collocationResolution = collocationResolution;
result.collocationCount = size(collocationPoints,1);
result.testResolution = testResolution;
result.testCount = size(testPoints,1);
result.boundaryDataChoice = boundaryDataChoice;
result.boundaryDataName = boundaryDataName;
result.boundaryDataFunction = boundaryDataFunction;
result.boundaryDataParameters = boundaryDataParameters;
result.boundaryDataScale = boundaryDataScale;
result.collocationData = collocationData;
result.testData = testData;
result.exactTotalCharge = exactTotalCharge;
result.svdCutoff = svdCutoff;
result.dipoleLengthScale = dipoleLengthScale;
result.useQuadrupoles = useQuadrupoles;
result.quadrupoleLengthScale = quadrupoleLengthScale;
result.quadrupoleBasisNames = ["(xx-yy)/sqrt(2)", ...
    "(xx+yy-2zz)/sqrt(6)","sqrt(2)xy", ...
    "sqrt(2)xz","sqrt(2)yz"];
result.methodNames = methodNames;
result.exactCapacity = exactCapacity;
result.collocationMaxResidual = collocationMaxResidual;
result.testMaxResidual = testMaxResidual;
result.collocationMaxRelativeResidual = collocationMaxRelativeResidual;
result.testMaxRelativeResidual = testMaxRelativeResidual;
result.totalMonopoleStrengths = totalMonopoleStrengths;
result.totalMonopoleStrengthError = totalMonopoleStrengthError;
result.totalChargeDiagnostic = chargeDiagnostic;
result.totalChargeDiagnosticLabel = totalChargeDiagnosticLabel;
result.coefficientNorms = coefficientNorms;
result.monopoleCoefficientNorms = monopoleCoefficientNorms;
result.scaledDipoleCoefficientNorms = scaledDipoleCoefficientNorms;
result.scaledQuadrupoleCoefficientNorms = ...
    scaledQuadrupoleCoefficientNorms;
result.retainedRanks = retainedRanks;
result.keptConditionNumbers = keptConditionNumbers;
result.assemblyTimes = assemblyTimes;
result.solveTimes = solveTimes;
result.finestSourcePoints = finestSourcePoints;
result.finestNormalizedLineCoordinates = finestNormalizedLineCoordinates;
result.finestMonopoleOnlyCoefficients = finestMonopoleOnlyCoefficients;
result.finestCombinedMonopoleCoefficients = ...
    finestCombinedMonopoleCoefficients;
result.finestCombinedScaledDipoleCoefficients = ...
    finestCombinedScaledDipoleCoefficients;
result.finestCombinedDipoleMoments = physicalDipoleMoments;
result.finestCombinedScaledQuadrupoleCoefficients = ...
    finestCombinedScaledQuadrupoleCoefficients;
result.finestCombinedQuadrupoleCoefficients = ...
    physicalQuadrupoleCoefficients;
result.finestTestResiduals = finestTestResiduals;
end

function [dataFunction,dataName,exactTotalCharge,parameters] = ...
    make_boundary_data_function(choice,semiAxes,exactCapacity, ...
    customDataName,customDataFunction,customExactTotalCharge)
%MAKE_BOUNDARY_DATA_FUNCTION Keep selectable Dirichlet data in one place.
choice = lower(strtrim(string(choice)));
if ~isscalar(choice) || ismissing(choice)
    error('boundaryDataChoice must be a nonmissing string scalar.');
end
parameters = struct;

switch choice
    case "constant"
        dataName = "unit constant data";
        dataFunction = @(points) ones(size(points,1),1);
        exactTotalCharge = exactCapacity;

    case "linear"
        % This is a harmonic polynomial.  Its odd part has zero total
        % charge on the centred spheroid, leaving the constant-data charge.
        coefficients = [0.35,-0.25,0.20];
        dataName = "constant plus linear harmonic data";
        dataFunction = @(points) 1+ ...
            coefficients(1)*points(:,1)/semiAxes(1)+ ...
            coefficients(2)*points(:,2)/semiAxes(2)+ ...
            coefficients(3)*points(:,3)/semiAxes(3);
        exactTotalCharge = exactCapacity;
        parameters.linearCoefficients = coefficients;

    case "quadratic"
        % x^2-y^2 and y*z are harmonic, so this exercises even and odd
        % angular structure without introducing surface discontinuities.
        coefficients = [0.30,-0.18,0.14,0.12];
        lengthScale = semiAxes(1);
        dataName = "mixed linear and quadratic harmonic data";
        dataFunction = @(points) 1+ ...
            coefficients(1)*points(:,1)/lengthScale+ ...
            coefficients(2)*points(:,3)/lengthScale+ ...
            coefficients(3)*(points(:,1).^2-points(:,2).^2)/ ...
                lengthScale^2+ ...
            coefficients(4)*points(:,2).*points(:,3)/lengthScale^2;
        exactTotalCharge = NaN;
        parameters.polynomialCoefficients = coefficients;
        parameters.lengthScale = lengthScale;

    case "mixed-multipoles"
        % These singularities lie strictly inside the spheroid but not all
        % on its focal line, creating smooth, challenging boundary data.
        chargePoints = [0.18,0.10,-0.08; -0.25,-0.12,0.14];
        chargeStrengths = [1.00;-0.35];
        dipolePoints = [0.05,-0.14,0.09; -0.10,0.16,-0.05];
        dipoleMoments = [0.25,-0.15,0.20; -0.10,0.18,0.12];
        allInteriorPoints = [chargePoints;dipolePoints];
        ellipsoidalRadiiSquared = sum( ...
            (allInteriorPoints./semiAxes).^2,2);
        assert(all(ellipsoidalRadiiSquared < 1), ...
            'Every boundary-data singularity must lie inside the spheroid.');
        dataName = "off-axis interior charges and dipoles";
        dataFunction = @(points) interior_multipole_boundary_data( ...
            points,chargePoints,chargeStrengths,dipolePoints,dipoleMoments);
        exactTotalCharge = sum(chargeStrengths);
        parameters.chargePoints = chargePoints;
        parameters.chargeStrengths = chargeStrengths;
        parameters.dipolePoints = dipolePoints;
        parameters.dipoleMoments = dipoleMoments;

    case "custom"
        if ~isa(customDataFunction,'function_handle')
            error('customBoundaryDataFunction must be a function handle.');
        end
        customDataName = string(customDataName);
        if ~isscalar(customDataName) || ismissing(customDataName) || ...
                strlength(customDataName) == 0
            error('customBoundaryDataName must be a nonempty string scalar.');
        end
        if ~isnumeric(customExactTotalCharge) || ...
                ~isscalar(customExactTotalCharge) || ...
                ~isreal(customExactTotalCharge) || ...
                ~(isfinite(customExactTotalCharge) || ...
                isnan(customExactTotalCharge))
            error(['customExactTotalCharge must be a finite real scalar ', ...
                'or NaN.']);
        end
        dataName = customDataName;
        dataFunction = customDataFunction;
        exactTotalCharge = customExactTotalCharge;

    otherwise
        error(['Unknown boundaryDataChoice "%s". Choose "constant", ', ...
            '"linear", "quadratic", "mixed-multipoles", or "custom".'], ...
            choice);
end
end

function values = interior_multipole_boundary_data( ...
    points,chargePoints,chargeStrengths,dipolePoints,dipoleMoments)
%INTERIOR_MULTIPOLE_BOUNDARY_DATA Evaluate the exact exterior test field.
chargeMatrix = lap3d_charge_potential_matrix(points,chargePoints);
dipoleMatrix = laplace_dipole_potential_matrix(points,dipolePoints);
values = chargeMatrix*chargeStrengths+dipoleMatrix*dipoleMoments(:);
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

function Q = laplace_quadrupole_potential_matrix( ...
    targets,sources,includeQuadrupoles)
%LAPLACE_QUADRUPOLE_POTENTIAL_MATRIX Build five trace-free Hessian columns.
% Columns are grouped by the orthonormal symmetric trace-free tensors
% (xx-yy)/sqrt(2), (xx+yy-2zz)/sqrt(6), sqrt(2)xy, sqrt(2)xz, sqrt(2)yz.
if ~includeQuadrupoles
    Q = zeros(size(targets,1),0);
    return
end
dx = targets(:,1)-sources(:,1).';
dy = targets(:,2)-sources(:,2).';
dz = targets(:,3)-sources(:,3).';
distanceSquared = dx.^2+dy.^2+dz.^2;
if any(distanceSquared == 0,'all')
    error('laplace_quadrupole_potential_matrix:CoincidentPoints', ...
        'A target coincides with a source.');
end
inverseDistanceFifth = distanceSquared.^(-5/2)/(4*pi);
hxx = (3*dx.^2-distanceSquared).*inverseDistanceFifth;
hyy = (3*dy.^2-distanceSquared).*inverseDistanceFifth;
hzz = (3*dz.^2-distanceSquared).*inverseDistanceFifth;
hxy = 3*dx.*dy.*inverseDistanceFifth;
hxz = 3*dx.*dz.*inverseDistanceFifth;
hyz = 3*dy.*dz.*inverseDistanceFifth;
Q = [(hxx-hyy)/sqrt(2),(hxx+hyy-2*hzz)/sqrt(6), ...
    sqrt(2)*hxy,sqrt(2)*hxz,sqrt(2)*hyz];
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
