% OCT14_2026_SINGLE_DROP_MFS  Dense, force-free MFS experiment for one drop.
% Run this script to solve the representative case first, then the sweep.
% The model is a fixed spherical two-fluid interface with a prescribed
% tangential traction jump. The normal points from the drop into the
% exterior fluid, whose viscosity is one.
% All vectors use [x1;y1;z1;x2;y2;z2;...] nodal ordering.

%% Parameters to change
dropRadius = 1;
innerProxyRadius = 0.4;
outerProxyRadius = 2;
S1 = 1;                         % prescribed tangential traction amplitude
interiorViscosity = 2;          % ratio to the exterior viscosity, which is 1
sourceResolutionList = [50 100 200 300 400 500 600];
collocationToSourceRatio = 1.2;
svdTolerance = 1e-13;           % retain singular values > tolerance*sigma_1
tractionRowWeight = 1;          % >1 prioritizes traction in the least-squares fit
scale_columns = false;          % normalize each row-scaled matrix column before SVD
validationResolution = 1000;   % minimum number of independent surface nodes
representativeResolution = 600;
activeDirection = [0 0 1];
planeGridResolution = 121;      % points along each axis of the flow slice
planeHalfWidthInRadii = 2.5;
planeQuiverStride = 7;

%% Paths, parameter checks, and independent validation design
experimentDirectory = fileparts(mfilename('fullpath'));
dropsDirectory = fileparts(experimentDirectory);
repositoryDirectory = fileparts(dropsDirectory);
addpath(dropsDirectory, fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'geometry','spheres'));

assert(0 < innerProxyRadius && innerProxyRadius < dropRadius && ...
    dropRadius < outerProxyRadius, ...
    'Proxy radii must satisfy 0 < inner < drop < outer.');
assert(isfinite(S1) && isscalar(S1));
assert(isfinite(interiorViscosity) && interiorViscosity > 0);
assert(isfinite(svdTolerance) && svdTolerance > 0 && svdTolerance < 1);
assert(isfinite(tractionRowWeight) && tractionRowWeight > 0);
assert(islogical(scale_columns) && isscalar(scale_columns), ...
    'scale_columns must be true or false.');
assert(isfinite(collocationToSourceRatio) && collocationToSourceRatio > 1);
assert(planeGridResolution >= 3 && planeGridResolution == round(planeGridResolution));
assert(planeHalfWidthInRadii > 1 && isfinite(planeHalfWidthInRadii));
assert(planeQuiverStride >= 1 && planeQuiverStride == round(planeQuiverStride));
assert(all(sourceResolutionList >= 3) && ...
    all(sourceResolutionList == round(sourceResolutionList)) && ...
    numel(unique(sourceResolutionList)) == numel(sourceResolutionList));
assert(validationResolution >= 3 && ...
    validationResolution == round(validationResolution));
assert(norm(activeDirection) > 0 && all(isfinite(activeDirection)));
activeDirection = activeDirection / norm(activeDirection);
representativeIndex = find(sourceResolutionList == ...
    representativeResolution, 1);
assert(~isempty(representativeIndex), ...
    'The representative resolution must be in sourceResolutionList.');

% get_sphdesign returns the largest design count no greater than its input.
% Read its catalogue to choose the first design with at least the desired
% collocation or validation count.
catalogueFile = fullfile(repositoryDirectory, 'geometry', 'spheres', ...
    'spherical_designs', 'filelist.txt');
catalogueID = fopen(catalogueFile, 'r');
assert(catalogueID >= 0, 'Could not read the spherical-design catalogue.');
catalogueData = textscan(catalogueID, '%*[^.].%d');
fclose(catalogueID);
availableCounts = catalogueData{1};
validationIndex = find(availableCounts >= validationResolution, 1);
assert(~isempty(validationIndex), ...
    'No spherical design has the requested validation count.');
[validationUnitPoints, ~] = get_sphdesign(availableCounts(validationIndex));

% Rotate the validation design so its points differ from every collocation
% design, even if their catalogued counts happen to coincide.
azimuth = 0.37;
elevation = 0.23;
rotationZ = [cos(azimuth), -sin(azimuth), 0; ...
             sin(azimuth), cos(azimuth), 0; 0, 0, 1];
rotationY = [cos(elevation), 0, sin(elevation); ...
             0, 1, 0; -sin(elevation), 0, cos(elevation)];
validationUnitPoints = (rotationZ*rotationY*validationUnitPoints.').';
validationPoints = dropRadius * validationUnitPoints;
validationNormals = validationUnitPoints;
validationCount = size(validationPoints,1);
fprintf('Independent validation design: requested >= %d, actual %d points.\n', ...
    validationResolution, validationCount);

% Natural physical scales for the two equation types. For S1=0, only
% absolute errors are meaningful, and unit scales keep the SVD defined.
if S1 == 0
    velocityScale = 1;
    tractionScale = 1;
    normalizedLabel = 'unit-scaled';
    fprintf('S1=0: normalized errors use unit scales, not relative scales.\n');
else
    velocityScale = abs(S1)*dropRadius; % exterior viscosity is one
    tractionScale = abs(S1);
    normalizedLabel = 'relative';
end
fprintf('Residual scales: velocity %.6g, traction %.6g.\n', ...
    velocityScale, tractionScale);
fprintf('Traction row weight in SVD: %.6g.\n', tractionRowWeight);
fprintf('Column scaling before SVD: %d.\n', scale_columns);

% U=-completionScale*K'*alpha. This derived scale gives the completion
% term the same velocity-per-force units as the Stokeslet block.
completionScale = 1/(8*pi*dropRadius);
caseOrder = [representativeIndex, ...
    setdiff(1:numel(sourceResolutionList), representativeIndex, 'stable')];
caseCount = numel(sourceResolutionList);
actualSourceCounts = zeros(caseCount,1);
targetCollocationCounts = zeros(caseCount,1);
actualCollocationCounts = zeros(caseCount,1);
retainedRanks = zeros(caseCount,1);
continuityAbsolute = zeros(caseCount,1);
normalAbsolute = zeros(caseCount,1);
tractionAbsolute = zeros(caseCount,1);
forceSumNorm = zeros(caseCount,1);
dropTranslations = zeros(caseCount,3);
maximumExteriorStrength = zeros(caseCount,1);
maximumInteriorStrength = zeros(caseCount,1);
maximumCoefficientComponent = zeros(caseCount,1);

%% Representative resolution first, followed by the remaining sweep
for caseIndex = caseOrder
    requestedN = sourceResolutionList(caseIndex);
    [sourceDirections, ~] = get_sphdesign(requestedN);
    sourceCount = size(sourceDirections,1);
    innerSources = innerProxyRadius * sourceDirections;
    outerSources = outerProxyRadius * sourceDirections;

    targetM = ceil(collocationToSourceRatio*sourceCount);
    collocationIndex = find(availableCounts >= targetM, 1);
    assert(~isempty(collocationIndex), ...
        'No spherical design has the required collocation count.');
    [collocationUnitPoints, ~] = ...
        get_sphdesign(availableCounts(collocationIndex));
    collocationPoints = dropRadius * collocationUnitPoints;
    collocationNormals = collocationUnitPoints;
    collocationCount = size(collocationPoints,1);
    assert(validationCount > collocationCount, ...
        'Increase validationResolution above the collocation count.');
    assert(max(validationNormals*collocationNormals.', [], 'all') < 1-1e-12, ...
        'The validation and collocation sets share a point.');

    fprintf('\nRequested N=%d, actual N=%d; target M=%d, actual M=%d; M/N=%.4f\n', ...
        requestedN, sourceCount, targetM, collocationCount, ...
        collocationCount/sourceCount);

    % Rows: 3M velocity-continuity, M normal-velocity, 2M tangential-traction.
    % Columns: 3N exterior alpha, 3N interior gamma. Force freedom is exact
    % through lambda=P*alpha; pressure gauges vanish under tangent projection.
    [A, tangentOne, tangentTwo] = drop_recompleted_matrix( ...
        innerSources, outerSources, collocationPoints, collocationNormals, ...
        interiorViscosity, completionScale);
    normalDotDirection = collocationNormals*activeDirection.';
    activeTangentOne = S1 * (tangentOne*activeDirection.') .* ...
        normalDotDirection;
    activeTangentTwo = S1 * (tangentTwo*activeDirection.') .* ...
        normalDotDirection;
    tractionData = reshape([activeTangentOne activeTangentTwo].', [], 1);
    rhs = [zeros(4*collocationCount,1); tractionData];

    % The SVD sees dimensionless equations. tractionRowWeight changes the
    % fitted least-squares balance; reported residuals retain physical units.
    A(1:4*collocationCount,:) = ...
        A(1:4*collocationCount,:) / velocityScale;
    A(4*collocationCount+1:end,:) = ...
        tractionRowWeight*A(4*collocationCount+1:end,:) / tractionScale;
    rhs(4*collocationCount+1:end) = ...
        tractionRowWeight*rhs(4*collocationCount+1:end) / tractionScale;

    % Let D contain the Euclidean column norms of the row-scaled matrix.
    % Solve (A/D)*y=rhs, then recover the original [alpha;gamma] as x=D\y.
    % Scaling changes the singular spectrum and hence the truncated solution.
    if scale_columns
        columnScales = vecnorm(A,2,1).'; %#ok<UNRCH> Top-level switch is editable.
        assert(all(isfinite(columnScales) & columnScales > 0), ...
            'Column scaling requires finite, nonzero matrix columns.');
        fprintf('Column-norm range before scaling: %.3e to %.3e.\n', ...
            min(columnScales), max(columnScales));
        A = A ./ columnScales.';
    else
        columnScales = ones(size(A,2),1);
    end

    [leftVectors, singularMatrix, rightVectors] = svd(A, 'econ');
    singularValues = diag(singularMatrix);
    retained = singularValues > svdTolerance*singularValues(1);
    retainedRank = nnz(retained);
    assert(retainedRank > 0, 'The SVD discarded every singular value.');
    scaledCoefficients = rightVectors(:,retained) * ...
        ((leftVectors(:,retained).'*rhs) ./ singularValues(retained));
    coefficients = scaledCoefficients ./ columnScales;
    fprintf('Retained SVD rank: %d of %d columns (cutoff %.3g * sigma_1).\n', ...
        retainedRank, 6*sourceCount, svdTolerance);

    % These are maximum pointwise collocation residuals divided by the same
    % physical scales as the independent validation residuals below.
    fitResidual = A*scaledCoefficients-rhs;
    fitContinuity = max(vecnorm(reshape( ...
        fitResidual(1:3*collocationCount),3,[]).',2,2));
    fitNormal = max(abs(fitResidual(3*collocationCount+1:4*collocationCount)));
    fitTraction = max(vecnorm(reshape( ...
        fitResidual(4*collocationCount+1:end)/tractionRowWeight,2,[]).',2,2));
    fprintf('Collocation normalized maxima: continuity %.3e, normal %.3e, traction %.3e.\n', ...
        fitContinuity, fitNormal, fitTraction);
    clear A leftVectors singularMatrix rightVectors

    alphaNodes = reshape(coefficients(1:3*sourceCount), 3, sourceCount);
    lambdaNodes = alphaNodes - mean(alphaNodes,2);
    lambda = lambdaNodes(:);
    gamma = coefficients(3*sourceCount+1:end);
    translation = -completionScale*sum(alphaNodes,2);
    exteriorStrengthMagnitudes = vecnorm(lambdaNodes,2,1).';
    interiorStrengthMagnitudes = vecnorm(reshape(gamma,3,sourceCount),2,1).';
    maximumExteriorStrength(caseIndex) = max(exteriorStrengthMagnitudes);
    maximumInteriorStrength(caseIndex) = max(interiorStrengthMagnitudes);
    maximumCoefficientComponent(caseIndex) = max(abs(coefficients));
    fprintf('Largest Stokeslet |strength|: exterior %.3e, interior %.3e; largest solved coefficient component %.3e.\n', ...
        maximumExteriorStrength(caseIndex), ...
        maximumInteriorStrength(caseIndex), ...
        maximumCoefficientComponent(caseIndex));

    % Independent surface evaluation of all three imposed conditions.
    exteriorVelocityBlock = drop_velocity_block( ...
        innerSources, validationPoints, 1);
    interiorVelocityBlock = drop_velocity_block( ...
        outerSources, validationPoints, interiorViscosity);
    exteriorVelocity = reshape(exteriorVelocityBlock*lambda, 3, []).';
    interiorVelocity = reshape(interiorVelocityBlock*gamma, 3, []).';
    localContinuityError = vecnorm( ...
        exteriorVelocity-interiorVelocity, 2, 2);
    localNormalError = abs(sum(exteriorVelocity.*validationNormals,2) - ...
        validationNormals*translation);
    clear exteriorVelocityBlock interiorVelocityBlock

    exteriorTractionBlock = drop_traction_block( ...
        innerSources, validationPoints, validationNormals);
    interiorTractionBlock = drop_traction_block( ...
        outerSources, validationPoints, validationNormals);
    tractionJump = reshape( ...
        exteriorTractionBlock*lambda-interiorTractionBlock*gamma, 3, []).';
    jumpNormal = sum(tractionJump.*validationNormals,2);
    projectedJump = tractionJump - jumpNormal.*validationNormals;
    validationDirectionDotNormal = validationNormals*activeDirection.';
    activeTraction = S1 * validationDirectionDotNormal .* ...
        (activeDirection - ...
         validationDirectionDotNormal.*validationNormals);
    localTractionError = vecnorm(projectedJump-activeTraction, 2, 2);
    clear exteriorTractionBlock interiorTractionBlock

    actualSourceCounts(caseIndex) = sourceCount;
    targetCollocationCounts(caseIndex) = targetM;
    actualCollocationCounts(caseIndex) = collocationCount;
    retainedRanks(caseIndex) = retainedRank;
    continuityAbsolute(caseIndex) = max(localContinuityError);
    normalAbsolute(caseIndex) = max(localNormalError);
    tractionAbsolute(caseIndex) = max(localTractionError);
    forceSumNorm(caseIndex) = norm(sum(lambdaNodes,2));
    dropTranslations(caseIndex,:) = translation.';
    fprintf('Validation max velocity continuity: %.3e (%s %.3e)\n', ...
        continuityAbsolute(caseIndex), normalizedLabel, ...
        continuityAbsolute(caseIndex)/velocityScale);
    fprintf('Validation max normal velocity:    %.3e (%s %.3e)\n', ...
        normalAbsolute(caseIndex), normalizedLabel, ...
        normalAbsolute(caseIndex)/velocityScale);
    fprintf('Validation max tangential traction: %.3e (%s %.3e)\n', ...
        tractionAbsolute(caseIndex), normalizedLabel, ...
        tractionAbsolute(caseIndex)/tractionScale);
    fprintf('Sum of exterior strengths (force-free check): %.3e\n', ...
        forceSumNorm(caseIndex));
    fprintf('Drop translation U = [%.6g %.6g %.6g]\n', translation);

    if caseIndex == representativeIndex
        representativeContinuityError = localContinuityError;
        representativeNormalError = localNormalError;
        representativeTractionError = localTractionError;
        representativeActualN = sourceCount;
        representativeInnerSources = innerSources;
        representativeOuterSources = outerSources;
        representativeLambda = lambda;
        representativeGamma = gamma;
        representativeExteriorStrengthMagnitudes = exteriorStrengthMagnitudes;
        representativeInteriorStrengthMagnitudes = interiorStrengthMagnitudes;
    end
end

%% Numerical summary and figures
summaryTable = table(sourceResolutionList(:), actualSourceCounts, ...
    targetCollocationCounts, actualCollocationCounts, retainedRanks, ...
    continuityAbsolute, continuityAbsolute/velocityScale, ...
    normalAbsolute, normalAbsolute/velocityScale, ...
    tractionAbsolute, tractionAbsolute/tractionScale, forceSumNorm, ...
    maximumExteriorStrength, maximumInteriorStrength, ...
    max(maximumExteriorStrength,maximumInteriorStrength), ...
    maximumCoefficientComponent, ...
    dropTranslations(:,1), dropTranslations(:,2), dropTranslations(:,3), ...
    'VariableNames', {'RequestedN','ActualN','TargetM','ActualM', ...
    'RetainedRank','VelocityContinuityAbs','VelocityContinuityNormalized', ...
    'NormalVelocityAbs','NormalVelocityNormalized', ...
    'TangentialTractionAbs','TangentialTractionNormalized', ...
    'ForceSumNorm','MaxExteriorStrength','MaxInteriorStrength', ...
    'MaxStokesletStrength','MaxSolvedCoefficientComponent', ...
    'Ux','Uy','Uz'});
disp(summaryTable);

figure('Name', 'Single-drop MFS convergence');
semilogy(actualSourceCounts, continuityAbsolute/velocityScale, 'o-', ...
    actualSourceCounts, normalAbsolute/velocityScale, 's-', ...
    actualSourceCounts, tractionAbsolute/tractionScale, '^-', ...
    'LineWidth', 1.5);
axis tight;
grid on;
xlabel('Actual source count N per proxy surface');
ylabel('Maximum pointwise error / stated physical scale');
legend('Velocity continuity', 'Normal velocity', ...
    'Tangential traction', 'Location', 'best');
title(sprintf('Single-drop validation on %d independent points', validationCount));

figure('Name', 'Representative velocity-condition errors');
tiledlayout(1,2);
velocityColourMaximum = max([representativeContinuityError; ...
    representativeNormalError]);
if velocityColourMaximum == 0
    velocityColourMaximum = 1;
end
nexttile;
scatter3(validationPoints(:,1), validationPoints(:,2), ...
    validationPoints(:,3), 22, representativeContinuityError, 'filled');
axis equal; grid on; view(3); colorbar;
%clim([0 velocityColourMaximum]);
xlabel('x'); ylabel('y'); zlabel('z');
title('Velocity continuity error');
nexttile;
scatter3(validationPoints(:,1), validationPoints(:,2), ...
    validationPoints(:,3), 22, representativeNormalError, 'filled');
axis equal; grid on; view(3); colorbar;
%clim([0 velocityColourMaximum]);
xlabel('x'); ylabel('y'); zlabel('z');
title('Normal velocity error');
sgtitle(sprintf('Requested N=%d, actual N=%d; velocity units', ...
    representativeResolution, representativeActualN));

figure('Name', 'Representative tangential-traction error');
scatter3(validationPoints(:,1), validationPoints(:,2), ...
    validationPoints(:,3), 22, representativeTractionError, 'filled');
axis equal; grid on; view(3); colorbar;
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('Tangential traction error; requested N=%d, actual N=%d', ...
    representativeResolution, representativeActualN));

%% Representative velocity field in a plane through the drop centre
% The plane contains the active direction. The full 3D speed sets colour;
% unit quiver arrows show the direction projected into this plane.
planeVertical = activeDirection;
referenceAxis = [0 1 0];
if abs(referenceAxis*planeVertical.') > 0.9
    referenceAxis = [1 0 0];
end
planeHorizontal = cross(referenceAxis, planeVertical);
planeHorizontal = planeHorizontal/norm(planeHorizontal);
planeAxis = linspace(-planeHalfWidthInRadii*dropRadius, ...
    planeHalfWidthInRadii*dropRadius, planeGridResolution);
[horizontalGrid, verticalGrid] = meshgrid(planeAxis, planeAxis);
planePoints = horizontalGrid(:)*planeHorizontal + ...
    verticalGrid(:)*planeVertical;
insideDrop = sum(planePoints.^2,2) <= dropRadius^2;
planeVelocity = zeros(size(planePoints));

% Evaluate only each field in its physical domain. Batching keeps the
% temporary dense velocity block small at the largest source resolution.
planeEvaluationBatchSize = 256;
for domain = 1:2
    if domain == 1
        pointIndices = find(insideDrop);
        proxyPoints = representativeOuterSources;
        strengths = representativeGamma;
        fluidViscosity = interiorViscosity;
    else
        pointIndices = find(~insideDrop);
        proxyPoints = representativeInnerSources;
        strengths = representativeLambda;
        fluidViscosity = 1;
    end
    for firstPoint = 1:planeEvaluationBatchSize:numel(pointIndices)
        selected = pointIndices(firstPoint:min( ...
            firstPoint+planeEvaluationBatchSize-1,numel(pointIndices)));
        velocityBlock = drop_velocity_block( ...
            proxyPoints, planePoints(selected,:), fluidViscosity);
        planeVelocity(selected,:) = reshape( ...
            velocityBlock*strengths, 3, []).';
    end
end
planeSpeed = reshape(vecnorm(planeVelocity,2,2), size(horizontalGrid));
horizontalVelocity = reshape(planeVelocity*planeHorizontal.', ...
    size(horizontalGrid));
verticalVelocity = reshape(planeVelocity*planeVertical.', ...
    size(horizontalGrid));
quiverNodes = 1:planeQuiverStride:planeGridResolution;
quiverHorizontal = horizontalVelocity(quiverNodes,quiverNodes);
quiverVertical = verticalVelocity(quiverNodes,quiverNodes);
quiverSpeed = hypot(quiverHorizontal,quiverVertical);
quiverHorizontal = quiverHorizontal./max(quiverSpeed,eps);
quiverVertical = quiverVertical./max(quiverSpeed,eps);

figure('Name', 'Representative drop flow in an active-axis plane');
imagesc(planeAxis, planeAxis, planeSpeed);
set(gca,'YDir','normal');
axis equal tight; hold on;
quiver(horizontalGrid(quiverNodes,quiverNodes), ...
    verticalGrid(quiverNodes,quiverNodes), ...
    quiverHorizontal, quiverVertical, 0.55, 'k', 'LineWidth', 0.8);
circleAngle = linspace(0,2*pi,300);
plot(dropRadius*cos(circleAngle),dropRadius*sin(circleAngle), ...
    'w-', 'LineWidth', 1.2);
colourBar = colorbar;
colourBar.Label.String = '|u| (full 3D speed)';
xlabel('Plane coordinate perpendicular to active direction');
ylabel('Plane coordinate along active direction');
title(sprintf('Drop flow; requested N=%d, actual N=%d', ...
    representativeResolution, representativeActualN));

%% Representative proxy-source strengths
% lambda is the physical exterior Stokeslet strength. gamma is the interior
% Stokeslet strength; alpha also contains the translation completion and is
% therefore not used for the source-strength colours.
largestStrength = max([representativeExteriorStrengthMagnitudes; ...
    representativeInteriorStrengthMagnitudes]);
if largestStrength == 0
    logFloor = 1e-12;
else
    logFloor = max(largestStrength*1e-12, realmin);
end
exteriorLogStrength = log10(max( ...
    representativeExteriorStrengthMagnitudes,logFloor));
interiorLogStrength = log10(max( ...
    representativeInteriorStrengthMagnitudes,logFloor));
colourLimits = [min([exteriorLogStrength;interiorLogStrength]), ...
    max([exteriorLogStrength;interiorLogStrength])];
if colourLimits(1) == colourLimits(2)
    colourLimits = colourLimits + [-0.5 0.5];
end

figure('Name', 'Representative Stokeslet strengths');
tiledlayout(1,2);
nexttile;
scatter3(representativeInnerSources(:,1), ...
    representativeInnerSources(:,2), representativeInnerSources(:,3), ...
    25, exteriorLogStrength, 'filled');
axis equal; grid on; view(3); %clim(colourLimits);
xlabel('x'); ylabel('y'); zlabel('z');
title('Exterior field: inner-proxy |\lambda_j|');
colourBar = colorbar;
colourBar.Label.String = 'log_{10} |\lambda_j|';
nexttile;
scatter3(representativeOuterSources(:,1), ...
    representativeOuterSources(:,2), representativeOuterSources(:,3), ...
    25, interiorLogStrength, 'filled');
axis equal; grid on; view(3); %clim(colourLimits);
xlabel('x'); ylabel('y'); zlabel('z');
title('Interior field: outer-proxy |\gamma_j|');
colourBar = colorbar;
colourBar.Label.String = 'log_{10} |\gamma_j|';
sgtitle(sprintf('Per-node Stokeslet strength; requested N=%d, actual N=%d', ...
    representativeResolution, representativeActualN));
