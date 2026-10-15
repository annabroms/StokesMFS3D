% OCT14_2026_TWO_DROP_GAP_SWEEP Dense SVD experiment for two spherical drops.
% The centres lie on the x axis, separated by 2*dropRadius + delta. Both
% drops have the same active direction. Run this script to sweep delta.
% The shared incident field follows MFS_drops-1.pdf, equations (41)-(43).
% N and M are counts PER DROP; vectors are interleaved by Cartesian node.

%% Parameters to change (matching the current single-drop experiment)
dropRadius = 1;
innerProxyRadius = 0.6;
outerProxyRadius = 2;
S1 = 1;
interiorViscosity = 2;          % ratio to exterior viscosity 1
requestedN = 600;
collocationToSourceRatio = 1.2;
svdTolerance = 1e-13;
tractionRowWeight = 1;
scale_columns = false;
validationResolution = 1000;   % minimum independent points per sphere
gapList = [1 0.5 0.2 0.1 0.05 0.02];
representativeGap = 0.1;
activeDirection = [0 0 1];      % same prescribed direction on both drops
planeGridResolution = [181 121]; % [horizontal, vertical]
planePaddingRadii = 1.5;
planeQuiverStride = 8;

%% Paths, geometry grids, and physical residual scales
experimentDirectory = fileparts(mfilename('fullpath'));
dropsDirectory = fileparts(experimentDirectory);
repositoryDirectory = fileparts(dropsDirectory);
addpath(dropsDirectory,fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'geometry','spheres'));

assert(0 < innerProxyRadius && innerProxyRadius < dropRadius && ...
    dropRadius < outerProxyRadius, ...
    'Proxy radii must satisfy 0 < inner < drop < outer.');
assert(isfinite(S1) && isscalar(S1));
assert(isfinite(interiorViscosity) && interiorViscosity > 0);
assert(requestedN >= 3 && requestedN == round(requestedN));
assert(collocationToSourceRatio > 1 && isfinite(collocationToSourceRatio));
assert(svdTolerance > 0 && svdTolerance < 1 && isfinite(svdTolerance));
assert(tractionRowWeight > 0 && isfinite(tractionRowWeight));
assert(islogical(scale_columns) && isscalar(scale_columns));
assert(validationResolution >= 3 && ...
    validationResolution == round(validationResolution));
assert(isvector(gapList) && all(gapList > 0) && ...
    all(isfinite(gapList)) && numel(unique(gapList)) == numel(gapList));
assert(any(gapList == representativeGap), ...
    'representativeGap must appear in gapList.');
assert(isequal(size(planeGridResolution),[1 2]) && ...
    all(planeGridResolution >= 3) && ...
    all(planeGridResolution == round(planeGridResolution)));
assert(planePaddingRadii > 0 && isfinite(planePaddingRadii));
assert(planeQuiverStride >= 1 && planeQuiverStride == round(planeQuiverStride));
assert(all(isfinite(activeDirection)) && norm(activeDirection) > 0);
activeDirection = activeDirection/norm(activeDirection);

catalogueFile = fullfile(repositoryDirectory,'geometry','spheres', ...
    'spherical_designs','filelist.txt');
catalogueID = fopen(catalogueFile,'r');
assert(catalogueID >= 0,'Could not read the spherical-design catalogue.');
catalogueData = textscan(catalogueID,'%*[^.].%d');
fclose(catalogueID);
availableCounts = catalogueData{1};

[sourceDirections,~] = get_sphdesign(requestedN);
sourceCount = size(sourceDirections,1);
targetM = ceil(collocationToSourceRatio*sourceCount);
collocationIndex = find(availableCounts >= targetM,1);
validationIndex = find(availableCounts >= validationResolution,1);
assert(~isempty(collocationIndex) && ~isempty(validationIndex), ...
    'The spherical-design catalogue cannot supply the requested counts.');
[collocationDirections,~] = get_sphdesign(availableCounts(collocationIndex));
[validationDirections,~] = get_sphdesign(availableCounts(validationIndex));
collocationCount = size(collocationDirections,1);
validationCount = size(validationDirections,1);
assert(validationCount > collocationCount, ...
    'The independent validation design must have more points than collocation.');

% Rotate the independent validation design away from the collocation nodes.
azimuth = 0.37;
elevation = 0.23;
rotationZ = [cos(azimuth) -sin(azimuth) 0; ...
             sin(azimuth) cos(azimuth) 0;0 0 1];
rotationY = [cos(elevation) 0 sin(elevation); ...
             0 1 0;-sin(elevation) 0 cos(elevation)];
validationDirections = (rotationZ*rotationY*validationDirections.').';
assert(max(validationDirections*collocationDirections.',[],'all') < 1-1e-12, ...
    'Validation and collocation designs share a point.');

if S1 == 0
    velocityScale = 1;
    tractionScale = 1;
    normalizedLabel = 'unit-scaled';
else
    velocityScale = abs(S1)*dropRadius; % exterior viscosity is one
    tractionScale = abs(S1);
    normalizedLabel = 'relative';
end
completionScale = 1/(8*pi*dropRadius);
fprintf('Per sphere: requested N=%d, actual N=%d; target M=%d, actual M=%d (M/N=%.4f).\n', ...
    requestedN,sourceCount,targetM,collocationCount, ...
    collocationCount/sourceCount);
fprintf('Independent validation: requested >=%d, actual %d points per sphere.\n', ...
    validationResolution,validationCount);
fprintf('Residual scales: velocity %.6g, traction %.6g; traction row weight %.6g; column scaling %d.\n', ...
    velocityScale,tractionScale,tractionRowWeight,scale_columns);
fprintf('Dense system per gap: %d rows x %d columns (%.2f GiB for A alone).\n', ...
    12*collocationCount,12*sourceCount, ...
    8*(12*collocationCount)*(12*sourceCount)/2^30);

caseCount = numel(gapList);
retainedRanks = zeros(caseCount,1);
continuityAbsolute = zeros(caseCount,2);
normalAbsolute = zeros(caseCount,2);
tractionAbsolute = zeros(caseCount,2);
forceSumNorm = zeros(caseCount,2);
maximumExteriorStrength = zeros(caseCount,2);
maximumInteriorStrength = zeros(caseCount,2);
maximumCoefficientComponent = zeros(caseCount,2);
translations = zeros(caseCount,3,2);
representativeIndex = find(gapList == representativeGap,1);
caseOrder = [representativeIndex, ...
    setdiff(1:caseCount,representativeIndex,'stable')];

%% Representative gap first, then the remaining gap sweep
for caseIndex = caseOrder
    delta = gapList(caseIndex);
    centreDistance = 2*dropRadius+delta;
    centres = [-centreDistance/2 0 0;centreDistance/2 0 0];
    innerSources = cell(1,2);
    outerSources = cell(1,2);
    collocationPoints = cell(1,2);
    validationPoints = cell(1,2);
    normals = cell(1,2);
    for body = 1:2
        innerSources{body} = centres(body,:) + ...
            innerProxyRadius*sourceDirections;
        outerSources{body} = centres(body,:) + ...
            outerProxyRadius*sourceDirections;
        collocationPoints{body} = centres(body,:) + ...
            dropRadius*collocationDirections;
        validationPoints{body} = centres(body,:) + ...
            dropRadius*validationDirections;
        normals{body} = collocationDirections;
    end
    fprintf('\nGap delta=%.6g (centre distance %.6g)\n',delta,centreDistance);

    % Body-block order [self 1, interaction 1<-2; interaction 2<-1, self 2].
    % Each body's columns are [alpha;gamma], and each body's rows are
    % [velocity continuity; normal velocity; two tangent traction values].
    [A,tangentOne,tangentTwo] = drop_two_body_matrix( ...
        innerSources,outerSources,collocationPoints,normals, ...
        interiorViscosity,completionScale);
    rhs = zeros(12*collocationCount,1);
    for body = 1:2
        normalDotDirection = normals{body}*activeDirection.';
        tangentData = S1 * [ ...
            (tangentOne{body}*activeDirection.').*normalDotDirection, ...
            (tangentTwo{body}*activeDirection.').*normalDotDirection];
        bodyOffset = (body-1)*6*collocationCount;
        rhs(bodyOffset+4*collocationCount+(1:2*collocationCount)) = ...
            reshape(tangentData.',[],1);

        velocityRows = bodyOffset+(1:4*collocationCount);
        tractionRows = bodyOffset+4*collocationCount+(1:2*collocationCount);
        A(velocityRows,:) = A(velocityRows,:)/velocityScale;
        A(tractionRows,:) = tractionRowWeight*A(tractionRows,:)/tractionScale;
        rhs(tractionRows) = tractionRowWeight*rhs(tractionRows)/tractionScale;
    end

    % Optional right scaling uses norms AFTER physical row scaling. The SVD
    % solves (A/D)*y=rhs, and original [alpha_1;gamma_1;alpha_2;gamma_2]
    % are recovered by x=D\y. Truncation acts on this chosen matrix.
    if scale_columns
        columnScales = vecnorm(A,2,1).'; %#ok<UNRCH> Top-level switch is editable.
        assert(all(isfinite(columnScales) & columnScales > 0), ...
            'Column scaling requires finite, nonzero matrix columns.');
        fprintf('Column-norm range before scaling: %.3e to %.3e.\n', ...
            min(columnScales),max(columnScales));
        A = A./columnScales.';
    else
        columnScales = ones(size(A,2),1);
    end

    [leftVectors,singularMatrix,rightVectors] = svd(A,'econ');
    singularValues = diag(singularMatrix);
    retained = singularValues > svdTolerance*singularValues(1);
    retainedRanks(caseIndex) = nnz(retained);
    assert(retainedRanks(caseIndex) > 0, ...
        'The SVD discarded every singular value.');
    scaledCoefficients = rightVectors(:,retained) * ...
        ((leftVectors(:,retained).'*rhs)./singularValues(retained));
    coefficients = scaledCoefficients./columnScales;
    fprintf('Retained SVD rank: %d of %d columns (cutoff %.3g * sigma_1).\n', ...
        retainedRanks(caseIndex),12*sourceCount,svdTolerance);

    fitResidual = A*scaledCoefficients-rhs;
    for body = 1:2
        bodyOffset = (body-1)*6*collocationCount;
        localFit = fitResidual(bodyOffset+(1:6*collocationCount));
        fitContinuity = max(vecnorm(reshape( ...
            localFit(1:3*collocationCount),3,[]).',2,2));
        fitNormal = max(abs(localFit(3*collocationCount+1:4*collocationCount)));
        fitTraction = max(vecnorm(reshape( ...
            localFit(4*collocationCount+1:end)/tractionRowWeight,2,[]).',2,2));
        fprintf('Sphere %d collocation normalized maxima: continuity %.3e, normal %.3e, traction %.3e.\n', ...
            body,fitContinuity,fitNormal,fitTraction);
    end
    clear A leftVectors singularMatrix rightVectors

    lambda = cell(1,2);
    gamma = cell(1,2);
    velocityU = zeros(2,3);
    exteriorMagnitudes = cell(1,2);
    interiorMagnitudes = cell(1,2);
    for body = 1:2
        coefficientOffset = (body-1)*6*sourceCount;
        bodyCoefficients = coefficients(coefficientOffset+(1:6*sourceCount));
        alphaNodes = reshape(bodyCoefficients(1:3*sourceCount),3,sourceCount);
        lambdaNodes = alphaNodes-mean(alphaNodes,2);
        lambda{body} = lambdaNodes(:);
        gamma{body} = bodyCoefficients(3*sourceCount+1:end);
        velocityU(body,:) = (-completionScale*sum(alphaNodes,2)).';
        forceSumNorm(caseIndex,body) = norm(sum(lambdaNodes,2));
        exteriorMagnitudes{body} = vecnorm(lambdaNodes,2,1).';
        interiorMagnitudes{body} = vecnorm( ...
            reshape(gamma{body},3,sourceCount),2,1).';
        maximumExteriorStrength(caseIndex,body) = max(exteriorMagnitudes{body});
        maximumInteriorStrength(caseIndex,body) = max(interiorMagnitudes{body});
        maximumCoefficientComponent(caseIndex,body) = max(abs(bodyCoefficients));
        translations(caseIndex,:,body) = velocityU(body,:);
        fprintf('Sphere %d: |sum lambda|=%.3e, max |lambda_j|=%.3e, max |gamma_j|=%.3e, max coefficient component=%.3e; U=[%.5g %.5g %.5g].\n', ...
            body,forceSumNorm(caseIndex,body), ...
            maximumExteriorStrength(caseIndex,body), ...
            maximumInteriorStrength(caseIndex,body), ...
            maximumCoefficientComponent(caseIndex,body),velocityU(body,:));
    end

    % At each independent surface, exterior velocity includes both lambda
    % fields. Interior velocity contains its own gamma and the other
    % sphere's incident lambda. The latter cancels in continuity but its
    % traction remains with coefficient (1-mu).
    for body = 1:2
        other = 3-body;
        targets = validationPoints{body};
        targetNormals = validationDirections;
        selfVelocityBlock = drop_velocity_block(innerSources{body},targets,1);
        incidentVelocityBlock = drop_velocity_block(innerSources{other},targets,1);
        interiorVelocityBlock = drop_velocity_block( ...
            outerSources{body},targets,interiorViscosity);
        selfVelocity = selfVelocityBlock*lambda{body};
        incidentVelocity = incidentVelocityBlock*lambda{other};
        exteriorVelocity = reshape(selfVelocity+incidentVelocity,3,[]).';
        interiorVelocity = reshape( ...
            interiorVelocityBlock*gamma{body}+incidentVelocity,3,[]).';
        localContinuityError = vecnorm(exteriorVelocity-interiorVelocity,2,2);
        localNormalError = abs(sum(exteriorVelocity.*targetNormals,2) - ...
            targetNormals*velocityU(body,:).');
        clear selfVelocityBlock incidentVelocityBlock interiorVelocityBlock

        selfTractionBlock = drop_traction_block( ...
            innerSources{body},targets,targetNormals);
        incidentTractionBlock = drop_traction_block( ...
            innerSources{other},targets,targetNormals);
        interiorTractionBlock = drop_traction_block( ...
            outerSources{body},targets,targetNormals);
        tractionJump = reshape( ...
            selfTractionBlock*lambda{body} + ...
            (1-interiorViscosity)*incidentTractionBlock*lambda{other} - ...
            interiorTractionBlock*gamma{body},3,[]).';
        projectedJump = tractionJump - ...
            sum(tractionJump.*targetNormals,2).*targetNormals;
        normalDotDirection = targetNormals*activeDirection.';
        activeTraction = S1*normalDotDirection.* ...
            (activeDirection-normalDotDirection.*targetNormals);
        localTractionError = vecnorm(projectedJump-activeTraction,2,2);
        clear selfTractionBlock incidentTractionBlock interiorTractionBlock

        continuityAbsolute(caseIndex,body) = max(localContinuityError);
        normalAbsolute(caseIndex,body) = max(localNormalError);
        tractionAbsolute(caseIndex,body) = max(localTractionError);
        fprintf('Sphere %d validation max: continuity %.3e (%s %.3e), normal %.3e (%s %.3e), traction %.3e (%s %.3e).\n', ...
            body,continuityAbsolute(caseIndex,body),normalizedLabel, ...
            continuityAbsolute(caseIndex,body)/velocityScale, ...
            normalAbsolute(caseIndex,body),normalizedLabel, ...
            normalAbsolute(caseIndex,body)/velocityScale, ...
            tractionAbsolute(caseIndex,body),normalizedLabel, ...
            tractionAbsolute(caseIndex,body)/tractionScale);

        if caseIndex == representativeIndex && body == 1
            representativeContinuityError = localContinuityError;
            representativeNormalError = localNormalError;
            representativeTractionError = localTractionError;
            representativeValidationPoints = targets;
        end
    end
    if caseIndex == representativeIndex
        representativeCentres = centres;
        representativeInnerSources = innerSources;
        representativeOuterSources = outerSources;
        representativeLambda = lambda;
        representativeGamma = gamma;
        representativeExteriorMagnitudes = exteriorMagnitudes{1};
        representativeInteriorMagnitudes = interiorMagnitudes{1};
    end
end

%% Gap trend and sphere-1 surface residual plots
summaryTable = table(gapList(:),repmat(sourceCount,caseCount,1), ...
    repmat(collocationCount,caseCount,1),retainedRanks, ...
    max(continuityAbsolute,[],2),max(continuityAbsolute,[],2)/velocityScale, ...
    max(normalAbsolute,[],2),max(normalAbsolute,[],2)/velocityScale, ...
    max(tractionAbsolute,[],2),max(tractionAbsolute,[],2)/tractionScale, ...
    max(maximumExteriorStrength,[],2),max(maximumInteriorStrength,[],2), ...
    max(maximumCoefficientComponent,[],2), ...
    'VariableNames',{'Delta','ActualNPerSphere','ActualMPerSphere', ...
    'RetainedRank','MaxContinuityAbs','MaxContinuityNormalized', ...
    'MaxNormalAbs','MaxNormalNormalized','MaxTractionAbs', ...
    'MaxTractionNormalized','MaxExteriorStrength', ...
    'MaxInteriorStrength','MaxCoefficientComponent'});
disp(summaryTable);

figure('Name','Two-drop gap sweep: independent surface conditions');
loglog(gapList,max(continuityAbsolute,[],2)/velocityScale,'o-', ...
    gapList,max(normalAbsolute,[],2)/velocityScale,'s-', ...
    gapList,max(tractionAbsolute,[],2)/tractionScale,'^-', ...
    'LineWidth',1.5);
set(gca,'XDir','normal');
axis tight; grid on;
xlabel('Gap delta');
ylabel('Largest pointwise error across both drops / stated physical scale');
legend('Velocity continuity','Normal velocity','Tangential traction', ...
    'Location','best');
title(sprintf('Two drops; actual N=%d, M=%d per drop; %d validation points per drop', ...
    sourceCount,collocationCount,validationCount));

figure('Name','Sphere 1 velocity-condition errors');
tiledlayout(1,2);
nexttile;
scatter3(representativeValidationPoints(:,1), ...
    representativeValidationPoints(:,2),representativeValidationPoints(:,3), ...
    22,representativeContinuityError,'filled');
axis equal; grid on; view(3); colorbar;
xlabel('x');ylabel('y');zlabel('z');
title('Velocity continuity error');
nexttile;
scatter3(representativeValidationPoints(:,1), ...
    representativeValidationPoints(:,2),representativeValidationPoints(:,3), ...
    22,representativeNormalError,'filled');
axis equal; grid on; view(3); colorbar;
xlabel('x');ylabel('y');zlabel('z');
title('Normal velocity error');
sgtitle(sprintf('Sphere 1, gap %.4g; velocity units',representativeGap));

figure('Name','Sphere 1 tangential-traction error');
scatter3(representativeValidationPoints(:,1), ...
    representativeValidationPoints(:,2),representativeValidationPoints(:,3), ...
    22,representativeTractionError,'filled');
axis equal; grid on; view(3); colorbar;
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('Sphere 1 tangential traction error; gap %.4g',representativeGap));

%% Sphere-1 proxy-source strengths at the representative gap
largestStrength = max([representativeExteriorMagnitudes; ...
    representativeInteriorMagnitudes]);
if largestStrength == 0
    logFloor = 1e-12;
else
    logFloor = max(largestStrength*1e-12,realmin);
end
exteriorLogStrength = log10(max(representativeExteriorMagnitudes,logFloor));
interiorLogStrength = log10(max(representativeInteriorMagnitudes,logFloor));
figure('Name','Sphere 1 Stokeslet strengths');
tiledlayout(1,2);
nexttile;
scatter3(representativeInnerSources{1}(:,1), ...
    representativeInnerSources{1}(:,2),representativeInnerSources{1}(:,3), ...
    25,exteriorLogStrength,'filled');
axis equal; grid on; view(3);
xlabel('x');ylabel('y');zlabel('z');
title('Exterior field: inner-proxy |\lambda_j|');
colourBar = colorbar;
colourBar.Label.String = 'log_{10} |\lambda_j|';
nexttile;
scatter3(representativeOuterSources{1}(:,1), ...
    representativeOuterSources{1}(:,2),representativeOuterSources{1}(:,3), ...
    25,interiorLogStrength,'filled');
axis equal; grid on; view(3);
xlabel('x');ylabel('y');zlabel('z');
title('Interior field: outer-proxy |\gamma_j|');
colourBar = colorbar;
colourBar.Label.String = 'log_{10} |\gamma_j|';
sgtitle(sprintf('Sphere 1, gap %.4g; per-node Stokeslet strength', ...
    representativeGap));

%% Physical flow in a plane through both sphere centres
% The plane contains the centreline and the active direction when those
% directions are not parallel. Colours show full 3D speed, while unit
% quiver arrows show the velocity direction projected into the plane.
planeHorizontal = [1 0 0];
planeVertical = activeDirection - ...
    (activeDirection*planeHorizontal.')*planeHorizontal;
if norm(planeVertical) < 1e-12
    planeVertical = [0 0 1];
end
planeVertical = planeVertical/norm(planeVertical);
horizontalExtent = abs(representativeCentres(2,1)) + ...
    (1+planePaddingRadii)*dropRadius;
verticalExtent = (1+planePaddingRadii)*dropRadius;
horizontalAxis = linspace(-horizontalExtent,horizontalExtent, ...
    planeGridResolution(1));
verticalAxis = linspace(-verticalExtent,verticalExtent, ...
    planeGridResolution(2));
[horizontalGrid,verticalGrid] = meshgrid(horizontalAxis,verticalAxis);
planePoints = horizontalGrid(:)*planeHorizontal + ...
    verticalGrid(:)*planeVertical;
insideOne = sum((planePoints-representativeCentres(1,:)).^2,2) <= ...
    dropRadius^2;
insideTwo = sum((planePoints-representativeCentres(2,:)).^2,2) <= ...
    dropRadius^2;
assert(~any(insideOne & insideTwo),'Positive gaps must separate the drops.');
planeVelocity = zeros(size(planePoints));
planeEvaluationBatchSize = 256;
for region = 0:2
    if region == 0
        pointIndices = find(~insideOne & ~insideTwo);
    elseif region == 1
        pointIndices = find(insideOne);
    else
        pointIndices = find(insideTwo);
    end
    for firstPoint = 1:planeEvaluationBatchSize:numel(pointIndices)
        selected = pointIndices(firstPoint:min( ...
            firstPoint+planeEvaluationBatchSize-1,numel(pointIndices)));
        targets = planePoints(selected,:);
        if region == 0
            velocityValues = ...
                drop_velocity_block(representativeInnerSources{1},targets,1)* ...
                    representativeLambda{1} + ...
                drop_velocity_block(representativeInnerSources{2},targets,1)* ...
                    representativeLambda{2};
        elseif region == 1
            velocityValues = ...
                drop_velocity_block( ...
                    representativeOuterSources{1},targets,interiorViscosity)* ...
                    representativeGamma{1} + ...
                drop_velocity_block(representativeInnerSources{2},targets,1)* ...
                    representativeLambda{2};
        else
            velocityValues = ...
                drop_velocity_block( ...
                    representativeOuterSources{2},targets,interiorViscosity)* ...
                    representativeGamma{2} + ...
                drop_velocity_block(representativeInnerSources{1},targets,1)* ...
                    representativeLambda{1};
        end
        planeVelocity(selected,:) = reshape(velocityValues,3,[]).';
    end
end
planeSpeed = reshape(vecnorm(planeVelocity,2,2),size(horizontalGrid));
horizontalVelocity = reshape(planeVelocity*planeHorizontal.',size(horizontalGrid));
verticalVelocity = reshape(planeVelocity*planeVertical.',size(horizontalGrid));
quiverRows = 1:planeQuiverStride:planeGridResolution(2);
quiverColumns = 1:planeQuiverStride:planeGridResolution(1);
quiverHorizontal = horizontalVelocity(quiverRows,quiverColumns);
quiverVertical = verticalVelocity(quiverRows,quiverColumns);
quiverSpeed = hypot(quiverHorizontal,quiverVertical);
quiverHorizontal = quiverHorizontal./max(quiverSpeed,eps);
quiverVertical = quiverVertical./max(quiverSpeed,eps);

figure('Name','Two-drop flow in a centreline cross section');
imagesc(horizontalAxis,verticalAxis,planeSpeed);
set(gca,'YDir','normal');axis equal tight;hold on;
quiver(horizontalGrid(quiverRows,quiverColumns), ...
    verticalGrid(quiverRows,quiverColumns), ...
    quiverHorizontal,quiverVertical,0.55,'k','LineWidth',0.8);
circleAngle = linspace(0,2*pi,300);
for body = 1:2
    plot(representativeCentres(body,1)+dropRadius*cos(circleAngle), ...
        dropRadius*sin(circleAngle),'w-','LineWidth',1.2);
end
colourBar = colorbar;
colourBar.Label.String = '|u| (full 3D speed)';
xlabel('Coordinate along sphere-centre line');
ylabel('Coordinate in centreline/active-direction plane');
title(sprintf('Two-drop flow; gap %.4g, actual N=%d per sphere', ...
    representativeGap,sourceCount));
