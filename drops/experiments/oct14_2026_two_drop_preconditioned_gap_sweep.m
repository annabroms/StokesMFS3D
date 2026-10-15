% OCT14_2026_TWO_DROP_PRECONDITIONED_GAP_SWEEP Compare unrestarted one-body
% GMRES with the dense SVD reference at every gap. Run this script to start
% the sweep. N, M, and validation counts below are PER DROP.

%% Parameters (matching oct14_2026_two_drop_gap_sweep.m)
clear;
close all;
dropRadius = 1;
innerProxyRadius = 0.6;
outerProxyRadius = 2;
S1 = 1;
interiorViscosity = 2;           % ratio to exterior viscosity mu_0=1
requestedN = 600;
collocationToSourceRatio = 1.2;
validationResolution = 1000;    % at least this many independent points
gapList = [1 0.5 0.2 0.1 0.05 0.02 0.01];
activeDirection = [0 0 1];       % identical on both drops
tractionRowWeight = 1;
scale_columns = false;
denseSvdTolerance = 1e-13;
oneBodySvdTolerance = 1e-13;
krylovTolerance = 1e-10;
residualIdentityReference = 1e-8; % plotted guide, never a stopping test
maxIterations = 150;             % total iterations; GMRES is unrestarted
% These are target points handled per chunk, not grid resolutions. Larger
% chunks use more temporary memory; they do not change N, M, or the model.
interactionTargetBatchSize = 128; % collocation targets per GMRES matvec chunk
validationTargetBatchSize = 128;  % independent check targets per chunk

%% Paths, designs, and physical scaling
experimentDirectory = fileparts(mfilename('fullpath'));
dropsDirectory = fileparts(experimentDirectory);
repositoryDirectory = fileparts(dropsDirectory);
addpath(dropsDirectory,fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'common','gmres'), ...
    fullfile(repositoryDirectory,'geometry','spheres'));

activeDirection = activeDirection/norm(activeDirection);
activeDirections = repmat(activeDirection,2,1);

catalogueFile = fullfile(repositoryDirectory,'geometry','spheres', ...
    'spherical_designs','filelist.txt');
catalogueID = fopen(catalogueFile,'r');
if catalogueID<0
    error('Could not read the spherical-design catalogue.');
end
catalogueData = textscan(catalogueID,'%*[^.].%d');
fclose(catalogueID);
availableCounts = catalogueData{1};
[sourceDirections,~] = get_sphdesign(requestedN);
sourceCount = size(sourceDirections,1);
requestedM = ceil(collocationToSourceRatio*sourceCount);
collocationIndex = find(availableCounts>=requestedM,1);
validationIndex = find(availableCounts>=validationResolution,1);
if isempty(collocationIndex) || isempty(validationIndex)
    error('No design supplies the requested collocation or validation count.');
end
[collocationDirections,~] = ...
    get_sphdesign(availableCounts(collocationIndex));
[validationDirections,~] = ...
    get_sphdesign(availableCounts(validationIndex));
collocationCount = size(collocationDirections,1);
validationCount = size(validationDirections,1);
validationCountShortfall = max(0,collocationCount+1-validationCount);

% Rotate the independent validation design away from collocation nodes.
azimuth = 0.37;
elevation = 0.23;
rotationZ = [cos(azimuth) -sin(azimuth) 0; ...
             sin(azimuth) cos(azimuth) 0;0 0 1];
rotationY = [cos(elevation) 0 sin(elevation); ...
             0 1 0;-sin(elevation) 0 cos(elevation)];
validationDirections = (rotationZ*rotationY*validationDirections.').';
validationMaximumDot = max(validationDirections*collocationDirections.',[],'all');
validationOverlapExcess = max(0,validationMaximumDot-(1-1e-12));

if S1==0
    velocityScale = 1;
    tractionScale = 1;
    residualLabel = 'unit-scaled';
else
    velocityScale = abs(S1)*dropRadius; % mu_0=1
    tractionScale = abs(S1);
    residualLabel = 'relative';
end
completionScale = 1/(8*pi*dropRadius);
physicalScales = [velocityScale velocityScale tractionScale];
fprintf('Per drop: requested N=%d, actual N=%d; requested M>=%d, actual M=%d (M/N=%.4f).\n', ...
    requestedN,sourceCount,requestedM,collocationCount, ...
    collocationCount/sourceCount);
fprintf('Independent validation: requested >=%d, actual %d per drop.\n', ...
    validationResolution,validationCount);
fprintf('Validation diagnostics: point-count shortfall %d; overlap excess %.3e (maximum direction dot product %.16g).\n', ...
    validationCountShortfall,validationOverlapExcess,validationMaximumDot);
fprintf('Scales: velocity %.6g, traction %.6g; traction row weight %.6g.\n', ...
    velocityScale,tractionScale,tractionRowWeight);
fprintf('Cutoffs: dense SVD %.3g, one-body SVD %.3g, GMRES %.3g; max %d unrestarted iterations.\n', ...
    denseSvdTolerance,oneBodySvdTolerance,krylovTolerance,maxIterations);
fprintf('Dense matrix per gap: %d by %d, %.2f GiB for A alone.\n', ...
    12*collocationCount,12*sourceCount, ...
    8*(12*collocationCount)*(12*sourceCount)/2^30);
if scale_columns
    fprintf(['Column scaling is local for one-body SVD and global for ', ...
        'dense SVD; compare their retained ranks when interpreting results.\n']); %#ok<UNRCH>
end

solverOpt.interiorViscosity = interiorViscosity;
solverOpt.completionScale = completionScale;
solverOpt.velocityScale = velocityScale;
solverOpt.tractionScale = tractionScale;
solverOpt.tractionRowWeight = tractionRowWeight;
solverOpt.oneBodySvdTolerance = oneBodySvdTolerance;
solverOpt.krylovTolerance = krylovTolerance;
solverOpt.maxIterations = maxIterations;
solverOpt.interactionBatchSize = interactionTargetBatchSize;
solverOpt.scale_columns = scale_columns;

caseCount = numel(gapList);
denseTranslations = zeros(caseCount,2,3);
iterativeTranslations = zeros(caseCount,2,3);
denseValidation = zeros(caseCount,2,3);
iterativeValidation = zeros(caseCount,2,3);
denseRanks = zeros(caseCount,1);
oneBodyRanks = zeros(caseCount,1);
iterations = zeros(caseCount,1);
residualHistories = cell(caseCount,1);
denseAssemblyTimes = zeros(caseCount,1);
denseSolveTimes = zeros(caseCount,1);
iterativeSetupTimes = zeros(caseCount,1);
iterativeSolveTimes = zeros(caseCount,1);
iterativePostprocessTimes = zeros(caseCount,1);
denseCollocationResidual = zeros(caseCount,1);
denseMaxCoefficient = zeros(caseCount,1);
iterativeMaxCoefficient = zeros(caseCount,1);
denseRhsVariation = zeros(caseCount,1);
denseRhsByGap = cell(caseCount,1);
iterativeRhsMismatch = zeros(caseCount,1);
denseZeroColumnCount = zeros(caseCount,1);
oneBodyZeroColumnCount = zeros(caseCount,1);
residualIdentityDefect = zeros(caseCount,1);
denseForceFreeDefect = zeros(caseCount,2);
iterativeForceFreeDefect = zeros(caseCount,2);
regularizedResidual = zeros(caseCount,1);
originalResidual = zeros(caseCount,1);
discardedRangeNorm = zeros(caseCount,1);
translationAbsoluteDifference = zeros(caseCount,2);
translationRelativeDifference = NaN(caseCount,2);

%% Dense reference sweep first: release large SVD arrays before GMRES setup
for caseIndex = 1:caseCount
    delta = gapList(caseIndex);
    centreDistance = 2*dropRadius+delta;
    centres = [-centreDistance/2 0 0;centreDistance/2 0 0];
    innerSources = cell(1,2); outerSources = cell(1,2);
    collocationPoints = cell(1,2); validationPoints = cell(1,2);
    normals = cell(1,2); validationNormals = cell(1,2);
    for body = 1:2
        innerSources{body} = centres(body,:)+ ...
            innerProxyRadius*sourceDirections;
        outerSources{body} = centres(body,:)+ ...
            outerProxyRadius*sourceDirections;
        collocationPoints{body} = centres(body,:)+ ...
            dropRadius*collocationDirections;
        validationPoints{body} = centres(body,:)+ ...
            dropRadius*validationDirections;
        normals{body} = collocationDirections;
        validationNormals{body} = validationDirections;
    end

    assemblyClock = tic;
    [denseMatrix,tangentOne,tangentTwo] = drop_two_body_matrix( ...
        innerSources,outerSources,collocationPoints,normals, ...
        interiorViscosity,completionScale);
    rhs = zeros(12*collocationCount,1);
    for body = 1:2
        normalDotDirection = normals{body}*activeDirection.';
        tangentData = S1*[ ...
            (tangentOne{body}*activeDirection.').*normalDotDirection, ...
            (tangentTwo{body}*activeDirection.').*normalDotDirection];
        offset = (body-1)*6*collocationCount;
        tractionRows = offset+4*collocationCount+(1:2*collocationCount);
        rhs(tractionRows) = tractionRowWeight/tractionScale* ...
            reshape(tangentData.',[],1);
        velocityRows = offset+(1:4*collocationCount);
        denseMatrix(velocityRows,:) = ...
            denseMatrix(velocityRows,:)/velocityScale;
        denseMatrix(tractionRows,:) = ...
            tractionRowWeight*denseMatrix(tractionRows,:)/tractionScale;
    end
    if caseIndex==1
        referenceRhs = rhs;
    else
        denseRhsVariation(caseIndex) = ...
            norm(rhs-referenceRhs,inf)/max(1,norm(referenceRhs,inf));
    end
    denseRhsByGap{caseIndex} = rhs;
    denseAssemblyTimes(caseIndex) = toc(assemblyClock);

    solveClock = tic;
    if scale_columns
        globalColumnScales = vecnorm(denseMatrix,2,1).'; %#ok<UNRCH>
        denseZeroColumnCount(caseIndex) = nnz(globalColumnScales==0);
        globalColumnScales(globalColumnScales==0) = 1;
        denseMatrix = denseMatrix./globalColumnScales.';
    else
        globalColumnScales = ones(size(denseMatrix,2),1);
    end
    [leftVectors,singularMatrix,rightVectors] = svd(denseMatrix,'econ');
    singularValues = diag(singularMatrix);
    retained = singularValues>denseSvdTolerance*singularValues(1);
    denseRanks(caseIndex) = nnz(retained);
    scaledCoefficients = rightVectors(:,retained)* ...
        ((leftVectors(:,retained).'*rhs)./singularValues(retained));
    denseCoefficients = scaledCoefficients./globalColumnScales;
    denseMaxCoefficient(caseIndex) = max(abs(denseCoefficients));
    denseCollocationResidual(caseIndex) = ...
        norm(denseMatrix*scaledCoefficients-rhs)/max(norm(rhs),realmin);
    denseSolveTimes(caseIndex) = toc(solveClock);
    clear denseMatrix leftVectors singularMatrix rightVectors

    denseLambda = cell(1,2); denseGamma = cell(1,2);
    denseU = zeros(2,3);
    for body = 1:2
        columns = (body-1)*6*sourceCount+(1:6*sourceCount);
        local = denseCoefficients(columns);
        alphaNodes = reshape(local(1:3*sourceCount),3,sourceCount);
        lambdaNodes = alphaNodes-mean(alphaNodes,2);
        denseLambda{body} = lambdaNodes(:);
        denseForceFreeDefect(caseIndex,body) = ...
            norm(sum(lambdaNodes,2))/ ...
            max(realmin,sum(vecnorm(lambdaNodes,2,1)));
        denseGamma{body} = local(3*sourceCount+1:end);
        denseU(body,:) = (-completionScale*sum(alphaNodes,2)).';
    end
    denseTranslations(caseIndex,:,:) = denseU;
    denseValidation(caseIndex,:,:) = drop_validate_surface_conditions( ...
        denseLambda,denseGamma,denseU,innerSources,outerSources, ...
        validationPoints,validationNormals,activeDirections,S1, ...
        interiorViscosity,validationTargetBatchSize);
    fprintf('Dense delta=%.5g: rank %d/%d, setup %.3fs, SVD solve %.3fs, scaled fit %.3e.\n', ...
        delta,denseRanks(caseIndex),12*sourceCount, ...
        denseAssemblyTimes(caseIndex),denseSolveTimes(caseIndex), ...
        denseCollocationResidual(caseIndex));
end

%% Iterative sweep: the first solve factors one sphere, later gaps reuse it
factorBank = {};
for caseIndex = 1:caseCount
    delta = gapList(caseIndex);
    centreDistance = 2*dropRadius+delta;
    centres = [-centreDistance/2 0 0;centreDistance/2 0 0];
    innerSources = cell(1,2); outerSources = cell(1,2);
    collocationPoints = cell(1,2); validationPoints = cell(1,2);
    normals = cell(1,2); validationNormals = cell(1,2);
    for body = 1:2
        innerSources{body} = centres(body,:)+ ...
            innerProxyRadius*sourceDirections;
        outerSources{body} = centres(body,:)+ ...
            outerProxyRadius*sourceDirections;
        collocationPoints{body} = centres(body,:)+ ...
            dropRadius*collocationDirections;
        validationPoints{body} = centres(body,:)+ ...
            dropRadius*validationDirections;
        normals{body} = collocationDirections;
        validationNormals{body} = validationDirections;
    end

    [iterativeU,~,iterativeGamma,iterativeLambda,info,factorBank] = ...
        solve_drops_preconditioned(innerSources,outerSources, ...
        collocationPoints,normals,S1,activeDirections,solverOpt,factorBank);
    iterativeRhsMismatch(caseIndex) = ...
        norm(info.rhsScaled-denseRhsByGap{caseIndex},inf)/ ...
        max(1,norm(denseRhsByGap{caseIndex},inf));
    expectedNewFactors = double(caseIndex==1);
    % Reuse is part of the timing experiment, unlike accuracy diagnostics.
    assert(isscalar(factorBank) && ...
        info.newFactorCount==expectedNewFactors, ...
        'Identical translated spheres must reuse one factorization across gaps.');
    iterativeTranslations(caseIndex,:,:) = iterativeU;
    iterativeValidation(caseIndex,:,:) = ...
        drop_validate_surface_conditions(iterativeLambda,iterativeGamma, ...
        iterativeU,innerSources,outerSources,validationPoints, ...
        validationNormals,activeDirections,S1,interiorViscosity, ...
        validationTargetBatchSize);
    iterations(caseIndex) = info.iterations;
    residualHistories{caseIndex} = info.residualHistory;
    oneBodyRanks(caseIndex) = info.retainedRanks(1);
    oneBodyZeroColumnCount(caseIndex) = max(info.zeroColumnCounts);
    iterativeSetupTimes(caseIndex) = info.setupTime;
    iterativeSolveTimes(caseIndex) = info.solveTime;
    iterativePostprocessTimes(caseIndex) = info.postprocessTime;
    regularizedResidual(caseIndex) = info.regularizedResidual;
    originalResidual(caseIndex) = info.originalResidual;
    discardedRangeNorm(caseIndex) = info.discardedRangeNorm;
    iterativeMaxCoefficient(caseIndex) = info.maxCoefficientComponent;
    residualIdentityDefect(caseIndex) = info.residualIdentityError;
    for body = 1:2
        lambdaNodes = reshape(iterativeLambda{body},3,sourceCount);
        iterativeForceFreeDefect(caseIndex,body) = ...
            norm(sum(lambdaNodes,2))/ ...
            max(realmin,sum(vecnorm(lambdaNodes,2,1)));
    end

    for body = 1:2
        referenceU = reshape(denseTranslations(caseIndex,body,:),1,3);
        translationAbsoluteDifference(caseIndex,body) = ...
            norm(iterativeU(body,:)-referenceU);
        referenceMagnitude = norm(referenceU);
        if referenceMagnitude>sqrt(eps)*velocityScale
            translationRelativeDifference(caseIndex,body) = ...
                translationAbsoluteDifference(caseIndex,body)/referenceMagnitude;
        end
    end
    fprintf(['Iterative delta=%.5g: rank %d/%d, %d iterations, ', ...
        'GMRES relres %.3e (recomputed %.3e), original %.3e, ', ...
        'discarded %.3e; setup %.3fs (SVD %.3fs), solve %.3fs, ', ...
        'recovery %.3fs.\n'],delta,oneBodyRanks(caseIndex), ...
        6*sourceCount,iterations(caseIndex),info.reportedResidual, ...
        regularizedResidual(caseIndex),originalResidual(caseIndex), ...
        discardedRangeNorm(caseIndex),info.setupTime, ...
        info.factorizationTime,info.solveTime,info.postprocessTime);
    if ~info.converged
        fprintf('  GMRES did not reach the configured tolerance.\n');
    end
    fprintf(['  Maximum coefficient component: dense %.3e, ', ...
        'iterative %.3e; force-free sum norms [%.3e %.3e].\n'], ...
        denseMaxCoefficient(caseIndex),iterativeMaxCoefficient(caseIndex), ...
        info.forceSumNorm);
    fprintf(['  Diagnostic defects: dense RHS variation %.3e, ', ...
        'iterative RHS mismatch %.3e, residual identity %.3e.\n'], ...
        denseRhsVariation(caseIndex),iterativeRhsMismatch(caseIndex), ...
        residualIdentityDefect(caseIndex));
    for body = 1:2
        referenceU = reshape(denseTranslations(caseIndex,body,:),1,3);
        denseErrors = squeeze(denseValidation(caseIndex,body,:)).';
        iterativeErrors = squeeze(iterativeValidation(caseIndex,body,:)).';
        fprintf(['  Drop %d: U_dense=[%.4g %.4g %.4g], ', ...
            'U_GMRES=[%.4g %.4g %.4g]; |Delta U|=%.3e, ', ...
            'relative=%.3e (NaN if reference U is negligible).\n'], ...
            body,referenceU,iterativeU(body,:), ...
            translationAbsoluteDifference(caseIndex,body), ...
            translationRelativeDifference(caseIndex,body));
        fprintf(['    Validation absolute [continuity normal traction]: ', ...
            'dense [%.3e %.3e %.3e], iterative [%.3e %.3e %.3e].\n'], ...
            denseErrors,iterativeErrors);
        fprintf(['    Validation %s / [%.3g %.3g %.3g]: ', ...
            'dense [%.3e %.3e %.3e], iterative [%.3e %.3e %.3e].\n'], ...
            residualLabel,physicalScales, ...
            denseErrors./physicalScales, ...
            iterativeErrors./physicalScales);
    end
end

%% Compact comparison tables and gap plots
denseMaximum = squeeze(max(denseValidation,[],2));
iterativeMaximum = squeeze(max(iterativeValidation,[],2));
summaryTable = table(gapList(:),repmat(sourceCount,caseCount,1), ...
    repmat(collocationCount,caseCount,1),iterations,oneBodyRanks, ...
    denseRanks,regularizedResidual,originalResidual,discardedRangeNorm, ...
    denseMaxCoefficient,iterativeMaxCoefficient, ...
    denseMaximum(:,1)/velocityScale, ...
    iterativeMaximum(:,1)/velocityScale, ...
    denseMaximum(:,2)/velocityScale, ...
    iterativeMaximum(:,2)/velocityScale, ...
    denseMaximum(:,3)/tractionScale, ...
    iterativeMaximum(:,3)/tractionScale, ...
    'VariableNames',{'Delta','ActualNPerDrop','ActualMPerDrop', ...
    'GMRESIterations','OneBodyRank','DenseRank','GMRESResidual', ...
    'OriginalResidual','DiscardedRange', ...
    'DenseMaxCoefficient','IterativeMaxCoefficient', ...
    'DenseContinuity','IterativeContinuity', ...
    'DenseNormal','IterativeNormal','DenseTraction','IterativeTraction'});
disp(summaryTable);
timingTable = table(gapList(:),denseAssemblyTimes,denseSolveTimes, ...
    iterativeSetupTimes,iterativeSolveTimes,iterativePostprocessTimes, ...
    'VariableNames',{'Delta','DenseSetupSeconds','DenseSolveSeconds', ...
    'IterativeSetupSeconds','IterativeSolveSeconds', ...
    'IterativePostprocessSeconds'});
disp(timingTable);

[sortedDelta,plotOrder] = sort(gapList(:),'ascend');
figure('Name','Two-drop one-body GMRES iterations against gap');
semilogx(sortedDelta,iterations(plotOrder),'o-','LineWidth',1.5);
set(gca,'XDir','normal');axis tight;grid on;
xlabel('Gap delta');ylabel('Unrestarted GMRES iterations');
title(sprintf('Two drops; actual N=%d, M=%d per drop', ...
    sourceCount,collocationCount));
legend('Total Krylov iterations for the regularized system H eta = b', ...
    'Location','best');

figure('Name','Two-drop translation differences against gap', ...
    'Position',[100 100 1300 520]);
tiledlayout(1,2);
nexttile;
loglog(sortedDelta,max(translationAbsoluteDifference(plotOrder,1),realmin), ...
    'o-',sortedDelta, ...
    max(translationAbsoluteDifference(plotOrder,2),realmin),'s-', ...
    'LineWidth',1.4);
axis tight;grid on;xlabel('Gap delta');ylabel('Velocity units');
title('Absolute translation difference');
legend('Drop 1: ||U_{GMRES}-U_{dense}||_2', ...
    'Drop 2: ||U_{GMRES}-U_{dense}||_2','Location','best');
nexttile;
loglog(sortedDelta,max(translationRelativeDifference(plotOrder,1),realmin), ...
    'o-',sortedDelta, ...
    max(translationRelativeDifference(plotOrder,2),realmin),'s-', ...
    'LineWidth',1.4);
axis tight;grid on;xlabel('Gap delta');ylabel('Dimensionless difference');
title('Relative translation difference');
legend('Drop 1: ||Delta U||_2 / ||U_{dense}||_2, if reference exceeds sqrt(eps) velocity scale', ...
    'Drop 2: ||Delta U||_2 / ||U_{dense}||_2, if reference exceeds sqrt(eps) velocity scale', ...
    'Location','best');

figure('Name','Two-drop unrestarted GMRES residual histories', ...
    'Position',[100 100 1250 620]);
historyLabels = cell(caseCount+1,1);
historyLines = gobjects(caseCount+1,1);
for caseIndex = 1:caseCount
    history = residualHistories{caseIndex};
    historyLines(caseIndex) = semilogy( ...
        1:numel(history),max(history,eps),'-', ...
        'LineWidth',1.1);
    hold on;
    historyLabels{caseIndex} = sprintf( ...
        'Gap %.3g: Arnoldi estimate of ||H eta_k-b_{scaled}||_2 / ||b_{scaled}||_2', ...
        gapList(caseIndex));
end
historyLines(end) = yline(krylovTolerance,'k--');
historyLabels{end} = 'Configured GMRES relative-residual tolerance';
axis tight;grid on;
xlabel('GMRES iteration');
ylabel('Estimated ||H eta_k-b_{scaled}||_2 / ||b_{scaled}||_2 (floor eps)');
legend(historyLines,historyLabels,'Location','best');

figure('Name','Two-drop scaled system residuals against gap', ...
    'Position',[100 100 1100 600]);
systemLines = loglog(sortedDelta,max(regularizedResidual(plotOrder),realmin),'o-', ...
    sortedDelta,max(originalResidual(plotOrder),realmin),'s-', ...
    sortedDelta,max(discardedRangeNorm(plotOrder),realmin),'^-', ...
    'LineWidth',1.5);
hold on;
systemToleranceLine = yline(krylovTolerance,'k--');
set(gca,'XDir','normal');axis tight;grid on;
xlabel('Gap delta');ylabel('Row-scaled Euclidean residual / ||b_{scaled}||_2');
legend([systemLines(:);systemToleranceLine], ...
    {'GMRES equation: ||H eta-b_{scaled}||_2 / ||b_{scaled}||_2', ...
    'Original collocation system: ||A_{scaled} z-b_{scaled}||_2 / ||b_{scaled}||_2', ...
    'Discarded local range: ||(I-J) eta||_2 / ||b_{scaled}||_2', ...
    'Configured GMRES relative-residual tolerance'}, ...
    'Location','best');

conditionNames = {'Velocity continuity','Normal velocity', ...
    'Tangential traction'};
if S1==0
    velocityScaleDescription = 'unit velocity scale';
    tractionScaleDescription = 'unit traction scale';
else
    velocityScaleDescription = '|S_1| a / mu_0';
    tractionScaleDescription = '|S_1|';
end
conditionMetrics = { ...
    ['max_{i,m} ||u_{ext}-u_{int}||_2 / ' velocityScaleDescription], ...
    ['max_{i,m} |n dot (u_{ext}-U_i)| / ' velocityScaleDescription], ...
    ['max_{i,m} ||Pi_{tan}(T_{ext}-T_{int})-s||_2 / ' ...
        tractionScaleDescription]};
figure('Name','Two-drop independent surface residuals against gap', ...
    'Position',[100 100 1150 950]);
tiledlayout(3,1);
for condition = 1:3
    nexttile;
    denseCurve = denseMaximum(plotOrder,condition)/ ...
        physicalScales(condition);
    iterativeCurve = iterativeMaximum(plotOrder,condition)/ ...
        physicalScales(condition);
    loglog(sortedDelta,max(denseCurve,realmin),'o-', ...
        sortedDelta,max(iterativeCurve,realmin),'s-', ...
        'LineWidth',1.5);
    set(gca,'XDir','normal');axis tight;grid on;
    xlabel('Gap delta');ylabel(['Maximum ' residualLabel ' residual']);
    title(conditionNames{condition});
    legend(['Dense SVD: ' conditionMetrics{condition}], ...
        ['One-body GMRES: ' conditionMetrics{condition}], ...
        'Location','best');
end

%% Numerical diagnostics are shown rather than used as stopping assertions
figure('Name','Two-drop numerical agreement and violation diagnostics', ...
    'Position',[100 100 1450 900]);
tiledlayout(3,2,'TileSpacing','compact');

nexttile;
semilogx(sortedDelta,denseRhsVariation(plotOrder),'o-', ...
    sortedDelta,iterativeRhsMismatch(plotOrder),'s-','LineWidth',1.4);
axis tight;grid on;xlabel('Gap delta');ylabel('Relative infinity-norm defect');
title('Boundary data and row scaling');
legend('Dense scaled RHS variation: ||b(delta)-b(first)||_inf / max(1,||b(first)||_inf)', ...
    'Dense vs GMRES scaled RHS: ||b_GMRES-b_dense||_inf / max(1,||b_dense||_inf)', ...
    'Location','best');

nexttile;
identityLine = semilogx( ...
    sortedDelta,residualIdentityDefect(plotOrder),'o-', ...
    'LineWidth',1.4);
hold on;
identityReferenceLine = yline(residualIdentityReference,'k--');
axis tight;grid on;xlabel('Gap delta');ylabel('Dimensionless identity defect');
title('Rectangular regularization identity');
legend([identityLine;identityReferenceLine], ...
    'Row-scaled identity: ||Az-b-[H eta-b-(I-J)eta]||_2 / max(1,||Az-b||_2)', ...
    'Visual reference only: 1e-8 (does not stop the sweep)', ...
    'Location','best');

nexttile;
semilogx(sortedDelta,max(denseForceFreeDefect(plotOrder,:),[],2),'o-', ...
    sortedDelta,max(iterativeForceFreeDefect(plotOrder,:),[],2),'s-', ...
    'LineWidth',1.4);
axis tight;grid on;xlabel('Gap delta');ylabel('Largest relative net strength');
title('Force-free projection');
legend('Dense: max_i ||sum_j lambda_{ij}||_2 / max(realmin,sum_j ||lambda_{ij}||_2)', ...
    'GMRES: max_i ||sum_j lambda_{ij}||_2 / max(realmin,sum_j ||lambda_{ij}||_2)', ...
    'Location','best');

nexttile;
semilogx(sortedDelta,repmat(validationCountShortfall,caseCount,1), ...
    'o-','LineWidth',1.4);
axis tight;grid on;xlabel('Gap delta');ylabel('Missing validation points');
title('Independent-grid size');
legend('max(0, M_collocation + 1 - M_validation) points per drop', ...
    'Location','best');

nexttile;
semilogx(sortedDelta,repmat(validationOverlapExcess,caseCount,1), ...
    'o-','LineWidth',1.4);
axis tight;grid on;xlabel('Gap delta');ylabel('Dot-product excess');
title('Validation-node separation');
legend('max(0, max_{v,c}(v dot c) - (1-1e-12)); 0 means no shared direction', ...
    'Location','best');

nexttile;
semilogx(sortedDelta,denseZeroColumnCount(plotOrder),'o-', ...
    sortedDelta,oneBodyZeroColumnCount(plotOrder),'s-', ...
    'LineWidth',1.4);
axis tight;grid on;xlabel('Gap delta');ylabel('Zero-norm column count');
title('Column scaling diagnostics');
legend('Dense matrix: number of zero-norm columns if scaling enabled', ...
    'One-body matrix: max_i zero-norm columns if scaling enabled', ...
    'Location','best');
