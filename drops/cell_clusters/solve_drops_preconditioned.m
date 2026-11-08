function [translations, alpha, gamma, lambda, diagnostics, factorBank] = ...
    solve_drops_preconditioned(innerSources, outerSources, ...
    surfacePoints, normals, S1, activeDirections, opt, factorBank)
%SOLVE_DROPS_PRECONDITIONED Solve the recompleted drop system by one-body GMRES.
% Geometry inputs are one cell per body. The physical strengths lambda_i
% sum to zero, and translations are U_i=-completionScale_i*sum(alpha_i).
% opt.completionScale may be one scalar or one value per body.
% The boundary-data unknown eta has [continuity;normal;tangent] rows.

if nargin == 0
    self_test();
    return
end
if nargin < 8
    factorBank = {};
end

bodyCount = numel(innerSources);
if isscalar(S1)
    S1 = repmat(S1,bodyCount,1);
end
activeDirections = activeDirections./vecnorm(activeDirections,2,2);
if isscalar(opt.interiorViscosity)
    opt.interiorViscosity = repmat(opt.interiorViscosity,bodyCount,1);
end
if isscalar(opt.completionScale)
    opt.completionScale = repmat(opt.completionScale,bodyCount,1);
end
assert(numel(opt.completionScale)==bodyCount && ...
    all(isfinite(opt.completionScale)) && all(opt.completionScale>0), ...
    'Provide one positive completion scale per body, or one shared scalar.');

sourceCounts = cellfun(@(points) size(points,1),innerSources);
targetCounts = cellfun(@(points) size(points,1),surfacePoints);
for body = 1:bodyCount
    assert(isequal(size(innerSources{body}),size(outerSources{body})) && ...
        isequal(size(surfacePoints{body}),size(normals{body})) && ...
        size(innerSources{body},2)==3 && size(surfacePoints{body},2)==3);
end

%% One-body preconditioning and scaled right-hand side
setupClock = tic;
factorizationTime = 0;
newFactorCount = 0;
factors = cell(1,bodyCount);
for body = 1:bodyCount
    anchor = innerSources{body}(1,:);
    innerRelative = innerSources{body}-anchor;
    outerRelative = outerSources{body}-anchor;
    surfaceRelative = surfacePoints{body}-anchor;
    selectedFactor = [];
    for candidateIndex = 1:numel(factorBank)
        candidate = factorBank{candidateIndex};
        if ~isequal(size(candidate.innerRelative),size(innerRelative)) || ...
                ~isequal(size(candidate.surfaceRelative),size(surfaceRelative)) || ...
                candidate.viscosity~=opt.interiorViscosity(body) || ...
                candidate.completionScale~=opt.completionScale(body) || ...
                candidate.velocityScale~=opt.velocityScale || ...
                candidate.tractionScale~=opt.tractionScale || ...
                candidate.tractionRowWeight~=opt.tractionRowWeight || ...
                candidate.svdTolerance~=opt.oneBodySvdTolerance || ...
                candidate.scaleColumns~=opt.scale_columns
            continue
        end
        geometryTolerance = 1e-12*max(1, ...
            max(abs(candidate.surfaceRelative),[],'all'));
        if max(abs(candidate.innerRelative-innerRelative),[],'all') < geometryTolerance && ...
                max(abs(candidate.outerRelative-outerRelative),[],'all') < geometryTolerance && ...
                max(abs(candidate.surfaceRelative-surfaceRelative),[],'all') < geometryTolerance && ...
                max(abs(candidate.normals-normals{body}),[],'all') < geometryTolerance
            selectedFactor = candidate;
            break
        end
    end
    if isempty(selectedFactor)
        factorClock = tic;
        selectedFactor = oneBodyPrecondDrop( ...
            innerSources{body},outerSources{body},surfacePoints{body}, ...
            normals{body},opt.interiorViscosity(body), ...
            opt.completionScale(body),opt.velocityScale,opt.tractionScale, ...
            opt.tractionRowWeight,opt.oneBodySvdTolerance,opt.scale_columns);
        factorizationTime = factorizationTime+toc(factorClock);
        factorBank{end+1} = selectedFactor; %#ok<AGROW>
        newFactorCount = newFactorCount+1;
    end
    factors{body} = selectedFactor;
end

rhs = zeros(6*sum(targetCounts),1);
rowOffset = 0;
for body = 1:bodyCount
    targetCount = targetCounts(body);
    direction = activeDirections(body,:).';
    normalDotDirection = normals{body}*direction;
    tangentData = S1(body)*[ ...
        (factors{body}.tangentOne*direction).*normalDotDirection, ...
        (factors{body}.tangentTwo*direction).*normalDotDirection];
    rhs(rowOffset+4*targetCount+(1:2*targetCount)) = ...
        opt.tractionRowWeight/opt.tractionScale* ...
        reshape(tangentData.',[],1);
    rowOffset = rowOffset+6*targetCount;
end
setupTime = toc(setupClock);

%% Solve the regularized identity-diagonal boundary-data system
% matvecDropMFS applies the local inverse and all cross-drop interactions.
operator = @(eta) matvecDropMFS( ...
    eta,innerSources,surfacePoints,normals,factors,opt);
solveClock = tic;
if norm(rhs)==0
    boundaryData = zeros(size(rhs));
    iterationCount = 0;
    residualHistory = 0;
    reportedResidual = 0;
else
    [boundaryData,iterationCount,residualHistory,reportedResidual] = ...
        helsing_gmres(operator,rhs,numel(rhs),opt.maxIterations, ...
        opt.krylovTolerance,0,0);
end
solveTime = toc(solveClock);

%% Map back to physical coefficients and translations
postprocessClock = tic;
alpha = cell(1,bodyCount);
gamma = cell(1,bodyCount);
lambda = cell(1,bodyCount);
translations = zeros(bodyCount,3);
forceSumNorm = zeros(bodyCount,1);
retainedRanks = zeros(bodyCount,1);
zeroColumnCounts = zeros(bodyCount,1);
coefficients = zeros(6*sum(sourceCounts),1);
discarded = zeros(size(boundaryData));
rowOffset = 0;
columnOffset = 0;
for body = 1:bodyCount
    targetCount = targetCounts(body);
    sourceCount = sourceCounts(body);
    rows = rowOffset+(1:6*targetCount);
    columns = columnOffset+(1:6*sourceCount);
    localCoefficients = applyOneBodyPinv(boundaryData(rows),factors{body});
    coefficients(columns) = localCoefficients;
    alpha{body} = localCoefficients(1:3*sourceCount);
    gamma{body} = localCoefficients(3*sourceCount+1:end);
    alphaNodes = reshape(alpha{body},3,sourceCount);
    lambdaNodes = alphaNodes-mean(alphaNodes,2);
    lambda{body} = lambdaNodes(:);
    translations(body,:) = ...
        (-opt.completionScale(body)*sum(alphaNodes,2)).';
    forceSumNorm(body) = norm(sum(lambdaNodes,2));
    retainedRanks(body) = factors{body}.rank;
    zeroColumnCounts(body) = factors{body}.zeroColumnCount;
    leftTranspose = factors{body}.leftTranspose;
    discarded(rows) = boundaryData(rows)- ...
        leftTranspose.'*(leftTranspose*boundaryData(rows));
    rowOffset = rowOffset+6*targetCount;
    columnOffset = columnOffset+6*sourceCount;
end
regularizedResidual = operator(boundaryData)-rhs;
originalResidual = drop_apply_A( ...
    coefficients,innerSources,surfacePoints,normals,factors,opt)-rhs;
normalizer = max(norm(rhs),realmin);
diagnostics.iterations = iterationCount;
diagnostics.residualHistory = residualHistory;
diagnostics.reportedResidual = reportedResidual;
diagnostics.regularizedResidual = norm(regularizedResidual)/normalizer;
diagnostics.originalResidual = norm(originalResidual)/normalizer;
diagnostics.discardedRangeNorm = norm(discarded)/normalizer;
diagnostics.residualIdentityError = ...
    norm(originalResidual-(regularizedResidual-discarded))/ ...
    max(1,norm(originalResidual));
diagnostics.converged = diagnostics.regularizedResidual <= opt.krylovTolerance;
diagnostics.rhsScaled = rhs;
diagnostics.retainedRanks = retainedRanks;
diagnostics.zeroColumnCounts = zeroColumnCounts;
diagnostics.forceSumNorm = forceSumNorm;
diagnostics.maxCoefficientComponent = max(abs(coefficients));
diagnostics.newFactorCount = newFactorCount;
diagnostics.factorizationTime = factorizationTime;
diagnostics.setupTime = setupTime;
diagnostics.solveTime = solveTime;
diagnostics.postprocessTime = toc(postprocessClock);
end

function self_test()
% Exercise factor reuse, recovery, and the zero-forcing GMRES guard.
axisNodes = [eye(3);-eye(3)];
corners = [1 1 1;1 1 -1;1 -1 1;-1 1 1; ...
           1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]/sqrt(3);
centres = [-1.6 0 0;1.6 0 0];
inner = cell(1,2); outer = cell(1,2);
surface = cell(1,2); normals = cell(1,2);
for body = 1:2
    inner{body} = centres(body,:)+0.6*axisNodes;
    outer{body} = centres(body,:)+2*axisNodes;
    surface{body} = centres(body,:)+corners;
    normals{body} = corners;
end
opt.interiorViscosity = 2;
opt.completionScale = 1/(8*pi);
opt.velocityScale = 1;
opt.tractionScale = 1;
opt.tractionRowWeight = 1;
opt.oneBodySvdTolerance = 1e-12;
opt.krylovTolerance = 1e-8;
opt.maxIterations = 20;
opt.interactionBatchSize = 3;
opt.scale_columns = false;
[U,alpha,gamma,lambda,info,bank] = solve_drops_preconditioned( ...
    inner,outer,surface,normals,1,repmat([0 0 1],2,1),opt);
assert(isequal(size(U),[2 3]) && numel(alpha)==2 && numel(gamma)==2);
assert(isscalar(bank) && info.newFactorCount==1);
assert(all(info.forceSumNorm < 1e-12) && ...
    all(cellfun(@(x) numel(x)==18,lambda)));
assert(info.converged && info.iterations>0);
assert(info.residualIdentityError < 1e-10);
[zeroU,~,~,~,zeroInfo,bank] = solve_drops_preconditioned( ...
    inner,outer,surface,normals,0,repmat([0 0 1],2,1),opt,bank);
assert(zeroInfo.iterations==0 && zeroInfo.newFactorCount==0 && ...
    norm(zeroU)==0 && isscalar(bank));
for body = 1:2
    shift = (-1)^body*0.1*[1 0 0];
    inner{body} = inner{body}+shift;
    outer{body} = outer{body}+shift;
    surface{body} = surface{body}+shift;
end
[~,~,~,~,shiftInfo,bank] = solve_drops_preconditioned( ...
    inner,outer,surface,normals,0,repmat([0 0 1],2,1),opt,bank);
assert(shiftInfo.newFactorCount==0 && isscalar(bank));

% Unequal radii require different local completions and separate factors.
radii = [0.8;1.2];
centres = [-2 0 0;2 0 0];
for body = 1:2
    inner{body} = centres(body,:)+0.6*radii(body)*axisNodes;
    outer{body} = centres(body,:)+2*radii(body)*axisNodes;
    surface{body} = centres(body,:)+radii(body)*corners;
    normals{body} = corners;
end
opt.interiorViscosity = [2;3];
opt.completionScale = 1./(8*pi*radii);
[heterogeneousU,heterogeneousAlpha,~,~,heterogeneousInfo,bank] = ...
    solve_drops_preconditioned(inner,outer,surface,normals, ...
    [1;0.7],repmat([0 0 1],2,1),opt,{});
assert(heterogeneousInfo.newFactorCount==2 && numel(bank)==2);
for body = 1:2
    expectedU = -opt.completionScale(body)* ...
        sum(reshape(heterogeneousAlpha{body},3,[]),2).';
    assert(norm(heterogeneousU(body,:)-expectedU)<1e-12);
    assert(bank{body}.completionScale==opt.completionScale(body));
end
assert(heterogeneousInfo.residualIdentityError<1e-10);
fprintf('solve_drops_preconditioned self-test passed.\n');
end
