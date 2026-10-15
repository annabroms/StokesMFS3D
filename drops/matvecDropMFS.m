function response = matvecDropMFS(boundaryData, innerSources, ...
    surfacePoints, normals, factors, opt)
%MATVECDROPMFS Apply the identity-diagonal one-body regularized system.
% For each body z_i=A_ii^+ eta_i; only j~=i interactions are evaluated.
% Thus H_ii=I rather than the rank-deficient projector A_ii*A_ii^+.

if nargin == 0
    self_test();
    return
end

bodyCount = numel(innerSources);
targetCounts = cellfun(@(points) size(points,1),surfacePoints);
assert(iscolumn(boundaryData) && ...
    numel(boundaryData)==6*sum(targetCounts) && ...
    numel(factors)==bodyCount && numel(normals)==bodyCount);

lambda = cell(1,bodyCount);
rowOffset = 0;
for body = 1:bodyCount
    targetCount = targetCounts(body);
    sourceCount = size(innerSources{body},1);
    localData = boundaryData(rowOffset+(1:6*targetCount));
    localCoefficients = applyOneBodyPinv(localData,factors{body});
    assert(numel(localCoefficients)==6*sourceCount);
    alphaNodes = reshape(localCoefficients(1:3*sourceCount),3,sourceCount);
    lambda{body} = reshape(alphaNodes-mean(alphaNodes,2),[],1);
    rowOffset = rowOffset+6*targetCount;
end
response = boundaryData+drop_apply_interactions( ...
    lambda,innerSources,surfacePoints,normals,factors,opt);
end

function self_test()
% Verify identity diagonal blocks and the rectangular residual identity.
axisNodes = [eye(3);-eye(3)];
corners = [1 1 1;1 1 -1;1 -1 1;-1 1 1; ...
           1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]/sqrt(3);
centres = [-1.6 0 0;1.6 0 0];
inner = cell(1,2); outer = cell(1,2);
surface = cell(1,2); normals = cell(1,2); factors = cell(1,2);
for body = 1:2
    inner{body} = centres(body,:)+0.6*axisNodes;
    outer{body} = centres(body,:)+2*axisNodes;
    surface{body} = centres(body,:)+corners;
    normals{body} = corners;
    factors{body} = oneBodyPrecondDrop(inner{body},outer{body}, ...
        surface{body},normals{body},2,1/(8*pi),3,5,1.7,1e-12,false);
end
opt.interiorViscosity = [2 2];
opt.velocityScale = 3;
opt.tractionScale = 5;
opt.tractionRowWeight = 1.7;
opt.interactionBatchSize = 3;
eta = sin((1:96).');
eta(49:end) = 0;
response = matvecDropMFS(eta,inner,surface,normals,factors,opt);
assert(norm(response(1:48)-eta(1:48),inf) < 1e-13);
assert(norm(response(49:end),inf) > 0);

eta = sin((1:96).');
rhs = cos((1:96).');
coefficients = zeros(72,1);
discarded = zeros(96,1);
for body = 1:2
    rows = (body-1)*48+(1:48);
    columns = (body-1)*36+(1:36);
    coefficients(columns) = applyOneBodyPinv(eta(rows),factors{body});
    leftTranspose = factors{body}.leftTranspose;
    discarded(rows) = eta(rows)- ...
        leftTranspose.'*(leftTranspose*eta(rows));
end
regularizedResidual = matvecDropMFS( ...
    eta,inner,surface,normals,factors,opt)-rhs;
originalResidual = drop_apply_A( ...
    coefficients,inner,surface,normals,factors,opt)-rhs;
assert(norm(originalResidual-(regularizedResidual-discarded),inf) < 2e-11);
fprintf('matvecDropMFS self-test passed.\n');
end
