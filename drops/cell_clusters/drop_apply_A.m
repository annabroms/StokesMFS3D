function response = drop_apply_A(coefficients, innerSources, ...
    surfacePoints, normals, factors, opt)
%DROP_APPLY_A Apply the row-scaled physical coefficient operator.
% Maps body-major [alpha_i;gamma_i] coefficients to continuity,
% normal-velocity, and tangential-traction residuals. Cached one-body
% matrices supply the diagonal blocks; direct cross-drop sums supply the
% off-diagonal blocks. This is the original operator, not the identity-
% diagonal boundary-data operator applied by matvecDropMFS.

if nargin == 0
    self_test();
    return
end

bodyCount = numel(innerSources);
assert(numel(surfacePoints)==bodyCount && numel(normals)==bodyCount && ...
    numel(factors)==bodyCount && iscolumn(coefficients));
sourceCounts = cellfun(@(points) size(points,1),innerSources);
targetCounts = cellfun(@(points) size(points,1),surfacePoints);
assert(numel(coefficients)==6*sum(sourceCounts));

response = zeros(6*sum(targetCounts),1);
lambda = cell(1,bodyCount);
columnOffset = 0;
rowOffset = 0;
for body = 1:bodyCount
    sourceCount = sourceCounts(body);
    targetCount = targetCounts(body);
    localCoefficients = coefficients(columnOffset+(1:6*sourceCount));
    assert(isequal(size(factors{body}.selfMatrix), ...
        [6*targetCount,6*sourceCount]));
    response(rowOffset+(1:6*targetCount)) = ...
        factors{body}.selfMatrix*localCoefficients;
    alphaNodes = reshape(localCoefficients(1:3*sourceCount),3,sourceCount);
    lambda{body} = reshape(alphaNodes-mean(alphaNodes,2),[],1);
    columnOffset = columnOffset+6*sourceCount;
    rowOffset = rowOffset+6*targetCount;
end
response = response+drop_apply_interactions( ...
    lambda,innerSources,surfacePoints,normals,factors,opt);
end

function self_test()
% Match a small explicitly assembled and row-scaled two-drop matrix.
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
matrix = drop_two_body_matrix(inner,outer,surface,normals,2,1/(8*pi));
for body = 1:2
    offset = (body-1)*48;
    matrix(offset+(1:32),:) = matrix(offset+(1:32),:)/3;
    matrix(offset+33:offset+48,:) = ...
        1.7*matrix(offset+33:offset+48,:)/5;
end
coefficients = sin((1:72).');
actual = drop_apply_A(coefficients,inner,surface,normals,factors,opt);
assert(isequal(size(actual),[96 1]));
assert(norm(actual-matrix*coefficients,inf) < 2e-11);
fprintf('drop_apply_A self-test passed.\n');
end
