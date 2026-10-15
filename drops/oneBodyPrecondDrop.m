function factors = oneBodyPrecondDrop(innerSources, outerSources, ...
    surfacePoints, normals, viscosity, completionScale, velocityScale, ...
    tractionScale, tractionRowWeight, svdTolerance, scaleColumns)
%ONEBODYPRECONDDROP Factor a row-scaled recompleted one-drop block.
% The block maps [alpha(3N);gamma(3N)] to [continuity(3M);
% normal velocity(M);tangential traction(2M)]. The retained factors
% represent the pseudoinverse of the scaled block without forming it.

if nargin == 0
    self_test();
    return
end

[selfMatrix, tangentOne, tangentTwo] = drop_recompleted_matrix( ...
    innerSources, outerSources, surfacePoints, normals, viscosity, ...
    completionScale);
targetCount = size(surfacePoints,1);
selfMatrix(1:4*targetCount,:) = ...
    selfMatrix(1:4*targetCount,:)/velocityScale;
selfMatrix(4*targetCount+1:end,:) = ...
    tractionRowWeight*selfMatrix(4*targetCount+1:end,:)/tractionScale;

if scaleColumns
    columnScales = vecnorm(selfMatrix,2,1).';
    % Zero columns carry no information; a unit scale leaves them zero.
    zeroColumnCount = nnz(columnScales==0);
    columnScales(columnScales==0) = 1;
    matrixForSVD = selfMatrix./columnScales.';
else
    zeroColumnCount = 0;
    columnScales = ones(size(selfMatrix,2),1);
    matrixForSVD = selfMatrix;
end
[leftVectors, singularMatrix, rightVectors] = svd(matrixForSVD,'econ');
singularValues = diag(singularMatrix);
retained = singularValues > svdTolerance*singularValues(1);

factors.leftTranspose = leftVectors(:,retained).';
factors.singularValues = singularValues(retained);
factors.rightVectors = rightVectors(:,retained);
factors.columnScales = columnScales;
factors.zeroColumnCount = zeroColumnCount;
factors.selfMatrix = selfMatrix;
factors.tangentOne = tangentOne;
factors.tangentTwo = tangentTwo;
factors.rank = nnz(retained);
factors.fullColumnCount = size(selfMatrix,2);
factors.svdTolerance = svdTolerance;
factors.scaleColumns = scaleColumns;
factors.viscosity = viscosity;
factors.completionScale = completionScale;
factors.velocityScale = velocityScale;
factors.tractionScale = tractionScale;
factors.tractionRowWeight = tractionRowWeight;

% Translation-invariant geometry is retained to validate factor reuse.
anchor = innerSources(1,:);
factors.innerRelative = innerSources-anchor;
factors.outerRelative = outerSources-anchor;
factors.surfaceRelative = surfacePoints-anchor;
factors.normals = normals;
end

function self_test()
% Check the rectangular dimensions, projection, and reuse metadata.
axisNodes = [eye(3);-eye(3)];
corners = [1 1 1;1 1 -1;1 -1 1;-1 1 1; ...
           1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]/sqrt(3);
factors = oneBodyPrecondDrop(0.6*axisNodes,2*axisNodes, ...
    corners,corners,2,1/(8*pi),1,1,1,1e-12,false);
assert(isequal(size(factors.selfMatrix),[48 36]));
assert(isequal(size(factors.leftTranspose,2),48));
assert(isequal(size(factors.rightVectors,1),36));
assert(factors.rank <= 36 && factors.rank > 0);
eta = sin((1:48).');
projected = factors.leftTranspose.'*(factors.leftTranspose*eta);
assert(norm(factors.leftTranspose.'*factors.leftTranspose*projected- ...
    projected) < 1e-11);
scaledFactors = oneBodyPrecondDrop(0.6*axisNodes,2*axisNodes, ...
    corners,corners,2,1/(8*pi),1,1,1,1e-12,true);
scaledProjection = scaledFactors.leftTranspose.'* ...
    (scaledFactors.leftTranspose*eta);
assert(norm(scaledFactors.selfMatrix* ...
    applyOneBodyPinv(eta,scaledFactors)-scaledProjection) < 1e-11);
fprintf('oneBodyPrecondDrop self-test passed (rank %d of 36).\n', ...
    factors.rank);
end
