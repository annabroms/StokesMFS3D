function [A, tangentOne, tangentTwo] = drop_recompleted_matrix( ...
    innerSources, outerSources, surfacePoints, normals, ...
    interiorViscosity, completionScale)
%DROP_RECOMPLETED_MATRIX Assemble one drop's dense boundary-condition block.
% A is 6M-by-6N, with columns [alpha; gamma]. Each source node stores
% [x;y;z]. Rows are 3M velocity-continuity components, M normal-velocity
% residuals, then two tangential-traction residuals per surface node.
% alpha is projected to the force-free exterior strength
% lambda=(I-K*K'/N)*alpha, K=kron(ones(N,1),eye(3)); its removed mean
% determines U=-completionScale*K'*alpha. The normal row includes the
% residual contribution from this translation. gamma generates the
% interior field; its velocity is divided by interiorViscosity. Traction
% uses the outward normal and is not divided by viscosity here.

if nargin == 0
    self_test();
    return
end

validateattributes(innerSources, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'}, mfilename, 'innerSources');
assert(isequal(size(innerSources), size(outerSources)), ...
    'Inner and outer proxy surfaces must use the same source count.');
assert(isequal(size(surfacePoints), size(normals)), ...
    'There must be one normal per surface point.');
validateattributes(completionScale, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'completionScale');

sourceCount = size(innerSources,1);
targetCount = size(surfacePoints,1);
referenceDirections = repmat([1 0 0], targetCount, 1);
nearXAxis = abs(normals(:,1)) > 0.9;
referenceDirections(nearXAxis,:) = repmat([0 1 0], sum(nearXAxis), 1);
tangentOne = cross(normals, referenceDirections, 2);
tangentOne = tangentOne ./ vecnorm(tangentOne, 2, 2);
tangentTwo = cross(normals, tangentOne, 2);

normalExtraction = zeros(targetCount, 3*targetCount);
tangentExtraction = zeros(2*targetCount, 3*targetCount);
for node = 1:targetCount
    columns = 3*node-2:3*node;
    normalExtraction(node,columns) = normals(node,:);
    tangentExtraction(2*node-1:2*node,columns) = ...
        [tangentOne(node,:); tangentTwo(node,:)];
end

G = drop_velocity_block(innerSources, surfacePoints, 1);
E = drop_velocity_block(outerSources, surfacePoints, interiorViscosity);
T = drop_traction_block(innerSources, surfacePoints, normals);
Q = drop_traction_block(outerSources, surfacePoints, normals);

constantSourceColumns = kron(ones(sourceCount,1), eye(3));
% Apply P=I-K*K'/N as a rank-three correction, avoiding a dense projector.
forceFreeVelocity = G - ...
    (G*constantSourceColumns)*constantSourceColumns.'/sourceCount;
forceFreeTraction = T - ...
    (T*constantSourceColumns)*constantSourceColumns.'/sourceCount;
translationColumns = completionScale * normals * constantSourceColumns.';

A = [forceFreeVelocity, -E; ...
     normalExtraction*forceFreeVelocity + translationColumns, ...
         zeros(targetCount, 3*sourceCount); ...
     tangentExtraction*forceFreeTraction, -tangentExtraction*Q];
end

function self_test()
%SELF_TEST Check block sizes and the force-free translation completion.
axisNodes = [eye(3); -eye(3)];
cornerNodes = [1 1 1; 1 1 -1; 1 -1 1; -1 1 1; ...
               1 -1 -1; -1 1 -1; -1 -1 1; -1 -1 -1] / sqrt(3);
sourceCount = size(axisNodes,1);
targetCount = size(cornerNodes,1);
completionScale = 0.3;
[A, tangentOne, tangentTwo] = drop_recompleted_matrix( ...
    0.68*axisNodes, 1.3*axisNodes, cornerNodes, cornerNodes, ...
    2, completionScale);
assert(isequal(size(A), [6*targetCount, 6*sourceCount]));
assert(max(abs(sum(tangentOne.*cornerNodes,2))) < 1e-14);
assert(max(abs(sum(tangentTwo.*cornerNodes,2))) < 1e-14);

constantCoefficient = [1; -2; 3];
alpha = repmat(constantCoefficient, sourceCount, 1);
response = A * [alpha; zeros(3*sourceCount,1)];
assert(norm(response(1:3*targetCount), inf) < 1e-12);
assert(norm(response(4*targetCount+1:end), inf) < 1e-12);
expectedNormal = completionScale*sourceCount * ...
    (cornerNodes*constantCoefficient);
assert(norm(response(3*targetCount+1:4*targetCount) - ...
    expectedNormal, inf) < 1e-12);
fprintf('drop_recompleted_matrix self-test passed.\n');
end
