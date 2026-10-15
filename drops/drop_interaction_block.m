function Ainteraction = drop_interaction_block( ...
    otherInnerSources, targetPoints, targetNormals, ...
    tangentOne, tangentTwo, targetViscosity)
%DROP_INTERACTION_BLOCK Recompleted target-from-other-drop MFS block.
% Columns are the other drop's [alpha;gamma], three components per node.
% Rows are target-drop [continuity(3M); normal velocity(M);
% tangent traction(2M, interleaved)]. With the shared incident field,
% continuity and all gamma columns vanish. The remaining rows are
% N*G*P and (1-targetViscosity)*D*T*P, respectively.

if nargin == 0
    self_test();
    return
end

validateattributes(otherInnerSources, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'}, mfilename, 'otherInnerSources');
validateattributes(targetPoints, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'}, mfilename, 'targetPoints');
assert(isequal(size(targetNormals),size(targetPoints)) && ...
    isequal(size(tangentOne),size(targetPoints)) && ...
    isequal(size(tangentTwo),size(targetPoints)), ...
    'Each target needs a normal and two tangent directions.');
validateattributes(targetViscosity, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'targetViscosity');

sourceCount = size(otherInnerSources,1);
targetCount = size(targetPoints,1);
G = drop_velocity_block(otherInnerSources,targetPoints,1);
T = drop_traction_block(otherInnerSources,targetPoints,targetNormals);
constantSourceColumns = kron(ones(sourceCount,1),eye(3));
forceFreeVelocity = G - ...
    (G*constantSourceColumns)*constantSourceColumns.'/sourceCount;
forceFreeTraction = T - ...
    (T*constantSourceColumns)*constantSourceColumns.'/sourceCount;

normalRows = zeros(targetCount,3*sourceCount);
tangentRows = zeros(2*targetCount,3*sourceCount);
for node = 1:targetCount
    targetRows = 3*node-2:3*node;
    normalRows(node,:) = targetNormals(node,:)*forceFreeVelocity(targetRows,:);
    tangentRows(2*node-1:2*node,:) = ...
        [tangentOne(node,:);tangentTwo(node,:)] * ...
        forceFreeTraction(targetRows,:);
end
Ainteraction = [zeros(3*targetCount,6*sourceCount); ...
    normalRows, zeros(targetCount,3*sourceCount); ...
    (1-targetViscosity)*tangentRows, ...
        zeros(2*targetCount,3*sourceCount)];
assert(isequal(size(Ainteraction),[6*targetCount,6*sourceCount]));
end

function self_test()
%SELF_TEST Check block layout, force projection, and viscosity coupling.
axisNodes = [eye(3);-eye(3)];
corners = [1 1 1;1 1 -1;1 -1 1;-1 1 1; ...
           1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]/sqrt(3);
targetCount = size(corners,1);
sourceCount = size(axisNodes,1);
tangentOne = cross(corners,repmat([0 0 1],targetCount,1),2);
tangentOne = tangentOne./vecnorm(tangentOne,2,2);
tangentTwo = cross(corners,tangentOne,2);
otherSources = [3.2 0 0] + 0.4*axisNodes;
block = drop_interaction_block(otherSources,corners,corners, ...
    tangentOne,tangentTwo,2);
assert(isequal(size(block),[6*targetCount,6*sourceCount]));
assert(~any(block(1:3*targetCount,:),'all'));
assert(~any(block(:,3*sourceCount+1:end),'all'));
assert(norm(block(4*targetCount+1:end,:),'fro') > 0);
constantAlpha = repmat([1;-2;3],sourceCount,1);
assert(norm(block*[constantAlpha;zeros(3*sourceCount,1)],inf) < 1e-12);

forceFreeAlpha = zeros(3*sourceCount,1);
forceFreeAlpha(1) = 1;
forceFreeAlpha(4) = -1;
response = block*[forceFreeAlpha;zeros(3*sourceCount,1)];
directVelocity = reshape( ...
    drop_velocity_block(otherSources,corners,1)*forceFreeAlpha,3,[]).';
directTraction = reshape( ...
    drop_traction_block(otherSources,corners,corners)*forceFreeAlpha,3,[]).';
expectedNormal = sum(corners.*directVelocity,2);
expectedTangents = reshape([ ...
    -sum(tangentOne.*directTraction,2), ...
    -sum(tangentTwo.*directTraction,2)].',[],1);
assert(norm(response(3*targetCount+1:4*targetCount)- ...
    expectedNormal,inf) < 1e-12);
assert(norm(response(4*targetCount+1:end)-expectedTangents,inf) < 1e-12);

matchingViscosityBlock = drop_interaction_block( ...
    otherSources,corners,corners,tangentOne,tangentTwo,1);
assert(~any(matchingViscosityBlock(4*targetCount+1:end,:),'all'));
assert(norm(matchingViscosityBlock(3*targetCount+1:4*targetCount,:),'fro') > 0);
fprintf('drop_interaction_block self-test passed.\n');
end
