function [A,tangentOne,tangentTwo] = drop_two_body_matrix( ...
    innerSources,outerSources,surfacePoints,normals, ...
    interiorViscosity,completionScale)
%DROP_TWO_BODY_MATRIX Assemble the dense two-drop reference matrix.
% Each geometry input is a two-cell array, one entry per drop. Rows and
% columns are body-major: body 1 precedes body 2. Each body's 6M rows are
% velocity continuity (3M), normal velocity (M), and two tangential
% traction values per surface node (2M). Its 6N columns are [alpha;gamma].
% Off-diagonal blocks contain only incident normal-velocity and
% viscosity-contrast-weighted tangential-traction terms.

if nargin == 0
    self_test();
    return
end

assert(iscell(innerSources) && iscell(outerSources) && ...
    iscell(surfacePoints) && iscell(normals) && ...
    numel(innerSources)==2 && numel(outerSources)==2 && ...
    numel(surfacePoints)==2 && numel(normals)==2, ...
    'Geometry inputs must be two-cell arrays.');

tangentOne = cell(1,2);
tangentTwo = cell(1,2);
[selfOne,tangentOne{1},tangentTwo{1}] = drop_recompleted_matrix( ...
    innerSources{1},outerSources{1},surfacePoints{1},normals{1}, ...
    interiorViscosity,completionScale);
[selfTwo,tangentOne{2},tangentTwo{2}] = drop_recompleted_matrix( ...
    innerSources{2},outerSources{2},surfacePoints{2},normals{2}, ...
    interiorViscosity,completionScale);
interactionOneFromTwo = drop_interaction_block( ...
    innerSources{2},surfacePoints{1},normals{1}, ...
    tangentOne{1},tangentTwo{1},interiorViscosity);
interactionTwoFromOne = drop_interaction_block( ...
    innerSources{1},surfacePoints{2},normals{2}, ...
    tangentOne{2},tangentTwo{2},interiorViscosity);
A = [selfOne,interactionOneFromTwo;interactionTwoFromOne,selfTwo];
end

function self_test()
%SELF_TEST Check the body-block order and the shared-field cancellation.
axisNodes = [eye(3);-eye(3)];
corners = [1 1 1;1 1 -1;1 -1 1;-1 1 1; ...
           1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]/sqrt(3);
centres = [-1.6 0 0;1.6 0 0];
inner = cell(1,2);
outer = cell(1,2);
surface = cell(1,2);
normals = cell(1,2);
for body = 1:2
    inner{body} = centres(body,:) + 0.4*axisNodes;
    outer{body} = centres(body,:) + 2*axisNodes;
    surface{body} = centres(body,:) + corners;
    normals{body} = corners;
end
sourceCount = size(axisNodes,1);
targetCount = size(corners,1);
A = drop_two_body_matrix(inner,outer,surface,normals,2,1/(8*pi));
assert(isequal(size(A),[12*targetCount,12*sourceCount]));
assert(~any(A(1:3*targetCount,6*sourceCount+1:end),'all'));
assert(~any(A(6*targetCount+1:9*targetCount,1:6*sourceCount),'all'));
assert(norm(A(3*targetCount+1:4*targetCount, ...
    6*sourceCount+1:9*sourceCount),'fro') > 0);
fprintf('drop_two_body_matrix self-test passed.\n');
end
