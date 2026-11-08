function [points,normals] = drop_image_validation_grid(grid,count,gapCount,angle)
%DROP_IMAGE_VALIDATION_GRID Independent rotated design plus gap-facing caps.
% angle and gapCount may be matching vectors for multiple concentric caps.
if nargin == 0
    self_test(); points={}; normals={}; return
end
angle = angle(:).'; gapCount = gapCount(:).';
if isscalar(gapCount), gapCount = repmat(gapCount,size(angle)); end
if numel(gapCount)~=numel(angle)
    error('drop_image_validation_grid:CapCounts','Supply one count per cap angle.');
end
[v,~] = get_sphdesign(count);
cz = cos(0.37); sz = sin(0.37);
cy = cos(0.23); sy = sin(0.23);
R = [cz -sz 0;sz cz 0;0 0 1]*[cy 0 sy;0 1 0;-sy 0 cy];
v = (R*v.').';
P = numel(grid.body);
points = cell(P,1); normals = cell(P,1);
for i = 1:P
    q = grid.centres(i,:); a = grid.radii(i);
    x = q+a*v;
    for k = 1:size(grid.body(i).contactDirections,1)
        d = grid.body(i).contactDirections(k,:);
        sphere = struct('axes',a*ones(3,1),'Q',eye(3),'c',q.');
        for capIndex = 1:numel(angle)
            capOpt = struct('collocationPatchMode',"fibonacci", ...
                'collocationPatchPointCount',gapCount(capIndex), ...
                'collocationPatchAngle',angle(capIndex));
            [cap,~] = collocation_surface_patch(sphere,q.'+a*d.',capOpt);
            x = [x;cap]; %#ok<AGROW>
        end
    end
    points{i} = x;
    normals{i} = (x-q)/a;
end
end

function self_test()
% Check independent surface placement for unequal radii.
grid.centres = [0 0 0;3.1 0 0]; grid.radii = [1;2];
grid.body(1).contactDirections = [1 0 0];
grid.body(2).contactDirections = [-1 0 0];
[x,n] = drop_image_validation_grid(grid,50,21,pi/5);
assert(all(cellfun(@(v) size(v,1)>50,x)));
for i = 1:2
    assert(max(abs(vecnorm(x{i}-grid.centres(i,:),2,2)-grid.radii(i)))<1e-12);
    assert(max(abs(vecnorm(n{i},2,2)-1))<1e-12);
end
[xx,nn] = drop_image_validation_grid(grid,50,[21 31],[pi/5 pi/60]);
for i = 1:2
    if size(xx{i},1)~=size(x{i},1)+31 || ...
            ~isequal(xx{i}(1:size(x{i},1),:),x{i}) || ...
            max(abs(vecnorm(nn{i},2,2)-1))>1e-12
        error('drop_image_validation_grid:SelfTest','Multiple validation caps failed.');
    end
end
fprintf('drop_image_validation_grid self-test passed.\n');
end
