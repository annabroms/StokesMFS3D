function directions = image_stresslet_directions(lineDirections,use_orth_dir_images,family)
%IMAGE_STRESSLET_DIRECTIONS Choose line-parallel or transverse source normals.
% T_n is normal to the line; T_t is tangent to it. Plain T uses the boolean.
if nargin == 0
    self_test(); directions = zeros(0,3); return
end
validateattributes(lineDirections,{'numeric'}, ...
    {'real','finite','2d','ncols',3});
validateattributes(use_orth_dir_images,{'logical'},{'scalar'});
if nargin < 3, family = "T"; end
switch string(family)
    case "T_n"
        use_orth_dir_images = true;
    case "T_t"
        use_orth_dir_images = false;
    case "T"
        % Preserve the original single-direction interface.
    otherwise
        error('image_stresslet_directions:BadFamily','Expected T, T_n or T_t.');
end
lengths = vecnorm(lineDirections,2,2);
if any(lengths==0)
    error('image_stresslet_directions:ZeroRay', ...
        'Every image line must have a nonzero direction.');
end
ray = lineDirections./lengths;
directions = ray;
if ~use_orth_dir_images || isempty(ray), return; end

% Cross each ray with the least-aligned Cartesian axis; opposite rays
% then receive opposite, deterministic transverse directions.
[~,axisIndex] = min(abs(ray),[],2);
reference = zeros(size(ray));
rows = (1:size(ray,1)).';
reference(sub2ind(size(ray),rows,axisIndex)) = 1;
directions = cross(ray,reference,2);
directions = directions./vecnorm(directions,2,2);
end

function self_test()
% Check orthogonality, unit length, pair reversal, and the legacy choice.
ray = [1 0 0;-1 0 0;0 0 1;1 2 3];
parallel = image_stresslet_directions(ray,false);
orthogonal = image_stresslet_directions(ray,true);
if ~isequal(orthogonal,image_stresslet_directions(ray,false,"T_n")) || ...
        ~isequal(parallel,image_stresslet_directions(ray,true,"T_t"))
    error('image_stresslet_directions:SelfTest','Explicit directions changed.');
end
if norm(parallel-ray./vecnorm(ray,2,2),'fro')>1e-14 || ...
        max(abs(sum(orthogonal.*ray,2)))>1e-14 || ...
        max(abs(vecnorm(orthogonal,2,2)-1))>1e-14 || ...
        norm(orthogonal(1,:)-[0 0 1])>1e-14 || ...
        norm(orthogonal(2,:)+orthogonal(1,:))>1e-14
    error('image_stresslet_directions:SelfTest', ...
        'The stresslet orientations do not match their specified rule.');
end
fprintf('image_stresslet_directions self-test passed.\n');
end
