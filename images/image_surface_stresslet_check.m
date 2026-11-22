function check = image_surface_stresslet_check( ...
    grid,lambda,fullCount,capCount)
%IMAGE_SURFACE_STRESSLET_CHECK Independently integrate all-family traction.
% Uses outward-from-body normals and nonoverlapping lumped surface areas.
% The solver's paper cap weights are not surface integration weights.

if nargin < 3, fullCount = 2000; end
if nargin < 4, capCount = 1000; end
q = grid.centres;
P = size(q,1);
[unitNodes,~] = get_sphdesign(fullCount);
S = zeros(3,3,P);
F = zeros(6*P,1);
pointCounts = zeros(P,1);
for p = 1:P
    x = q(p,:)+unitNodes;
    for contact = 1:size(grid.body(p).contactDirections,1)
        direction = grid.body(p).contactDirections(contact,:);
        sphere = struct('axes',ones(3,1),'Q',eye(3), ...
            'c',q(p,:).');
        patchOpt = struct('collocationPatchMode',"fibonacci", ...
            'collocationPatchPointCount',capCount, ...
            'collocationPatchAngle',pi/5);
        [cap,~] = collocation_surface_patch( ...
            sphere,q(p,:).'+direction.',patchOpt);
        for k = 1:size(cap,1)
            if ~any(sum((x-cap(k,:)).^2,2)<1e-20)
                x(end+1,:) = cap(k,:); %#ok<AGROW>
            end
        end
    end
    pointCounts(p) = size(x,1);
    n = x-q(p,:);
    n = n./vecnorm(n,2,2);
    sphere = struct('axes',ones(3,1),'Q',eye(3), ...
        'c',q(p,:).');
    weights = collocation_patch_area( ...
        x,ones(size(x,1),1),{sphere});
    t = zeros(size(x));
    for sourceBody = 1:P
        t = t+direct_image_traction( ...
            grid.body(sourceBody).basis,lambda{sourceBody},x,n);
    end
    moment = (x-q(p,:)).'*(weights.*t);
    S(:,:,p) = (moment+moment.')/2-trace(moment)*eye(3)/3;
    F(6*(p-1)+(1:3)) = sum(weights.*t,1).';
    F(6*(p-1)+(4:6)) = sum( ...
        weights.*cross(x-q(p,:),t,2),1).';
end
[sourceStress,~,~,~] = apply_KS_transpose(grid,lambda);
difference = S-sourceStress;
check = struct('tractionStresslet',S,'sourceStresslet',sourceStress, ...
    'relativeDifference',norm(difference(:))/ ...
        max(norm(sourceStress(:)),1e-12), ...
    'tractionForceTorque',F,'pointCounts',pointCounts, ...
    'fullCountRequested',fullCount,'capCountRequested',capCount);
end
