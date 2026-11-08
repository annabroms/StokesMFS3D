function velocityXZ = drop_plane_velocity(pointsXZ, interiorBody, ...
    centres, radii, innerSources, outerSources, lambda, gamma, ...
    viscosities, targetBatchSize)
%DROP_PLANE_VELOCITY Evaluate the supplied drop flow on the y=0 plane.
%
% INPUTS
%   pointsXZ       Q-by-2 target coordinates [x,z].
%   interiorBody   0: select the physical region at each target; k>0:
%                  use drop k's interior formula even just outside its
%                  interface. The latter lets ode45 locate that interface.
%   centres,radii  P-by-3 centres and P radii; spheres do not overlap.
%   innerSources,outerSources,lambda,gamma
%                  One cell per body. Strengths are interleaved [x;y;z]
%                  at each source. lambda is the force-free exterior data.
%   viscosities    Interior viscosity ratios; exterior viscosity is one.
%   targetBatchSize Maximum targets included in one Stokeslet block.
%
% OUTPUT
%   velocityXZ     Q-by-2 physical lab-frame velocity components [u_x,u_z].
%
% FIELD FORMULAS
%   Exterior: sum_j S(inner_j)*lambda_j.
%   In drop i: S(outer_i)*gamma_i/mu_i plus the incident exterior fields
%              from all other drops.
% A streamline in the global frame may cross interfaces. This function
% chooses the appropriate formula anew at every ode45 evaluation point.

if nargin == 0
    self_test();
    return
end

bodyCount = numel(innerSources);
pointCount = size(pointsXZ,1);
assert(size(pointsXZ,2)==2 && isequal(size(centres),[bodyCount 3]) && ...
    numel(radii)==bodyCount && numel(viscosities)==bodyCount);
velocityXZ = zeros(pointCount,2);

if interiorBody == 0
    targetRegions = zeros(pointCount,1);
    for body = 1:bodyCount
        distanceSquared = (pointsXZ(:,1)-centres(body,1)).^2 + ...
            centres(body,2)^2 + ...
            (pointsXZ(:,2)-centres(body,3)).^2;
        targetRegions(distanceSquared < radii(body)^2) = body;
    end
else
    targetRegions = interiorBody*ones(pointCount,1);
end

for region = 0:bodyCount
    indices = find(targetRegions==region);
    for first = 1:targetBatchSize:numel(indices)
        selected = indices(first:min(first+targetBatchSize-1,numel(indices)));
        targets = [pointsXZ(selected,1),zeros(numel(selected),1), ...
            pointsXZ(selected,2)];
        velocity = zeros(3*numel(selected),1);
        if region > 0
            velocity = drop_velocity_block( ...
                outerSources{region},targets,viscosities(region))*gamma{region};
        end
        for sourceBody = 1:bodyCount
            if sourceBody ~= region
                velocity = velocity + drop_velocity_block( ...
                    innerSources{sourceBody},targets,1)*lambda{sourceBody};
            end
        end
        values = reshape(velocity,3,[]).';
        velocityXZ(selected,:) = values(:,[1 3]);
    end
end
end

function self_test()
% Compare exterior and both interior formulas with direct Stokeslet sums.
nodes = [eye(3);-eye(3)];
centres = [-2 0 0;2 0 0];
radii = [1;1.3];
viscosities = [2;3];
innerSources = cell(1,2); outerSources = cell(1,2);
lambda = cell(1,2); gamma = cell(1,2);
for body = 1:2
    innerSources{body} = centres(body,:)+0.6*radii(body)*nodes;
    outerSources{body} = centres(body,:)+2*radii(body)*nodes;
    values = reshape(sin((1:18).'+body),3,[]);
    lambda{body} = reshape(values-mean(values,2),[],1);
    gamma{body} = cos((1:18).'+body);
end
pointsXZ = [0 0;0 1.4;-2.3 0.2;2.2 0.1];
actual = drop_plane_velocity(pointsXZ,0,centres,radii, ...
    innerSources,outerSources,lambda,gamma,viscosities,2);
expected = zeros(4,2);
regions = [0;0;1;2];
for point = 1:4
    target = [pointsXZ(point,1) 0 pointsXZ(point,2)];
    region = regions(point);
    field = zeros(3,1);
    if region > 0
        field = drop_velocity_block(outerSources{region},target, ...
            viscosities(region))*gamma{region};
    end
    for body = 1:2
        if body ~= region
            field = field + drop_velocity_block( ...
                innerSources{body},target,1)*lambda{body};
        end
    end
    expected(point,:) = field([1 3]).';
end
assert(isequal(size(actual),[4 2]) && ...
    norm(actual-expected,inf)<1e-13);

% The forced interior formula remains smooth across the physical surface.
outside = [centres(1,1)+1.01*radii(1),0];
forced = drop_plane_velocity(outside,1,centres,radii, ...
    innerSources,outerSources,lambda,gamma,viscosities,1);
target = [outside(1) 0 outside(2)];
reference = drop_velocity_block(outerSources{1},target,viscosities(1))* ...
    gamma{1} + drop_velocity_block(innerSources{2},target,1)*lambda{2};
assert(norm(forced-reference([1 3]).',inf)<1e-13);
fprintf('drop_plane_velocity self-test passed.\n');
end
