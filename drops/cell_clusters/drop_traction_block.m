function H = drop_traction_block(sources, targets, normals)
%DROP_TRACTION_BLOCK Stokeslet traction with the drop-to-fluid normal.
% The 3M-by-3N block maps interleaved source forces to interleaved
% tractions sigma*n at targets. For r=x-y, H=-3(r*n)rr'/(4*pi*|r|^5).
% This is the sign in the supplied drop formulation. Traction for a fixed
% Stokeslet strength is independent of the fluid viscosity.

if nargin == 0
    self_test();
    return
end

validateattributes(sources, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'}, mfilename, 'sources');
validateattributes(targets, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'}, mfilename, 'targets');
validateattributes(normals, {'numeric'}, ...
    {'real','finite','2d','ncols',3}, mfilename, 'normals');
assert(isequal(size(normals), size(targets)), ...
    'There must be one normal per target.');
assert(all(abs(vecnorm(normals, 2, 2) - 1) < 1e-12), ...
    'Target normals must have unit length.');

rx = targets(:,1) - sources(:,1).';
ry = targets(:,2) - sources(:,2).';
rz = targets(:,3) - sources(:,3).';
distanceSquared = rx.^2 + ry.^2 + rz.^2;
assert(all(distanceSquared > 0, 'all'), ...
    'A source coincides with a target.');
radialNormal = rx.*normals(:,1) + ry.*normals(:,2) + ...
    rz.*normals(:,3);
scale = -3/(4*pi) * radialNormal ./ distanceSquared.^(5/2);

targetCount = size(targets,1);
sourceCount = size(sources,1);
H = zeros(3*targetCount, 3*sourceCount);
xRows = 1:3:3*targetCount;
yRows = 2:3:3*targetCount;
zRows = 3:3:3*targetCount;
xColumns = 1:3:3*sourceCount;
yColumns = 2:3:3*sourceCount;
zColumns = 3:3:3*sourceCount;
H(xRows,xColumns) = scale.*rx.*rx;
H(xRows,yColumns) = scale.*rx.*ry;
H(xRows,zColumns) = scale.*rx.*rz;
H(yRows,xColumns) = scale.*ry.*rx;
H(yRows,yColumns) = scale.*ry.*ry;
H(yRows,zColumns) = scale.*ry.*rz;
H(zRows,xColumns) = scale.*rz.*rx;
H(zRows,yColumns) = scale.*rz.*ry;
H(zRows,zColumns) = scale.*rz.*rz;
end

function self_test()
%SELF_TEST Check dimensions, sign, and a tangential reference value.
H = drop_traction_block([0 0 0], [1 0 0; 0 1 0], ...
    [1 0 0; 0 1 0]);
assert(isequal(size(H), [6 3]));
assert(norm(H(1:3,:) - diag([-3/(4*pi), 0, 0]), 'fro') < 1e-14);
assert(norm(H(4:6,:) - diag([0, -3/(4*pi), 0]), 'fro') < 1e-14);
assert(norm(H(1:3,:) * [0; 1; 0]) < 1e-14);
fprintf('drop_traction_block self-test passed.\n');
end
