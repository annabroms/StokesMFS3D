function [points,meta] = exterior_image_line_points( ...
    centre,direction,radius,outerRadius,accumulationRadius,N,offset)
%EXTERIOR_IMAGE_LINE_POINTS Kelvin-mirror an interior accumulation point.
% Nodes start at a^2/r* and cluster there on a line toward the outer proxy.
if nargin == 0
    self_test(); points = zeros(0,3); meta = struct(); return
end
validateattributes(centre,{'numeric'},{'size',[1,3],'finite','real'});
validateattributes(direction,{'numeric'},{'size',[1,3],'finite','real'});
validateattributes(radius,{'numeric'},{'scalar','positive','finite'});
validateattributes(outerRadius,{'numeric'},{'scalar','>',radius,'finite'});
validateattributes(N,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(offset,{'numeric'},{'scalar','>',0,'<',1});
direction = direction/norm(direction);
reflected = radius^2/accumulationRadius;
endpoint = outerRadius*(1-offset);
if N > 0 && ~(radius < reflected && reflected < endpoint && endpoint < outerRadius)
    error('exterior_image_line_points:NoRoom', ...
        'Kelvin point and endpoint must lie strictly between boundary and outer proxy.');
end
if N == 0
    r = zeros(0,1);
elseif N == 1
    r = reflected;
else
    k = (0:N-1).';
    r = endpoint-(endpoint-reflected)*cos(pi*k/(2*(N-1)));
end
points = centre+r*direction;
meta.reflectedRadius = reflected;
meta.endpointRadius = endpoint;
meta.radialCoordinates = r;
end

function self_test()
% Check Kelvin inversion, node count, and both strict radial bounds.
[x,m] = exterior_image_line_points([0 0 0],[1 0 0],1,2,0.8,8,0.05);
assert(size(x,1)==8 && abs(x(1,1)-1.25)<1e-14);
assert(all(diff(m.radialCoordinates)>0) && all(x(:,1)>1) && all(x(:,1)<2));
assert(abs(x(end,1)-1.9)<1e-14);
[x,m] = exterior_image_line_points([0 0 0],[1 0 0],1,2,0.8,1,0.01);
assert(size(x,1)==1 && abs(x(1,1)-m.reflectedRadius)<1e-14);
fprintf('exterior_image_line_points self-test passed.\n');
end
