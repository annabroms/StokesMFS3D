function [points,metadata] = image_line_points( ...
    centre,neighborCentre,radius,proxyRadius,N_image,offsetFraction)
%IMAGE_LINE_POINTS Place lubrication-adapted sources inside one sphere.
%   [points,metadata] = image_line_points(c1,c2,R,Rp,N_image,offsetFraction)
%   image_line_points() runs a small geometry check and illustration.
%
% Inputs: centres are 1-by-3 Cartesian rows; R and Rp are physical radii,
% 0<Rp<R. N_image is the number of nodes on this pair's line, including the
% accumulation point. offsetFraction defaults to 0.05. The sphere gap is
% norm(c2-c1)-2R and must be positive. Outputs are N_image-by-3 locations ordered
% from the accumulation point towards the proxy sphere, and metadata with
% gap, unit line direction, accumulation radius and radial coordinates.
% Uses Eqs. (24), (26), and (27) of A. Broms, A. H. Barnett, and
% A.-K. Tornberg, "Accurate close interactions of Stokes spheres using
% lubrication-adapted image systems," J. Comput. Phys. 523 (2025), 113636.

if nargin == 0
    self_test();
    points = zeros(0,3);
    metadata = struct();
    return
end
if nargin < 6 || isempty(offsetFraction)
    offsetFraction = 0.05;
end
validateattributes(centre,{'numeric'},{'real','finite','size',[1,3]});
validateattributes(neighborCentre,{'numeric'},{'real','finite','size',[1,3]});
validateattributes(radius,{'numeric'},{'scalar','positive','finite'});
validateattributes(proxyRadius,{'numeric'},{'scalar','positive','<',radius});
validateattributes(N_image,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(offsetFraction,{'numeric'},{'scalar','nonnegative','finite'});
separation = norm(neighborCentre-centre);
gap = separation-2*radius;
assert(gap > 0,'image_line_points:OverlappingSpheres');
direction = (neighborCentre-centre)/separation;
relativeGap = gap/radius;
accumulationRadius = radius*(1+relativeGap/2- ...
    sqrt(relativeGap+relativeGap^2/4));
if N_image > 0
    innerRadius = proxyRadius*(1+offsetFraction);
    assert(innerRadius < accumulationRadius, ...
        'image_line_points:NoRoomBetweenProxyAndAccumulation');
    index = (0:N_image-1).';
    radialCoordinates = innerRadius+(accumulationRadius-innerRadius)* ...
        cos(pi*index/(2*N_image));
    points = centre+radialCoordinates*direction;
else
    radialCoordinates = zeros(0,1);
    points = zeros(0,3);
end
metadata = struct('gap',gap,'direction',direction, ...
    'accumulationRadius',accumulationRadius, ...
    'radialCoordinates',radialCoordinates);
end

function self_test()
%SELF_TEST checks the endpoints, rotation and interior membership.
centre = [0,0,0];
neighbor = [2.01,0,0];
[points,info] = image_line_points(centre,neighbor,1,0.63,20,0.05);
assert(size(points,1) == 20);
assert(abs(points(1,1)-info.accumulationRadius) < 1e-14);
assert(all(diff(info.radialCoordinates) < 0));
assert(all(vecnorm(points-centre,2,2) < 1));
[rotated,~] = image_line_points(centre,[0,2.01,0],1,0.63,20,0.05);
assert(norm(rotated(:,2)-points(:,1),inf) < 1e-14);
figure('Name','Image line geometry self-test');
plot(points(:,1),zeros(size(points,1),1),'ro'); hold on;
xline(0.63,'b--','proxy radius');
xline(info.accumulationRadius,'k--','accumulation point');
xlabel('Distance from sphere centre'); ylabel('Transverse coordinate');
title('Half-Chebyshev image line towards the gap'); grid on;
fprintf('image_line_points self-test passed.\n');
end
