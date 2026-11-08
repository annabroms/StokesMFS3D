function [Y,ringId,nr,theta] = ring_gl_sphere(radius,centre,Nu,Nv)
%RING_GL_SPHERE Gauss-Legendre-in-cos(theta) sphere with ring membership.
% Wraps setupsurfquad (periodic trapezoid in azimuth on each Gauss-Legendre
% latitude) and returns the nodes together with their ring bookkeeping.
%
%   [Y,ringId,nr,theta] = ring_gl_sphere(radius,centre,Nu,Nv)
%
% Inputs:  radius (scalar), centre (1x3), Nu azimuth-resolution scale,
%          Nv number of latitude rings (rounded up to even by setupsurfquad).
% Outputs: Y      N-by-3 node coordinates, ordered ring by ring (z increasing,
%                 azimuth psi = 2*pi*j/n_r increasing within a ring, j=0 at
%                 the x>0, y=0 meridian). Rows i=1..N are the node index used
%                 by every nodal [x1;y1;z1;x2;...] interleaved vector.
%          ringId N-by-1 ring index r(i) in 1..R.
%          nr     R-by-1 nodes per ring.
%          theta  R-by-1 colatitudes acos(v_r) in (0,pi); poles are never used.
% Cartesian coordinates, z is the polar axis.
if nargin == 0
    self_test();
    return
end
s = setupsurfquad(ellipsoid_param(radius,radius,radius),[Nu,Nv]);
nr = s.Nu(:);
Y = s.x.'+centre(:).';
ringId = repelem((1:numel(nr)).',nr);
theta = acos(s.v(:));
end

function self_test()
%SELF_TEST Check ring bookkeeping, exact radius, latitudes and first meridian.
radius = 2; centre = [0.1 -0.2 0.3];
[Y,ringId,nr,theta] = ring_gl_sphere(radius,centre,12,6);
R = numel(nr);
assert(R==6 && size(Y,1)==sum(nr) && numel(ringId)==sum(nr));
assert(all(nr>=3) && isequal(unique(ringId).',1:R) && issorted(ringId));
assert(max(abs(vecnorm(Y-centre,2,2)-radius))<1e-13);
assert(max(abs(Y(:,3)-centre(3)-radius*cos(theta(ringId))))<1e-13);
first = [1;cumsum(nr(1:end-1))+1];
assert(max(abs(Y(first,2)-centre(2)))<1e-13 && all(Y(first,1)>centre(1)));
assert(all(theta>0 & theta<pi));
fprintf('ring_gl_sphere self-test passed.\n');
end
