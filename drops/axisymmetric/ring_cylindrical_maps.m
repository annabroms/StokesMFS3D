function [C,W,KR] = ring_cylindrical_maps(ringId,psi)
%RING_CYLINDRICAL_MAPS Axisymmetric ring map: radial + axial amplitude per ring.
%
%   [C,W,KR] = ring_cylindrical_maps(ringId,psi)
%
% Node j of ring r carries lambda_j = a_rho,r*(cos psi_j, sin psi_j, 0) + a_z,r*(0,0,1).
% Inputs:  ringId N-by-1 ring index (1..R); psi N-by-1 azimuth of each node (rad,
%          about the z axis through the sphere centre).
% Outputs: C  3N-by-2R sparse, lambda = C*b with b = [a_rho,1;a_z,1;...;a_rho,R;a_z,R]
%             and lambda nodal-interleaved [x1;y1;z1;x2;...];
%          W  2R-by-1, W'*b = total axial force sum_r n_r a_z,r (the radial
%             parts cancel by azimuthal symmetry and give no net force);
%          KR 2R-by-1 (a_z entries 1), so W'*KR = N.
% Same roles as in ring_constant_maps, with the z force as the only constraint.
if nargin == 0
    self_test();
    return
end
ringId = ringId(:); psi = psi(:);
N = numel(ringId);
R = max(ringId);
nr = accumarray(ringId,1);
j = (1:N).';
rows = [3*j-2; 3*j-1; 3*j];
cols = [2*ringId-1; 2*ringId-1; 2*ringId];
vals = [cos(psi); sin(psi); ones(N,1)];
C = sparse(rows,cols,vals,3*N,2*R);
W = zeros(2*R,1);  W(2:2:end) = nr;
KR = zeros(2*R,1); KR(2:2:end) = 1;
end

function self_test()
%SELF_TEST Radial parts rotate with psi and carry no net force.
ringId = [1;1;1;1;2;2;2;2];
psi = [0;pi/2;pi;3*pi/2;0;pi/2;pi;3*pi/2];
[C,W,KR] = ring_cylindrical_maps(ringId,psi);
b = [2;5;-1;3];                          % [a_rho1;a_z1;a_rho2;a_z2]
lam = reshape(C*b,3,[]);
assert(norm(lam(:,2)-[0;2;5])<1e-14 && norm(lam(:,7)-[1;0;3])<1e-14);
F = sum(lam,2);
assert(norm(F-[0;0;4*5+4*3])<1e-13 && abs(W.'*b-F(3))<1e-13);
assert(W.'*KR==8);
fprintf('ring_cylindrical_maps self-test passed.\n');
end
