function A = lap3d_charge_potential_matrix(targets,sources,sourceNormals)
%LAP3D_CHARGE_POTENTIAL_MATRIX Dense 3D Laplace charge-to-potential matrix.
%
%   A = lap3d_charge_potential_matrix(targets,sources)
%
% targets is M-by-3, sources is N-by-3, and A is M-by-N with
% A(i,j) = 1/(4*pi*norm(targets(i,:)-sources(j,:))).
% With a third input of unit sourceNormals (N-by-3), instead return DLP:
% dG/dn_y = n_y dot (x-y)/(4*pi*|x-y|^3), the SOURCE normal derivative.
% These are point monopole/dipole strengths, not quadrature-weighted layer
% densities.  Targets lie off the proxy surface, so no jump term is added.
%
% The normalization and target-from-source convention follow MFS3D.jl in
% https://github.com/ahbarnett/mfs-mobility

if nargin == 0
    self_test;
    A = [];
    return
end

if size(targets,2) ~= 3 || size(sources,2) ~= 3
    error('lap3d_charge_potential_matrix:BadCoordinates', ...
        'targets and sources must have three columns.');
end
if ~isreal(targets) || ~isreal(sources) || ...
        any(~isfinite(targets),'all') || any(~isfinite(sources),'all')
    error('lap3d_charge_potential_matrix:BadCoordinates', ...
        'targets and sources must contain finite real coordinates.');
end

dx = targets(:,1) - sources(:,1).';
dy = targets(:,2) - sources(:,2).';
dz = targets(:,3) - sources(:,3).';
distanceSquared = dx.^2 + dy.^2 + dz.^2;

if any(distanceSquared == 0,'all')
    error('lap3d_charge_potential_matrix:CoincidentPoints', ...
        'A target coincides with a source.');
end

if nargin < 3
    A = 1./(4*pi*sqrt(distanceSquared));
else
    if ~isequal(size(sourceNormals),size(sources)) || ...
            ~isreal(sourceNormals) || any(~isfinite(sourceNormals),'all') || ...
            any(abs(vecnorm(sourceNormals,2,2)-1) > 1e-10)
        error('lap3d_charge_potential_matrix:BadSourceNormals', ...
            'sourceNormals must contain one finite real unit normal per source.');
    end
    normalDisplacement = dx.*sourceNormals(:,1).' + ...
        dy.*sourceNormals(:,2).' + dz.*sourceNormals(:,3).';
    A = normalDisplacement./(4*pi*distanceSquared.*sqrt(distanceSquared));
end
end


function self_test()
%SELF_TEST checks the normalization and target/source ordering.

targets = [2,0,0; 0,3,4];
sources = [0,0,0];
A = lap3d_charge_potential_matrix(targets,sources);
expected = [1/(8*pi); 1/(20*pi)];
assert(norm(A-expected,inf) < 20*eps);
sourceNormal = [1,0,0];
D = lap3d_charge_potential_matrix(targets,sources,sourceNormal);
assert(norm(D-[1/(16*pi);0],inf) < 20*eps);
h = 1e-5;
finiteDifference = (lap3d_charge_potential_matrix(targets,sources+h*sourceNormal)- ...
    lap3d_charge_potential_matrix(targets,sources-h*sourceNormal))/(2*h);
assert(norm(D-finiteDifference,inf) < 1e-10);
fprintf('lap3d_charge_potential_matrix self-test passed.\n');
end
