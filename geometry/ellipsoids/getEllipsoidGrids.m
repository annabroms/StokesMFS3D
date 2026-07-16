function [rin, rout, qvec, R, E, w, n_out, n_in] = getEllipsoidGrids(E0, P, delta, N1, N2, sep, R, qvec)
%GETELLIPSOIDGRIDS Generate MFS source and collocation grids on multiple ellipsoids.
%
%   [rin, rout, qvec, R, E, w, n_out, n_in] = GETELLIPSOIDGRIDS(E0, P, delta, N1, N2, sep, R, qvec)
%
%   Constructs source (proxy) points `rin` and collocation points `rout` on the
%   boundaries of `P` ellipsoidal particles, with prescribed geometric and spatial
%   configuration. The function supports ellipsoids of equal shape, oriented and
%   placed in 3D space, with minimal separation control.
%
%   INPUTS:
%       E0    - Vector with 3 components that sets the semi axis of the
%               ellipsoid
%       P     - Number of ellipsoids to generate.
%       delta - Minimum surface-to-surface separation between ellipsoids.
%               The cluster of particles is generated such that all
%               ellipsoids have at least one close neighbour
%       N1    - Number of uniform nodes in the s direction
%       N2    - Number of GL nodes in the t direction
%       sep   - Separation between prxoy and collocation surface, projected
%       in the normal direction 
%       R     - (Optional) 1 x P cell array of 3x3 rotation matrices for orientation.
%               If omitted or empty, positions and orientations are
%               generated at random.
%       qvec  - (Optional) P x 3 matrix of translation vectors for each
%               ellipsoid. A 1 x P cell array of 3D vectors is also accepted.
%               If omitted or empty, positions and orientations are
%               generated at random.
%
%   OUTPUTS:
%       rin   - concatenated source points for all particles.
%       rout  - concatenated collocation points on all surfaces.
%       qvec  - P x 3 matrix of translation vectors (used or generated).
%       R     - 1 x P cell array of 3x3 rotation matrices (used or generated).
%       E     - 1 x P cell array of aspect ratios for the P particles.
%       w     - PM x 1 vector of collocation quadrature weights, formed by
%               stacking one copy of b_outer.w for each ellipsoid.
%       n_out - PM x 3 array of collocation normals, stacked in the same
%               order as rout.
%       n_in  - PN x 3 array of inner proxy normals, stacked in the same
%               order as rin.
%
%   NOTES:
%       - The generated `rin` and `rout` can be used in MFS solvers for
%         mobility/resistance problems involving ellipsoidal particles.
%       - Each ellipsoid is assumed to have the same shape but may differ in
%         position and orientation.
%       - Both the proxy and collocation surfaces are discretized using the
%         quasi-uniform ellipsoid grid described by Stein and Barnett
%         in [Stein2022: Quadrature by Fundamental solutions]. Each ellipsoid 
%         is parameterized in Cartesian coordinates as
%             (a*sqrt(1 - t^2)*cos(s), b*sqrt(1 - t^2)*sin(s), c*t),
%         where (s, t) ∈ [0, 2π] × [-1, 1]. The parameter t is discretized
%         using Gauss–Legendre nodes, while s uses a
%         periodic trapezoidal rule.
%
%
%  

if nargin == 0
    self_test;
    return;
end

%Create source points
if nargin < 7 || isempty(R) || nargin < 8 || isempty(qvec)
    [E, R, t] = ellipsoid_cluster(E0,P,delta);
    qvec = cell_centers_to_matrix(t);
else
    E = cell(1,P);
    if iscell(qvec)
        qvec = cell_centers_to_matrix(qvec);
    end
end

if size(qvec,1) ~= P || size(qvec,2) ~= 3
    error('getEllipsoidGrids:qvec_size', 'qvec must be of size P x 3.');
end

%Create proxy grid
b_inner = ellipsoid_param(E0(1),E0(2),E0(3));   % baseline object at the origin, aligned
b_inner = setupsurfquad(b_inner,[ceil(N1),ceil(N2)]);

%Create grid of collocation points
b_outer = ellipsoid_param(E0(1),E0(2),E0(3));   % baseline object at the oridin, aligned
b_outer = setupsurfquad(b_outer,[ceil(N1*1.15),ceil(N2*1.15)]);
w_i = b_outer.w;

Nin = size(b_inner.x,2);
Nout = size(b_outer.x,2);
rin = zeros(P*Nin,3);
rout = zeros(P*Nout,3);
w = repmat(w_i(:),P,1);
n_out = zeros(P*Nout,3);
n_in = zeros(P*Nin,3);

for k = 1:P
    ii_in = (k-1)*Nin + (1:Nin);
    ii_out = (k-1)*Nout + (1:Nout);

    R_k = R{k};
    R{k} = R_k; 
    t_k = qvec(k,:)';
    y = b_inner.x-b_inner.nx*sep;
    y = t_k+R_k*y;
    n_inner_k = R_k * b_inner.nx;
    rin(ii_in,:) = y';
    n_in(ii_in,:) = n_inner_k';
    % 
    %% Create collocation points

    x = t_k + R_k * b_outer.x;    % rot then transl, b just for vis
    n_k = R_k * b_outer.nx;
    E{k} = E0; 
    rout(ii_out,:) = x';
    n_out(ii_out,:) = n_k';
end


end

function self_test()
rng(1);

E0 = [0.7 0.45 1.1];
P = 3;
delta = 0.125;
N1 = 28;
N2 = 18;
sep = 0.06;

[rin, rout, qvec, ~, ~, ~, n_out] = getEllipsoidGrids(E0,P,delta,N1,N2,sep);

Nin = size(rin,1) / P;
Nout = size(rout,1) / P;
cols = lines(P);

figure('Name','getEllipsoidGrids self-test');
hold on
h_out = [];
h_in = [];
for k = 1:P
    ii_in = (k-1)*Nin + (1:Nin);
    ii_out = (k-1)*Nout + (1:Nout);
    h1 = scatter3(rout(ii_out,1),rout(ii_out,2),rout(ii_out,3),18,cols(k,:),'filled');
    h2 = scatter3(rin(ii_in,1),rin(ii_in,2),rin(ii_in,3),8,cols(k,:));
    if k == 1
        h_out = h1;
        h_in = h2;
    end
    plot3(qvec(k,1),qvec(k,2),qvec(k,3),'ko','MarkerFaceColor','k','MarkerSize',5);
end
qstep = max(floor(size(rout,1)/40),1);
quiver3(rout(1:qstep:end,1),rout(1:qstep:end,2),rout(1:qstep:end,3), ...
    n_out(1:qstep:end,1),n_out(1:qstep:end,2),n_out(1:qstep:end,3),0.2,'k');
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
title('getEllipsoidGrids self-test: 3-ellipsoid cluster')
legend([h_out, h_in],{'collocation','source'},'Location','bestoutside')

end

function qvec = cell_centers_to_matrix(qcell)
qvec = reshape(cell2mat(qcell(:).'),3,[]).';
end
