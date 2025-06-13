function [rin, rout, q, R, E] = getEllipsoidGrids(E0, P, delta, N1, N2, sep, R, q)
%GETELLIPSOIDGRIDS Generate MFS source and collocation grids on multiple ellipsoids.
%
%   [rin, rout, t, R, E] = GETELLIPSOIDGRIDS(E0, P, delta, N1, N2, sep, R, q)
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
%       q     - (Optional) 1 x P cell array of 3D translation vectors for each ellipsoid.
%               If omitted or empty, positions and orientations are
%               generated at random. 
%
%   OUTPUTS:
%       rin   - 3*P*N1 x 1 vector: concatenated source points for all particles.
%       rout  - 3*P*N2 x 1 vector: concatenated collocation points on all surfaces.
%       t     - P x 3 matrix of translation vectors (used or generated).
%       R     - 1 x P cell array of 3x3 rotation matrices (used or generated).
%       E     - 1 x P cell array of aspect ratios for the P particles.
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

% E0 = [.5 .5 1]; %set size /aspect ratio
% P=1;   % build cluster of P same shape ellipsoids
% delta = 0.01;    % target dist of each to prev ellipsoids
% rng(0);    % seed
% tic; 

%Create source points
if nargin<7 || nargin<8
    [E R q xnear] = ellipsoid_cluster(E0,P,delta);
end
%N1 = 40; 
%N2 = 30; 
rin = [];
rout = []; 

%Create proxy grid
b_inner = ellipsoid_param(E0(1),E0(2),E0(3));   % baseline object at the origin, aligned
b_inner = setupsurfquad(b_inner,[ceil(N1),ceil(N2)]);

%Create grid of collocation points
b_outer = ellipsoid_param(E0(1),E0(2),E0(3));   % baseline object at the oridin, aligned
b_outer = setupsurfquad(b_outer,[ceil(N1*1.15),ceil(N2*1.15)]);

for k = 1:P
    
    R_k = R{k};
    R{k} = R_k; 
    t_k = q{k}; 
    x = t_k + R_k * b_inner.x;    % rot then transl, b just for vis
    nx = R_k*b_inner.nx;
    
    y = b_inner.x-b_inner.nx*sep;
    y = t_k+R_k*y;
    
    rin = [rin; y']; 

%     figure()
% %   plot source points
%     plot3(y(1,:),y(2,:),y(3,:),'r.');
%     hold on
%     plot3(x(1,:),x(2,:),x(3,:),'b.');
%     quiver3(y(1,:),y(2,:),y(3,:),nx(1,:),nx(2,:),nx(3,:))
%     axis equal
    % 
    %% Create collocation points

    x = t_k + R_k * b_outer.x;    % rot then transl, b just for vis
    %Plot collocation points
    % plot3(b.x(1,:),b.x(2,:),b.x(3,:),'b.');
    % axis equal vis3d   % gonna need always
    E{k} = E0; 
    rout = [rout; x'];
end


end