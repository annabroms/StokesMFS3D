function [rbase_in, rbase_out] = getDesignGrid(Rp, opt)
%GETDESIGNGRID Constructs inner and outer spherical design grids for MFS
%spheres
%
%   [rbase_in, rbase_out] = GETDESIGNGRID(Rp, opt) generates two point sets 
%   on the sphere, one inner scaled by a given radius < 1 and one outer on
%   the unit sphere.
%
%   INPUTS:
%     Rp   - Scalar radius of the inner proxy surface.
%     opt  - Struct with the following fields:
%            opt.des_n   : Integer, controls the number of spherical design points on the proxy sphere 
%            opt.a_glob  : Scalar > 1, factor to determine the number of 
%                          outer grid points, M, relative to the inner proxy
%                          points, N, Typically, M \approx 1.2N. 
%
%   NOTE: Spherical designs do not exist for all integers, so the true
%         number of proxy points N will be close to des_n. The same is true
%         for choosing M.
%
%   OUTPUTS:
%     rbase_in  - (N x 3) matrix of 3D Cartesian coordinates for the inner spherical 
%                 design grid, scaled by Rp.
%     rbase_out - (M x 3) matrix of 3D Cartesian coordinates for the outer spherical 
%                 design grid on the unit sphere.
%                 
%
%   INTERNALS:
%     Uses get_sphdesign(n), assumed to return:
%        X - an n x 3 matrix of spherical design points on the unit sphere.
%        w - associated quadrature weights (not used here).
%
%   EXAMPLE USAGE:
%     opt.des_n = 700;
%     opt.a_glob = 1.2;
%     Rp = 0.68;
%     [r1, r2] = getDesignGrid(Rp, opt);
%
%   See also: GET_SPHDESIGN
%
% Anna Broms June 12, 2025

if nargin<1
    self_test();
    return;
end

[X, ~] = get_sphdesign(opt.des_n);
rbase_in = Rp * X;

[X, ~] = get_sphdesign(ceil(opt.a_glob * opt.des_n));
rbase_out = X;

end

function self_test()

% Parameters
Rp = 0.6;
opt.des_n = 100;
opt.a_glob = 1.2;

% Get grids
[rbase_in, rbase_out] = getDesignGrid(Rp, opt);

figure;
scatter3(rbase_in(:,1), rbase_in(:,2), rbase_in(:,3), 36, 'b', 'filled'); hold on;
scatter3(rbase_out(:,1), rbase_out(:,2), rbase_out(:,3), 36, 'r', 'filled');

axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
legend('Inner grid (radius 0.6)', 'Outer grid (radius 1)');
title('Spherical design grids: close to uniform distributions');
view(3);
grid on;

end



