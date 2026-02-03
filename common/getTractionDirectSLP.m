function t = getTractionDirectSLP(r, rout, nt, f, mu)
%getTractionDirectSLP  Matrix-free traction of 3D Stokes SLP (direct summation).
%
%   t = getTractionDirectSLP(r, rout, nt, f, mu)
%
% Computes traction t(x) = sigma(u,p)(x) * n(x) induced by Stokeslets of
% strengths f at sources r, evaluated at targets rout with TARGET normals nt.
%
% INPUTS:
%   r     - Nsrc×3 source points
%   rout  - Ntar×3 target points
%   nt    - Ntar×3 unit normals at targets (outward)
%   f     - (3*Nsrc)×1 vector of Stokeslet strengths [fx;fy;fz] stacked
%           or Nsrc×3 array (either is accepted)
%   mu    - viscosity (optional, default 1)
%
% OUTPUT:
%   t     - (3*Ntar)×1 traction vector stacked by target
%
% NOTES:
%   - Assumes no coincident source/target points.
%
% See also: getTractionFast
%
% Anna Broms, Jan 26, 2026

if nargin<1
    self_test;
    return;
elseif nargin < 5 
    mu = 1; 
end 

Nsrc = size(r,1);
Ntar = size(rout,1);

% reshape strengths to Nsrc×3
if isvector(f)
    f = reshape(f, 3, Nsrc).';   % Nsrc×3
elseif isequal(size(f), [Nsrc,3])
    % ok
else
    error('f must be (3*Nsrc)×1 or Nsrc×3.');
end
fx = f(:,1); fy = f(:,2); fz = f(:,3);

t = zeros(3*Ntar,1);

prefac = 6/(8*pi);   % traction kernel prefactor (no /mu)

for j = 1:Ntar
    % r = x - y, where x is target j, y are sources
    rx = rout(j,1) - r(:,1);
    ry = rout(j,2) - r(:,2);
    rz = rout(j,3) - r(:,3);

    r2 = rx.^2 + ry.^2 + rz.^2;
    %r5 = r2.^(5/2);
    r5 = sqrt(r2).^5;

    % dot with TARGET normal at x
    rdotn = rx*nt(j,1) + ry*nt(j,2) + rz*nt(j,3);  % Nsrc×1

    % s = r · f (Nsrc×1)
    rdotf = rx.*fx + ry.*fy + rz.*fz;

    % scalar weight for each source: prefac * (r·n)(r·f)/|r|^5
    w = prefac * (rdotn .* rdotf) ./ r5;   % Nsrc×1

    % traction vector at target j: sum_k w_k * r_k
    tj = [sum(w .* rx);
          sum(w .* ry);
          sum(w .* rz)];

    t(3*j-2:3*j) = tj;
end
end

function self_test()

ns = 2; 
nt = 1; 
rvec_in = rand(ns,3); 
rvec_out = rand(nt,3); 
nn = rand(nt,3);
mu = rand(3*ns,1); 

T = getTractionMat(rvec_in,rvec_out,nn); 
res1 = T*mu;

    res2 = getTractionDirectSLP(rvec_in, rvec_out, nn, mu);

fprintf('Traction error is %1.3e\n',norm(res1-res2,inf)); 

end
