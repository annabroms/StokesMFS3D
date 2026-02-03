function T = getTractionMat(r, rout, nn)
%GETTRACTIONMAT 3D traction matrix for the Stokes single-layer potential.
%
%   T = GETTRACTIONMAT(r, rout, nn) returns the target-from-source matrix for the
%   free-space 3D Stokeslet traction operator. The matrix maps point forces located
%   at source points R to surface tractions (stress dotted with the outward normal)
%   evaluated at target points ROUT with unit normals NN. The viscosity is taken as
%   mu = 1.
%
%   INPUTS:
%     r    - Ns-by-3 real array of source locations s_k = (s_xk, s_yk, s_zk)
%     rout - Nt-by-3 real array of target locations t_m = (t_xm, t_ym, t_zm)
%     nn   - Nt-by-3 real array of outward unit normals n_m at target points
%
%   OUTPUT:
%     T    - (3*Nt)-by-(3*Ns) real matrix mapping point forces
%            f = [f_x1; f_y1; f_z1; f_x2; f_y2; f_z2; ...; f_xNs; f_yNs; f_zNs]
%            to tractions
%            tau = [tau_x1; tau_y1; tau_z1; tau_x2; tau_y2; tau_z2; ...
%                   tau_xNt; tau_yNt; tau_zNt]
%            via
%
%                tau = T * f .
%
%            Block ordering is x-, y-, then z-components. Each 3-by-3 block
%            corresponds to one target/source pair (t_m, s_k).
%
%   NOTES:
%     * This routine assumes Stokes flow in 3D, with viscosity mu = 1.
%     * No special treatment of singular self-interactions (targets coinciding
%       with sources) is performed; the routine assumes source and target sets
%       are disjoint.
%     * The dense matrix is expensive to form and store for large Nt or Ns; it is
%       mainly useful for verification or small problems.
%
%   See also: GETTRACTIONDIRECTSLP, GETTRACTIONFAST
%
%   Anna Broms, Feb 3, 2026
%
%   Note: There is also FMM code to compare to (e.g. FMM3D, Flatiron).

 
r = r';
rout = rout';

xsrc = r(1,:);
ysrc = r(2,:);
zsrc = r(3,:);
xtar = rout(1,:);
ytar = rout(2,:);
ztar = rout(3,:);
nx = nn(:,1);
ny = nn(:,2);
nz = nn(:,3);
mu = 1;

Nsrc = numel(xsrc);
Ntar = numel(xtar);
T = zeros(3*Ntar,3*Nsrc);

prefac = 6 / (8*pi*mu);

for j = 1:Ntar
    %Vector differences r = x - y_j
    rx = xtar(j) - xsrc;
    ry = ytar(j) - ysrc;
    rz = ztar(j) - zsrc;
    r2 = rx.^2 + ry.^2 + rz.^2;
    r5 = sqrt(r2).^5;

    %Dot with normal at target
    rdotn = rx*nx(j) + ry*ny(j) + rz*nz(j);

    %Compute 3x3 block per target
    Txx = prefac * rdotn .* (rx.*rx) ./ r5;
    Txy = prefac * rdotn .* (rx.*ry) ./ r5;
    Txz = prefac * rdotn .* (rx.*rz) ./ r5;
    Tyx = prefac * rdotn .* (ry.*rx) ./ r5;
    Tyy = prefac * rdotn .* (ry.*ry) ./ r5;
    Tyz = prefac * rdotn .* (ry.*rz) ./ r5;
    Tzx = prefac * rdotn .* (rz.*rx) ./ r5;
    Tzy = prefac * rdotn .* (rz.*ry) ./ r5;
    Tzz = prefac * rdotn .* (rz.*rz) ./ r5;

    %Assemble into block structure
    T(3*j-2:3*j,:) = reshape( ...
    [Txx; Txy; Txz;  ...
     Tyx; Tyy; Tyz;  ...
     Tzx; Tzy; Tzz], ...
    3, 3*Nsrc);
end

end




 