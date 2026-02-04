function T = stokes_DLP_mat(source, target, dir)
%STOKES_DLP_mat 3D Stokes double-layer potential (stresslet) matrix.
%
%   T = STOKES_DLP_mat(source, target, dir) returns the target-from-source
%   matrix for the free-space 3D Stokes double-layer potential (stresslet
%   kernel). The matrix maps force dipole strengths located at SOURCE points,
%   with dipole directions given by DIR, to velocities evaluated at TARGET
%   points. The viscosity is taken as mu = 1.
%
%   INPUTS:
%     source - Ns-by-3 real array of source locations s_k = (s_xk, s_yk, s_zk)
%     target - Nt-by-3 real array of target locations t_m = (t_xm, t_ym, t_zm)
%     dir    - Ns-by-3 real array of dipole directions d_k at source points
%
%   OUTPUT:
%     T      - (3*Nt)-by-(3*Ns) real matrix mapping dipole strengths
%              f = [f_x1; f_y1; f_z1; f_x2; f_y2; f_z2; ...; f_xNs; f_yNs; f_zNs]
%              to velocities
%              u = [u_x1; u_y1; u_z1; u_x2; u_y2; u_z2; ...; u_xNt; u_yNt; u_zNt]
%              via
%
%                  u = T * f .
%
%              Each 3-by-3 block corresponds to one target/source pair (t_m, s_k).
%              Block ordering is x-, y-, then z-components.
%
%   KERNEL (stresslet / double-layer kernel):
%     Let r = t - s be the separation vector between a target t and a source s,
%     with r = |r|. For a prescribed dipole direction d, the (velocity) stresslet
%     kernel used here is
%
%       T_ij(r, d) = -6/(8*pi*mu) * ( (r·d) r_i r_j ) / r^5 .
%
%     This routine assembles the dense matrix such that
%
%       u_i(t_m) = sum_k T_ij(t_m - s_k, d_k) f_j(s_k).
%
%     (Here mu = 1 is assumed.)
%
%   SINGULARITY / SELF TERMS:
%     The kernel is singular like 1/|r|^2 as r -> 0. No special treatment of
%     diagonal/self interactions is performed; the routine assumes that
%     source and target sets are disjoint.
%
%   NOTES:
%     This matrix is expensive to form and store for large Nt or Ns; it is
%     mainly useful for verification or small problems. Fast evaluation can
%     be achieved using FMM-based implementations.
%
%
%   Anna Broms, Feb 3, 2026

r    = source.';
rout = target.';

xsrc = r(1,:);
ysrc = r(2,:);
zsrc = r(3,:);
xtar = rout(1,:);
ytar = rout(2,:);
ztar = rout(3,:);

dx = dir(:,1);
dy = dir(:,2);
dz = dir(:,3);

mu = 1;

Nsrc = numel(xsrc);
Ntar = numel(xtar);
T = zeros(3*Ntar, 3*Nsrc);

prefac = -6 / (8*pi*mu);

for j = 1:Nsrc
    % Differences r = x - y_j (targets minus jth source)
    rx = xtar - xsrc(j);
    ry = ytar - ysrc(j);
    rz = ztar - zsrc(j);

    r2 = rx.^2 + ry.^2 + rz.^2;
    r5 = sqrt(r2).^5;

    % Dot with dipole direction at source j
    rdotd = rx*dx(j) + ry*dy(j) + rz*dz(j);

    % 3x3 block entries (each is Nt-by-1)
    Txx = prefac * rdotd .* (rx.*rx) ./ r5;
    Txy = prefac * rdotd .* (rx.*ry) ./ r5;
    Txz = prefac * rdotd .* (rx.*rz) ./ r5;

    Tyx = prefac * rdotd .* (ry.*rx) ./ r5;
    Tyy = prefac * rdotd .* (ry.*ry) ./ r5;
    Tyz = prefac * rdotd .* (ry.*rz) ./ r5;

    Tzx = prefac * rdotd .* (rz.*rx) ./ r5;
    Tzy = prefac * rdotd .* (rz.*ry) ./ r5;
    Tzz = prefac * rdotd .* (rz.*rz) ./ r5;

    % Assemble into block structure
    T(:, 3*j-2:3*j) = reshape( ...
        [Txx; Txy; Txz; ...
         Tyx; Tyy; Tyz; ...
         Tzx; Tzy; Tzz], ...
        3, 3*Ntar).';
end
end
