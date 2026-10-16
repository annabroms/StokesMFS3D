function A = lap3d_mfs_potential_matrix(targets,grid)
%LAP3D_MFS_POTENTIAL_MATRIX Shared assembly for fitting and checking Laplace MFS.
% Coefficient columns are [all source monopoles; enhancing normal dipoles].
% The monopole order is exactly grid.sources.  Dipole order follows
% grid.dipoleSourceIndices, with the corresponding grid.dipoleNormals.
% L = grid.dipoleLengthScale multiplies DLP columns: moment = L*d*n.
% Hence fitted q and d have common charge units.  There are no quadrature
% weights or fixed coupling between the independently fitted q and d.
% Builder/solver supplies empty dipole metadata for monopole-only grids.

A = lap3d_charge_potential_matrix(targets,grid.sources);
if isempty(grid.dipoleSourceIndices)
    return
end
dipoleSources = grid.sources(grid.dipoleSourceIndices,:);
D = lap3d_charge_potential_matrix(targets,dipoleSources,grid.dipoleNormals);
A = [A,grid.dipoleLengthScale*D];
end
