function A = stokes_mfs_velocity_matrix(basis,targets)
%STOKES_MFS_VELOCITY_MATRIX Assemble a combination of existing Stokes kernels.
%
% Each basis block has family, points (N-by-3), normalDirections (N-by-3),
% and scale.  Rows and each block's columns are interleaved by node:
% [ux1;uy1;uz1;ux2;...], [cx1;cy1;cz1;cx2;...].  Blocks are concatenated.
% S is a Stokeslet, R a rotlet, T a stresslet with the specified source
% normal, and D a potential dipole.  All existing kernels use viscosity 1.

ncol = 3*sum(arrayfun(@(block) size(block.points,1),basis));
A = zeros(3*size(targets,1),ncol);
col0 = 0;

for k = 1:numel(basis)
    block = basis(k);
    ns = size(block.points,1);
    cols = col0+(1:3*ns);
    switch block.family
        case "S"
            K = stokes_SLP_mat(block.points,targets);
        case "R"
            [K,~,~] = stokes_singularity_matrices( ...
                block.points,targets,eye(3),logical([1,0,0]));
        case "T"
            % Assemble the per-source normal contraction directly, keeping
            % one 3N-column block rather than three Cartesian blocks.
            assert(isequal(size(block.normalDirections),[ns,3]));
            [~,K,~] = stokes_singularity_matrices( ...
                block.points,targets,eye(3),logical([0,1,0]), ...
                block.normalDirections);
        case "D"
            [~,~,K] = stokes_singularity_matrices( ...
                block.points,targets,eye(3),logical([0,0,1]));
        otherwise
            error('stokes_mfs_velocity_matrix:BadFamily', ...
                'Unknown fundamental-solution family "%s".',block.family);
    end
    A(:,cols) = block.scale*K;
    col0 = col0+3*ns;
end
end
