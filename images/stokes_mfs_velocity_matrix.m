function A = stokes_mfs_velocity_matrix(basis,targets)
%STOKES_MFS_VELOCITY_MATRIX Assemble a combination of existing Stokes kernels.
%
% Each basis block has family, points (N-by-3), normalDirections (N-by-3),
% and scale.  Rows and each block's columns are interleaved by node:
% [ux1;uy1;uz1;ux2;...], [cx1;cy1;cz1;cx2;...].  Blocks are concatenated.
% S is a Stokeslet, R a rotlet, T a stresslet with the specified source
% normal, and D a potential dipole.  All existing kernels use viscosity 1.
% T_n and T_t label independent stresslet blocks with supplied normals.

if nargin == 0
    self_test(); A = []; return
end
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
        case {"T","T_n","T_t"}
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

function self_test()
% Verify both directional blocks and their exact shared stresslet mode.
b = struct('family',"T_n",'points',[0.1 0.2 0.3], ...
    'normalDirections',[0 0 1],'scale',1);
c = b; c.family = "T_t"; c.normalDirections = [1 0 0];
x = [1.2 -0.2 0.8;0.4 1.1 -0.6];
n = x./vecnorm(x,2,2);
basis = [b;c];
A = stokes_mfs_velocity_matrix(basis,x);
z = [c.normalDirections.';-b.normalDirections.'];
grid.body.basis = basis;
if norm(A*z)>1e-14 || norm(image_pressure(grid,{z},x,1))>1e-14 || ...
        norm(image_traction(grid,{z},x,n,1),'fro')>1e-14
    error('stokes_mfs_velocity_matrix:SelfTest','Shared stresslet mode is not null.');
end
basis(1).family = "T"; basis(2).family = "T";
if norm(A-stokes_mfs_velocity_matrix(basis,x),'fro')>1e-14
    error('stokes_mfs_velocity_matrix:SelfTest','Directional labels changed the kernel.');
end
fprintf('stokes_mfs_velocity_matrix self-test passed (joint stresslet null mode).\n');
end
