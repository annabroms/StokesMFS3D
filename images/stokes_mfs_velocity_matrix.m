function A = stokes_mfs_velocity_matrix(basis,targets)
%STOKES_MFS_VELOCITY_MATRIX Assemble a combination of existing Stokes kernels.
%
% Each basis block has family, points (N-by-3), normalDirections (N-by-3),
% and scale.  Rows and each block's columns are interleaved by node:
% [ux1;uy1;uz1;ux2;...], [cx1;cy1;cz1;cx2;...].  Blocks are concatenated.
% S is a Stokeslet, R a rotlet, T a stresslet with the specified source
% normal, and D a potential dipole.  All existing kernels use viscosity 1.
% T_n, T_n2, T_t label stresslet directions; Q_s uses the image-line axis.

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
        case {"T","T_n","T_n2","T_t"}
            % Assemble the per-source normal contraction directly, keeping
            % one 3N-column block rather than three Cartesian blocks.
            assert(isequal(size(block.normalDirections),[ns,3]));
            [~,K,~] = stokes_singularity_matrices( ...
                block.points,targets,eye(3),logical([0,1,0]), ...
                block.normalDirections);
        case "D"
            [~,~,K] = stokes_singularity_matrices( ...
                block.points,targets,eye(3),logical([0,0,1]));
        case "Q_s"
            K = stokes_axial_quadrupole_matrix( ...
                block.points,targets,block.normalDirections,"velocity");
        otherwise
            error('stokes_mfs_velocity_matrix:BadFamily', ...
                'Unknown fundamental-solution family "%s".',block.family);
    end
    A(:,cols) = block.scale*K;
    col0 = col0+3*ns;
end
end

function self_test()
% Verify three exact skew-strength modes and the six physical tensor modes.
b = struct('family',"T_t",'points',[0.1 0.2 0.3], ...
    'normalDirections',[1 0 0],'scale',1);
basis = repmat(b,3,1);
basis(2).family = "T_n";
basis(2).normalDirections = [0 1 0];
basis(3).family = "T_n2";
basis(3).normalDirections = [0 0 1];
x = [1.2 -0.2 0.8;0.4 1.1 -0.6;-0.8 0.3 1.1; ...
    0.7 -0.9 -0.5;1.1 0.9 0.4;-0.5 -1.3 0.8];
n = x./vecnorm(x,2,2);
A = stokes_mfs_velocity_matrix(basis,x);
grid.body.basis = basis;
for i = 1:2
    for j = i+1:3
        z = zeros(9,1);
        z(3*i-2:3*i) = basis(j).normalDirections.';
        z(3*j-2:3*j) = -basis(i).normalDirections.';
        if norm(A*z)>1e-14 || norm(image_pressure(grid,{z},x,1))>1e-14 || ...
                norm(image_traction(grid,{z},x,n,1),'fro')>1e-14
            error('stokes_mfs_velocity_matrix:SelfTest', ...
                'A skew-strength stresslet mode is not null.');
        end
    end
end
singularValues = svd(A);
if nnz(singularValues>1e-11*singularValues(1))~=6
    error('stokes_mfs_velocity_matrix:SelfTest', ...
        'The three-direction stresslet does not have rank six.');
end
legacy = basis;
legacy(1).family = "T"; legacy(2).family = "T"; legacy(3).family = "T";
if norm(A-stokes_mfs_velocity_matrix(legacy,x),'fro')>1e-14
    error('stokes_mfs_velocity_matrix:SelfTest','Directional labels changed the kernel.');
end
fprintf('stokes_mfs_velocity_matrix self-test passed (three stresslet null modes).\n');
end
