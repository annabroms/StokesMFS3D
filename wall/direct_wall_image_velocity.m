function [velocity,roundoffScale] = direct_wall_image_velocity( ...
    basis,coefficients,targets,kernelOptions)
%DIRECT_WALL_IMAGE_VELOCITY Independent sum of complete wall source fields.
%   coefficients{b} is an N_b-by-3 physical block, already multiplied by
%   basis(b).scale by image_source_coefficients. Do not scale it again.
%   velocity is target-by-3. Optional roundoffScale is the per-target sum
%   of absolute constituents reported by the optional kernels, plus a
%   conservative source/wall-distance scale for the Stokeslet blocks.
%   The returned scale is diagnostic, not an error bound on the solution.
%   No-argument invocation checks wall cancellation and free-space signs.

if nargin == 0
    self_test(); velocity = []; roundoffScale = []; return
end
if numel(basis) ~= numel(coefficients)
    error('wall:ImageCoefficients','One coefficient array is needed per basis block.');
end
velocity = zeros(size(targets));
roundoffScale = zeros(size(targets,1),1);
for b = 1:numel(basis)
    block = basis(b);
    strength = coefficients{b};
    if ~isequal(size(strength),size(block.points))
        error('wall:ImageCoefficients','Coefficient dimensions differ in block %d.',b);
    end
    if isempty(strength), continue; end
    options = kernelOptions;
    options.diagnostics = nargout > 1;
    switch block.family
        case "S"
            u = blake.velocity(block.points.',strength.',targets.',options);
            if nargout > 1
                % Heights bound inverse reflected distances, including PN terms.
                height = block.points(:,3)-options.wall_height;
                targetHeight = targets(:,3)-options.wall_height;
                for j = 1:size(strength,1)
                    distance = vecnorm(targets-block.points(j,:),2,2);
                    scale = 1./distance+1/height(j)+ ...
                        targetHeight/height(j)^2;
                    roundoffScale = roundoffScale+ ...
                        16*norm(strength(j,:))*scale/(8*pi*options.mu);
                end
            end
        case "R"
            [u,info] = blake.rotlet(block.points.',strength.',targets.',options);
        case "T"
            options.normals = block.normalDirections.';
            [u,info] = blake.stresslet(block.points.',-strength.',targets.',options);
        case "D"
            [u,info] = blake.potential_dipole(block.points.',-2*strength.',targets.',options);
        otherwise
            error('wall:ImageFamily','Unknown image family %s.',block.family);
    end
    velocity = velocity+u.';
    if nargout > 1 && block.family ~= "S"
        roundoffScale = roundoffScale+info.constituent_scale(:);
    end
end
end

function self_test()
% Expose sign/prefactor errors through existing independent free-space kernels.
options = struct('method','blake','mu',1.3,'wall_height',0.35);
source = [0.2,-0.1,1.1];
target = [0.4,0.3,1.6;-0.2,0.5,0.7];
wall = [0.2,-0.1,0.35;0.7,0.3,0.35];
families = ["S","R","T","T","T","D"];
directions = [0,-1,0;-1,0,0;0,0,-1];
maxWall = 0; maxFree = 0;
for b = 1:6
    block = struct('family',families(b),'points',source, ...
        'normalDirections',zeros(0,3),'scale',3);
    if families(b) == "T", block.normalDirections = directions(b-2,:); end
    coeff = {[0.4,-0.7,0.2]};
    [u,scale] = direct_wall_image_velocity(block,coeff,wall,options);
    maxWall = max(maxWall,max(vecnorm(u,2,2)));
    if max(vecnorm(u,2,2)./max(scale,realmin)) > 1e-12
        error('wall:ImageWallTest','Wall cancellation failed for block %d.',b);
    end
    % A common large height leaves the free-space field for fixed separation.
    far = block; far.points = source+[0,0,1e5];
    value = direct_wall_image_velocity(far,coeff,target+[0,0,1e5],options);
    if families(b) == "S"
        reference = stokes_SLP_mat(source,target)/options.mu;
    else
        n = [1,0,0];
        if families(b) == "T", n = block.normalDirections; end
        [R,T,D] = stokes_singularity_matrices(source,target,n, ...
            logical([families(b)=="R",families(b)=="T",families(b)=="D"]));
        reference = [R,T,D]/options.mu;
    end
    expected = reference*coeff{1}.';
    freeError = norm(reshape(value.',[],1)-expected)/norm(expected);
    maxFree = max(maxFree,freeError);
    if freeError > 1e-4
        error('wall:ImageFreeTest','Free-space convention failed for block %d.',b);
    end
    zero = direct_wall_image_velocity(block,{zeros(1,3)},target,options);
    if any(zero(:)), error('wall:ImageZeroTest','Zero strength gave nonzero flow.'); end
end
fprintf('direct_wall_image_velocity PASS: wall %.3e; far-field convention error %.3e.\n', ...
    maxWall,maxFree);
end
