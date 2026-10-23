function result = image_resistance_matvec(mu,grid,factors,batchSize)
%IMAGE_RESISTANCE_MATVEC Apply the row-scaled one-body-preconditioned operator.
%   result = image_resistance_matvec(mu,grid,factors,batchSize)
%   image_resistance_matvec() checks a tiny explicit two-body operator.
%
% mu/result have one interleaved [ux1;uy1;uz1;...] block per sphere, with
% block k of length 3*M(k). grid is from build_image_discretization and
% factors(k) is from one_body_image_preconditioner. The diagonal action is
% exactly identity; each off-diagonal contribution is W_k*G_kj*B_j^+*mu_j.
% W_k repeats the effective rowScale for x,y,z at each target. Each source-body
% interaction sums the selected image kernels directly in target batches, with
% no 3M-by-3N interaction matrix, global matrix, or FMM approximation.
% The source strengths and target velocities are interleaved by node.
% batchSize defaults to 128 target nodes.

if nargin == 0
    self_test();
    result = [];
    return
end
if nargin < 4 || isempty(batchSize)
    batchSize = 128;
end
P = numel(grid.body);
counts = arrayfun(@(body) size(body.collocation,1),grid.body);
counts = counts(:).';
offsets = [0,cumsum(3*counts)];
assert(numel(mu) == offsets(end));
result = mu(:);
for sourceBody = 1:P
    sourceRows = offsets(sourceBody)+(1:3*counts(sourceBody));
    sourceStrength = factors(sourceBody).Y* ...
        (factors(sourceBody).U'*mu(sourceRows));
    basis = grid.body(sourceBody).basis;
    coeff = image_source_coefficients(basis,sourceStrength);
    for targetBody = 1:P
        if targetBody == sourceBody
            continue
        end
        target = grid.body(targetBody);
        for first = 1:batchSize:counts(targetBody)
            nodes = first:min(first+batchSize-1,counts(targetBody));
            targetRows = offsets(targetBody)+ ...
                reshape(3*nodes-2+(0:2).',[],1);
            velocity = direct_image_velocity( ...
                basis,coeff,target.collocation(nodes,:));
            velocity = target.rowScale(nodes).*velocity;
            result(targetRows) = result(targetRows)+reshape(velocity.',[],1);
        end
    end
end
end

function self_test()
%SELF_TEST compares the matrix-free action with explicit cross blocks.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','dense'));
locations = [0,0,0;3,0,0];
base = [1,0,0;0,1,0;0,0,1;-1,0,0];
for k = 1:2
    body(k).collocation = base+locations(k,:); %#ok<AGROW>
    body(k).weights = (1:4).'; %#ok<AGROW>
    body(k).basis = struct('family',"S",'points',locations(k,:), ...
        'normalDirections',zeros(0,3),'scale',1); %#ok<AGROW>
end
mu = (1:24).'/24;
G12 = stokes_mfs_velocity_matrix(body(2).basis,body(1).collocation);
G21 = stokes_mfs_velocity_matrix(body(1).basis,body(2).collocation);
for scaled = [true,false]
    for k = 1:2
        if scaled
            body(k).rowScale = sqrt(body(k).weights);
        else
            body(k).rowScale = ones(4,1);
        end
        factors(k) = one_body_image_preconditioner( ...
            body(k),1e-12); %#ok<AGROW>
    end
    grid.body = body;
    actual = image_resistance_matvec(mu,grid,factors,2);
    expected = mu;
    expected(1:12) = expected(1:12)+ ...
        repelem(body(1).rowScale,3).*G12* ...
        (factors(2).Y*(factors(2).U'*mu(13:24)));
    expected(13:24) = expected(13:24)+ ...
        repelem(body(2).rowScale,3).*G21* ...
        (factors(1).Y*(factors(1).U'*mu(1:12)));
    assert(norm(actual-expected,inf) < 1e-13);
end
% Check every direct kernel, source scaling, and reuse of image-node geometry.
template = struct('family',"",'points',zeros(0,3), ...
    'normalDirections',zeros(0,3),'scale',1);
basis = repmat(template,4,1);
families = ["S","R","T","D"];
scale = [1,0.7,1.2,0.4];
proxy = [0,0,0;0.1,0.2,0.1];
image = [0.5,0,0;0.4,0.1,0];
coeff = cell(4,1);
strength = zeros(0,1);
for b = 1:4
    basis(b).family = families(b);
    basis(b).points = image;
    if b == 1
        basis(b).points = proxy;
    elseif b == 3
        basis(b).normalDirections = [1,0,0;0,1,0];
    end
    basis(b).scale = scale(b);
    raw = reshape(sin((1:6).'+b),3,[]).';
    coeff{b} = scale(b)*raw;
    strength = [strength;reshape(raw.',[],1)]; %#ok<AGROW>
end
targets = [1,0.2,0.1;1.1,-0.3,0.4;0.8,0.4,-0.5];
actual = reshape(direct_image_velocity(basis,coeff,targets).',[],1);
expected = stokes_mfs_velocity_matrix(basis,targets)*strength;
assert(norm(actual-expected,inf) < 1e-12*max(1,norm(expected,inf)));
fprintf('image_resistance_matvec self-test passed.\n');
end
