function result = image_mobility_matvec(mu,grid,factors,batchSize)
%IMAGE_MOBILITY_MATVEC Apply the identity-diagonal mobility operator.
%   image_mobility_matvec() checks an explicit unequal two-body system.
%
% mu has one interleaved 3*M(k) block per body. Off-diagonal blocks are
% W_k*G_kj*(I-C_j*H_j)*A_j^+, with W_k from body.rowScale; diagonal
% blocks are identity. The complete
% global matrix and pair interaction matrices are never assembled here.

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
for j = 1:P
    rows = offsets(j)+(1:3*counts(j));
    z = factors(j).Y*(factors(j).U'*mu(rows));
    projected = z-factors(j).C*(factors(j).H*z);
    basis = grid.body(j).basis;
    coeff = image_source_coefficients(basis,projected);
    for k = 1:P
        if k == j
            continue
        end
        target = grid.body(k);
        for first = 1:batchSize:counts(k)
            nodes = first:min(first+batchSize-1,counts(k));
            targetRows = offsets(k)+reshape(3*nodes-2+(0:2).',[],1);
            velocity = direct_image_velocity( ...
                basis,coeff,target.collocation(nodes,:));
            velocity = target.rowScale(nodes).*velocity;
            result(targetRows) = result(targetRows)+reshape(velocity.',[],1);
        end
    end
end
end

function self_test()
%SELF_TEST compares direct interactions with explicit projected blocks.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','dense'));
base = [1,0,0;0,1,0;0,0,1;-1,0,0;0,-1,0;0,0,-1; ...
    1,1,1;-1,-1,-1];
centres = [0,0,0;3,0,0];
body = repmat(struct('collocation',[],'weights',[],'basis',[]),2,1);
factors = repmat(struct('H',[],'C',[],'Y',[],'U',[], ...
    'rank',0,'singularValues',[],'sourceCount',0, ...
    'collocationCount',0),2,1);
for k = 1:2
    body(k).collocation = centres(k,:)+base./vecnorm(base,2,2);
    body(k).weights = ones(8,1);
    if k == 2
        body(k).collocation(end+1,:) = centres(k,:)+[1,-1,0]/sqrt(2);
        body(k).weights(end+1) = 1;
    end
    body(k).basis = struct('family',"S", ...
        'points',centres(k,:)+0.4*base(1:4,:), ...
        'normalDirections',zeros(0,3),'scale',1);
    rotlet = body(k).basis;
    rotlet.family = "R";
    rotlet.points = centres(k,:)+[0.2,0,0];
    if k == 2
        rotlet.points = centres(k,:)+[0.2,0,0;0,0.2,0];
    end
    body(k).basis(2) = rotlet;
end
grid.body = body;
counts = arrayfun(@(b) 3*size(b.collocation,1),body);
counts = counts(:).';
offsets = [0,cumsum(counts)];
mu = sin((1:offsets(end)).');
for scaled = [true,false]
    for k = 1:2
        body(k).weights = (1:numel(body(k).weights)).';
        if scaled
            body(k).rowScale = sqrt(body(k).weights);
        else
            body(k).rowScale = ones(size(body(k).weights));
        end
        factors(k) = one_body_image_mobility_preconditioner( ...
            body(k),centres(k,:),1e-11); %#ok<AGROW>
    end
    grid.body = body;
    actual = image_mobility_matvec(mu,grid,factors,3);
    A = eye(offsets(end));
    for j = 1:2
        k = 3-j;
        sourceRows = offsets(j)+(1:counts(j));
        targetRows = offsets(k)+(1:counts(k));
        G = stokes_mfs_velocity_matrix(body(j).basis,body(k).collocation);
        A(targetRows,sourceRows) = repelem(body(k).rowScale,3).*G* ...
            (eye(factors(j).sourceCount)-factors(j).C*factors(j).H)* ...
            (factors(j).Y*factors(j).U');
    end
    expected = A*mu;
    assert(norm(actual-expected,inf) < 1e-11*max(1,norm(expected,inf)));
end
fprintf('image_mobility_matvec self-test passed.\n');
end
