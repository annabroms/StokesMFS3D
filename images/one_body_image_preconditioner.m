function factor = one_body_image_preconditioner(body,cutoff)
%ONE_BODY_IMAGE_PRECONDITIONER Factor one complete enhanced self block.
%   factor = one_body_image_preconditioner(body,cutoff)
%   one_body_image_preconditioner() runs a small factor-action check.
%
% body is one entry from build_image_discretization: N-by-3 source blocks,
% M-by-3 collocation, physical area weights, and effective rowScale.
% cutoff is a relative singular-value threshold, default 5e-12. The block is
% B = diag(repelem(rowScale,3))*Gself, sized 3M-by-3N. factor.U and factor.Y
% implement B^+*mu = Y*(U'*mu), without forming B^+. factor.rank records
% the retained rank. Components are interleaved at each spatial node.

if nargin == 0
    self_test();
    factor = struct();
    return
end
if nargin < 2 || isempty(cutoff)
    cutoff = 5e-12;
end
validateattributes(cutoff,{'numeric'},{'scalar','positive','<',1});
assert(all(isfinite(body.weights) & body.weights > 0));
assert(numel(body.weights) == size(body.collocation,1));
assert(numel(body.rowScale) == size(body.collocation,1));
matrix = stokes_mfs_velocity_matrix(body.basis,body.collocation);
rowScale = repelem(body.rowScale,3);
matrix = rowScale.*matrix;
[Y,U,singularValues] = getPseudoFactors(matrix,cutoff,0);
assert(~isempty(singularValues),'The SVD cutoff removed every mode.');
factor.Y = Y;
factor.U = U;
factor.rank = numel(singularValues);
factor.singularValues = singularValues;
factor.N_total = size(matrix,2)/3;
factor.collocationCount = size(matrix,1)/3;
end

function self_test()
%SELF_TEST compares the factor application to a dense reference SVD.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','dense'));
body.collocation = [1,0,0;0,1,0;0,0,1;-1,0,0];
body.weights = [1;2;3;4];
body.basis = struct('family',"S",'points',[0,0,0], ...
    'normalDirections',zeros(0,3),'scale',1);
for scaled = [true,false]
    if scaled
        body.rowScale = sqrt(body.weights);
    else
        body.rowScale = ones(size(body.weights));
    end
    factor = one_body_image_preconditioner(body,1e-12);
    A = repelem(body.rowScale,3).* ...
        stokes_mfs_velocity_matrix(body.basis,body.collocation);
    rhs = (1:size(A,1)).';
    expected = pinv(A,1e-12*norm(A,2))*rhs;
    actual = factor.Y*(factor.U'*rhs);
    assert(norm(expected-actual) < 1e-12);
    assert(factor.rank == rank(A));
end
fprintf('one_body_image_preconditioner self-test passed.\n');
end
