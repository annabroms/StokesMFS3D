function factor = drop_image_one_body_preconditioner(body,mu,c,opt)
%DROP_IMAGE_ONE_BODY_PRECONDITIONER Factor the scaled local drop block.
if nargin == 0
    self_test(); factor = struct(); return
end
factor = drop_image_scaled_self_block(body,mu,c,opt);
A = factor.selfMatrix;
if opt.scale_columns
    scales = vecnorm(A,2,1).';
    zeroCount = nnz(scales==0);
    scales(scales==0) = 1;
else
    scales = ones(size(A,2),1); zeroCount = 0;
end
S = A./scales.';
[m,n] = size(S);
storage = whos('S');
r = min(m,n);
factorMiB = (m+n+r)*r*(storage.bytes/numel(S))/2^20;
fprintf(['drop one-body SVD: %d x %d, matrix %.2f MiB, ', ...
    'economy factors %.2f MiB estimated\n'], ...
    m,n,storage.bytes/2^20,factorMiB);
% Only the singular vectors paired with min(m,n) singular values enter the inverse.
[L,D,R] = svd(S,'econ');
s = diag(D);
keep = s>opt.svdCutoff*s(1);
factor.leftTranspose = L(:,keep).';
factor.singularValues = s(keep);
factor.rightVectors = R(:,keep);
factor.columnScales = scales;
factor.zeroColumnCount = zeroCount;
factor.rank = nnz(keep);
factor.svdTolerance = opt.svdCutoff;
fprintf('drop one-body retained rank %d/%d at relative cutoff %.3e\n', ...
    factor.rank,numel(s),opt.svdCutoff);
end

function self_test()
% Verify local factor action against the explicit truncated SVD projection.
v = [eye(3);-eye(3)]; x = [1 1 1;1 -1 -1;-1 1 -1;-1 -1 1]/sqrt(3);
b.family = "S"; b.points = 0.6*v; b.normalDirections = zeros(0,3); b.scale = 1;
body.basis = b; b.points = 2*v; body.interiorBasis = b;
body.collocation = x; body.normals = x; body.rowScale = ones(4,1);
opt = drop_image_options(struct('N',30),1,1);
f = drop_image_one_body_preconditioner(body,2,1/(8*pi),opt);
y = sin((1:24).');
z = applyOneBodyPinv(y,f);
projected = f.leftTranspose.'*(f.leftTranspose*y);
assert(norm(f.selfMatrix*z-projected)<1e-10);
fprintf('drop_image_one_body_preconditioner self-test passed.\n');
end
