function y = apply_normal_rank1(x, nin, scale)
%APPLY_NORMAL_RANK1 Apply blockdiag(scale*n_i*n_i') without forming it.
%
%   y = APPLY_NORMAL_RANK1(x, nin, scale) applies the per-node projector
%   (x_i dot n_i) n_i to each stacked 3-vector in x. Multiple right hand
%   sides are supported when x has size 3N-by-nrhs.

if nargin == 0
    self_test;
    return;
end
if nargin < 3 || isempty(scale)
    scale = 1;
end

nin = validate_normals(nin);
if mod(size(x,1),3) ~= 0
    error('apply_normal_rank1:badInput', ...
        'x must have 3N rows.');
end

N = size(x,1) / 3;
if size(nin,1) ~= N
    error('apply_normal_rank1:badNormals', ...
        'nin must have one row for each 3-vector in x.');
end

nrhs = size(x,2);
x3 = reshape(x,3,N,nrhs);
n3 = reshape(nin.',3,N,1);
dots = squeeze(sum(bsxfun(@times,x3,n3),1));
if nrhs == 1
    dots = dots(:);
end
y3 = bsxfun(@times,n3,reshape(dots,1,N,nrhs));
y = scale * reshape(y3,3*N,nrhs);

end

function nin = validate_normals(nin)
nin = full(nin);
if size(nin,2) ~= 3
    error('apply_normal_rank1:badNormals', ...
        'nin must be an N-by-3 array.');
end
nrms = vecnorm(nin,2,2);
if any(nrms <= eps)
    error('apply_normal_rank1:badNormals', ...
        'nin contains a zero normal.');
end
nin = nin ./ nrms;
end

function self_test()
nin = [1 0 0; 0 1 0; 1 1 1];
nin = nin ./ vecnorm(nin,2,2);
x = randn(9,4);
B = normal_rank1_block(nin,2.5);
y_dense = B * x;
y_fast = apply_normal_rank1(x,nin,2.5);
rel_err = norm(y_fast-y_dense,'fro') / max(norm(y_dense,'fro'),eps);

fprintf('apply_normal_rank1 self-test:\n');
fprintf('  dense action relative error = %.3e\n',rel_err);

if rel_err > 1e-14
    error('apply_normal_rank1:self_test', ...
        'Matrix-free normal rank-one action does not match dense action.');
end
end
