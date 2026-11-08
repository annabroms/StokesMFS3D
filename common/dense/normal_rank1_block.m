function B = normal_rank1_block(nin, scale)
%NORMAL_RANK1_BLOCK Build blockdiag(scale*n_i*n_i') from node normals.
%
%   B = NORMAL_RANK1_BLOCK(nin, scale) returns a 3N-by-3N matrix with
%   3-by-3 diagonal block scale*nin(i,:)'*nin(i,:) for each node i.

if nargin == 0
    self_test;
    return;
end
if nargin < 2 || isempty(scale)
    scale = 1;
end

nin = validate_normals(nin);
N = size(nin,1);
B = zeros(3*N);

for i = 1:N
    rows = 3*(i-1) + (1:3);
    n = nin(i,:).';
    B(rows,rows) = scale * (n * n.');
end

end

function nin = validate_normals(nin)
nin = full(nin);
if size(nin,2) ~= 3
    error('normal_rank1_block:badNormals', ...
        'nin must be an N-by-3 array.');
end
nrms = vecnorm(nin,2,2);
if any(nrms <= eps)
    error('normal_rank1_block:badNormals', ...
        'nin contains a zero normal.');
end
nin = nin ./ nrms;
end

function self_test()
nin = [1 0 0; 0 1 0; 1 1 1];
N = size(nin,1);
B = normal_rank1_block(nin,1);

sym_err = norm(B - B','fro');
idem_err = norm(B*B - B,'fro');
rank_err = rank(B) - N;

fprintf('normal_rank1_block self-test:\n');
fprintf('  symmetry error = %.3e\n',sym_err);
fprintf('  idempotence error = %.3e\n',idem_err);
fprintf('  rank error = %d\n',rank_err);

if sym_err > 1e-14 || idem_err > 1e-14 || rank_err ~= 0
    error('normal_rank1_block:self_test', ...
        'Rank-one normal block failed its self-test.');
end
end
