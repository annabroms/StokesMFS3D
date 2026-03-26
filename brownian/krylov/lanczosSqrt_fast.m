function [y,iter_errv,iters, exact_errv] = lanczosSqrt_fast(M,z,tol,max_it,pre_cond,yexact)
%LANCZOSSQRT_FAST Lanczos sqrt using tridiagonal eig for speed.
%
%   [y,iter_errv,iters,exact_errv] = lanczosSqrt_fast(M,z,tol,max_it,pre_cond,yexact)
%
%   Same interface as lanczosSqrt, but exploits the tridiagonal structure
%   of the Lanczos matrix to avoid full sqrtm.
%
%   See also: lanczosSqrt

if nargin < 5
    pre_cond = [];
end

if nargin < 1
    self_test;
    return;
end

if isa(M,'function_handle')
   afun = M;
else
   afun = @(x) M*x;
end

iter_err = 1;
v = z/norm(z);
y = v;

j = 1;
% Store tridiagonal entries only
alpha = zeros(max_it,1); % diagonal
beta = zeros(max_it,1);  % sub/super diagonal

while (iter_err > tol) && (j<max_it)
    % Lanczos step
    w = afun(v(:,j));
    alpha(j) = w' * v(:,j);

    if j>1
       w = w - beta(j-1) * v(:,j-1);
    end
    w = w - alpha(j) * v(:,j);

    beta(j) = norm(w);
    v(:,j+1) = w / beta(j);

    % Build tridiagonal Hj and apply sqrt(Hj) to e1 via eig
    d = alpha(1:j);
    e = beta(1:j-1);
    Hj = diag(d) + diag(e,1) + diag(e,-1);
    [Q, Lambda] = eig(Hj);
    s = Q' * [1; zeros(j-1,1)];
    hsqrt_e1 = Q * (sqrt(Lambda) * s);

    yold = y;
    y = norm(z) * (v(:,1:j) * hsqrt_e1);
    %y = real(y);

    if ~isempty(pre_cond)
        if isa(pre_cond,'function_handle')
            y = pre_cond(y);
        else
            y = pre_cond*y;
        end
    end

    iter_err = norm(y-yold)/norm(yold);
    iter_errv(j) = iter_err;

    if nargin == 6
        exact_errv(j) = norm(y-yexact)/norm(yexact);
    end

    j = j+1;
end

% if ~isempty(pre_cond)
%     if isa(pre_cond,'function_handle')
%         y = pre_cond(y);
%     else
%         y = pre_cond*y;
%     end
% end

if nargin<6
    exact_errv = [];
end

iters = j;

end

function self_test()
% Compare lanczosSqrt_fast against lanczosSqrt and sqrtm(M)*z.

rng(7);

n = 80;
[Q,~] = qr(randn(n));
d = linspace(0.5,3,n)';
M = Q*diag(d)*Q'; % SPD test matrix with controlled spectrum
Mfun = @(x) M*x;
z = randn(n,1);
y_dense = sqrtm(M)*z;

tol = 0;
max_it = n+1;

tic;
[y_ref, err_ref, iters_ref] = lanczosSqrt(Mfun, z, tol, max_it, [], y_dense);
t_ref = toc;

tic;
[y_fast, err_fast, iters_fast, exact_errv] = lanczosSqrt_fast(Mfun, z, tol, max_it, [], y_dense);
t_fast = toc;

rel_ref = norm(y_fast - y_ref) / norm(y_ref);
rel_dense = norm(y_fast - y_dense) / norm(y_dense);

fprintf('lanczosSqrt_fast self-test (n=%d):\n', n);
fprintf('  iters ref:  %d\n', numel(err_ref));
fprintf('  iters fast: %d\n', numel(err_fast));
fprintf('  iters ref flag:  %d\n', iters_ref);
fprintf('  iters flag: %d\n', iters_fast);
fprintf('  time ref:   %.3g s\n', t_ref);
fprintf('  time fast:  %.3g s\n', t_fast);
fprintf('  rel diff ref:   %.3e\n', rel_ref);
fprintf('  rel diff sqrtm: %.3e\n', rel_dense);
fprintf('  final exact err: %.3e\n', exact_errv(end));

if rel_ref >= 1e-12 || rel_dense >= 1e-12
    error('lanczosSqrt_fast:self_test', ...
        'Lanczos result did not match the SPD reference solutions on the self-test.');
end

fprintf('  self-test passed.\n');

end
