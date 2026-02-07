function [y,iter_errv,iters, exact_errv] = KrylovSqrtMsing_fast(M,z,tol,max_it,pre_cond,yexact)
%KRYLOVSQRTMSING_FAST Lanczos sqrt using tridiagonal eig for speed.
%
%   [y,iter_errv,iters,exact_errv] = KrylovSqrtMsing_fast(M,z,tol,max_it,pre_cond,yexact)
%
%   Same interface as KrylovSqrtMsing, but exploits the tridiagonal structure
%   of the Lanczos matrix to avoid full sqrtm.
%
%   See also: KrylovSqrtMsing

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
    y = real(y);

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

if nargin<6
    exact_errv = [];
end

iters = j;

end

function self_test()
% Compare KrylovSqrtMsing_fast vs KrylovSqrtMsing

rng(7);

n = 5000;
A = randn(n);
M = A'*A + 0.5*eye(n); % SPD test matrix
z = randn(n,1);

tol = 1e-10;
max_it = 100;

tic;
[y_ref, err_ref] = KrylovSqrtMsing(M, z, tol, max_it);
t_ref = toc;

tic;
[y_fast, err_fast] = KrylovSqrtMsing_fast(M, z, tol, max_it);
t_fast = toc;

rel_y = norm(y_fast - y_ref) / norm(y_ref);

fprintf('KrylovSqrtMsing_fast self-test (n=%d):\n', n);
fprintf('  iters ref:  %d\n', numel(err_ref));
fprintf('  iters fast: %d\n', numel(err_fast));
fprintf('  time ref:   %.3g s\n', t_ref);
fprintf('  time fast:  %.3g s\n', t_fast);
fprintf('  rel diff y: %.3e\n', rel_y);

if rel_y < 1e-10
    fprintf('  y reproduced within tolerance (1e-10).\n');
else
    fprintf('  y NOT reproduced within tolerance (1e-10).\n');
end

end
