function [y, iter_errv, iters, precond] = lanczosSqrt_precond(M, z, tol, max_it, k_pre, precond_in)
%LANCZOSSQRT_PRECOND Krylov sqrt with a low-rank recycling preconditioner.
%
%   [y, iter_errv, iters, precond] = lanczosSqrt_precond(M, z, tol, max_it, k_pre, precond_in)
%
%   Builds a low-rank preconditioner from the first k_pre Lanczos vectors,
%   then applies it inside lanczosSqrt. The preconditioner has the form
%   P(x) = x + U*(Q0'*x), similar in spirit to the recycling preconditioner
%   updates in precond_gmres.m (using U = Q*inv(H) - Q0).
%
%   INPUTS:
%     M       - matrix or function handle for the linear operator.
%     z       - RHS vector.
%     tol     - stopping tolerance for lanczosSqrt.
%     max_it  - maximum number of Krylov iterations.
%     k_pre     - number of Lanczos steps used to build the preconditioner.
%     precond_in- (optional) preconditioner struct to reuse (with Apply handle).
%
%   OUTPUTS:
%     y       - approximate sqrt(M)*z with preconditioning.
%     iter_errv - relative change per iteration from lanczosSqrt.
%     iters   - number of iterations performed.
%     precond - struct with Apply handle and stored U,Q0.
%
%   See also: lanczosSqrt, precond_gmres

if nargin < 1
    self_test;
    return;
end

if nargin < 5
    k_pre = 10;
end
if nargin < 6
    precond_in = [];
end

if isa(M,'function_handle')
    afun = M;
else
    afun = @(x) M*x;
end

% Build low-rank preconditioner from k_pre Lanczos steps (if not provided)
if isempty(precond_in)
    [Q, H] = lanczos_tridiag(afun, z, k_pre);
    Hk = H(1:k_pre, 1:k_pre);
    Q0 = Q(:,1:k_pre);
    U = Q(:,1:k_pre) * (Hk \ eye(k_pre)) - Q0;

    precond.U = U;
    precond.Q0 = Q0;
    precond.Apply = @(x) x + precond.U * (precond.Q0' * x);
else
    precond = precond_in;
end

% Run Krylov sqrt with the preconditioner
[y, iter_errv, iters] = lanczosSqrt(M, z, tol, max_it, precond.Apply);

end

function [Q, H] = lanczos_tridiag(afun, z, k)
%LANCZOS_TRIDIAG Build k-step Lanczos basis and tridiagonal H.

n = numel(z);
Q = zeros(n, k+1);
H = zeros(k+1, k);

q = z / norm(z);
Q(:,1) = q;

for j = 1:k
    w = afun(Q(:,j));
    H(j,j) = w' * Q(:,j);
    if j > 1
        w = w - H(j-1,j) * Q(:,j-1);
    end
    w = w - H(j,j) * Q(:,j);
    H(j+1,j) = norm(w);
    if H(j+1,j) == 0
        Q = Q(:,1:j);
        H = H(1:j,1:j);
        return;
    end
    Q(:,j+1) = w / H(j+1,j);
    H(j,j+1) = H(j+1,j);
end
end

function self_test()
% Compare baseline lanczosSqrt vs preconditioned variant.

rng(3);

n = 120;
A = randn(n);
M = A'*A + 0.5*eye(n); % SPD test matrix
z = randn(n,1);

tol = 1e-8;
max_it = 300;
k_pre = 5;

% Baseline (no preconditioner)
[y0, err0] = lanczosSqrt(M, z, tol, max_it);

% First run: build and return preconditioner
[y1, err1, ~, precond] = lanczosSqrt_precond(M, z, tol, max_it, k_pre);

% Second run: reuse the preconditioner inside lanczosSqrt_precond
[y2, err2] = lanczosSqrt_precond(M, z, tol, max_it, k_pre, precond);

it0 = numel(err0);
it1 = numel(err1);
it2 = numel(err2);

rel_y1 = norm(y1 - y2) / norm(y2);
rel_y0 = norm(y0 - y2) / norm(y2);

fprintf('lanczosSqrt_precond self-test (n=%d):\n', n);
fprintf('  iters no precond: %d\n', it0);
fprintf('  iters precond   : %d\n', it1);
fprintf('  iters reuse     : %d\n', it2);
fprintf('  rel diff (precond build vs reuse) = %.3e\n', rel_y1);
fprintf('  rel diff (new vs baseline)= %.3e\n', rel_y0);

% Quick convergence summary
fprintf('  final iter_err: baseline %.3e, precond %.3e, reuse %.3e\n', ...
    err0(end), err1(end), err2(end));

if rel_y1 < 1e-10
    fprintf('  y reproduced within tolerance (1e-10) vs precond run.\n');
else
    fprintf('  y NOT reproduced within tolerance (1e-10) vs precond run.\n');
end

end
