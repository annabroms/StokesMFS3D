function [y,iter_errv,iters,exact_errv] = arnoldiSqrt_fast(M,z,tol,max_it,pre_cond,yexact,orth_method)
%ARNOLDISQRT_FAST Arnoldi approximation of sqrt(M)*z.
%
%   [y,iter_errv,iters,exact_errv] = arnoldiSqrt_fast(M,z,tol,max_it,pre_cond,yexact,orth_method)
%
%   Same interface as lanczosSqrt_fast, but uses full Arnoldi instead
%   of Lanczos so that M need not be symmetric positive definite. The
%   projected square root is evaluated with sqrtm on the Hessenberg matrix.
%
%   orth_method is optional and may be 'householder' (default) or 'mgs'.
%   For backward compatibility, a string passed in place of yexact is
%   interpreted as orth_method.

if nargin < 5
    pre_cond = [];
end

if nargin < 1
    self_test;
    return;
end

if nargin < 6
    yexact = [];
end

if nargin < 7 || isempty(orth_method)
    orth_method = 'householder';
end

if is_text_scalar(yexact)
    orth_method = yexact;
    yexact = [];
end

orth_method = parse_orth_method(orth_method);
track_exact = ~isempty(yexact);

if isa(M,'function_handle')
    afun = M;
else
    afun = @(x) M*x;
end

znorm = norm(z);
if znorm == 0
    y = z;
    iter_errv = zeros(0,1);
    iters = 0;
    if track_exact
        exact_errv = zeros(0,1);
    else
        exact_errv = [];
    end

    if ~isempty(pre_cond)
        if isa(pre_cond,'function_handle')
            y = pre_cond(y);
        else
            y = pre_cond*y;
        end
    end
    return;
end

n = length(z);
max_it = min(max_it,n+1);
iter_err = 1;
V = zeros(n,max_it);
H = zeros(max_it,max_it-1);
V(:,1) = z/znorm;
y = V(:,1);
iter_errv = zeros(max_it-1,1);

if track_exact
    exact_errv = zeros(max_it-1,1);
else
    exact_errv = [];
end

if strcmp(orth_method,'householder')
    orth_state = init_householder_state(V(:,1),max_it);
else
    orth_state = [];
end

j = 1;
step_count = 0;
while (iter_err > tol) && (j < max_it)
    w = afun(V(:,j));

    if strcmp(orth_method,'householder')
        [H,V,orth_state,happy_breakdown] = householder_arnoldi_step(w,V,H,orth_state,j);
    else
        [H,V,happy_breakdown] = mgs_arnoldi_step(w,V,H,j);
    end

    Hj = H(1:j,1:j);
    hsqrt = sqrtm(Hj);

    yold = y;
    y = znorm * (V(:,1:j) * hsqrt(:,1));

    if ~isempty(pre_cond)
        if isa(pre_cond,'function_handle')
            y = pre_cond(y);
        else
            y = pre_cond*y;
        end
    end

    iter_err = norm(y-yold) / max(norm(yold),eps); %This is how convergence is measured in the Chow and Saad paper.
    iter_errv(j) = iter_err;

    if track_exact
        exact_errv(j) = norm(y-yexact) / max(norm(yexact),eps);
    end

    step_count = j;

    if happy_breakdown
        break;
    end

    j = j+1;
end

iter_errv = iter_errv(1:step_count);
if track_exact
    exact_errv = exact_errv(1:step_count);
end



iters = j;

end

function [H,V,happy_breakdown] = mgs_arnoldi_step(w,V,H,j)
for i = 1:j
    H(i,j) = V(:,i)' * w;
    w = w - H(i,j) * V(:,i);
end

H(j+1,j) = norm(w);
happy_breakdown = H(j+1,j) <= eps(max(1,norm(H(1:j+1,j),inf)));

if ~happy_breakdown
    V(:,j+1) = w / H(j+1,j);
end
end

function state = init_householder_state(q1,max_it)
n = length(q1);
state.u = zeros(n,max_it);
state.tau = zeros(max_it,1);
state.signs = ones(max_it,1);

[state.u(:,1), state.tau(1), alpha1] = make_householder_from_tail(q1,1,n);
state.signs(1) = sign_nonzero(alpha1);
end

function [H,V,state,happy_breakdown] = householder_arnoldi_step(w,V,H,state,j)
transformed = w;
for k = 1:j
    transformed = apply_householder(transformed,state.u(:,k),state.tau(k));
end

H(1:j,j) = state.signs(1:j) .* transformed(1:j);
tail = transformed(j+1:end);
tail_norm = norm(tail);
probe_len = min(j+1,length(transformed));
happy_breakdown = tail_norm <= eps(max(1,norm(transformed(1:probe_len),inf)));

if happy_breakdown
    H(j+1,j) = 0;
    return;
end

H(j+1,j) = tail_norm;
[state.u(:,j+1), state.tau(j+1), alpha] = make_householder_from_tail(tail,j+1,size(V,1));
state.signs(j+1) = sign_nonzero(alpha);
V(:,j+1) = state.signs(j+1) * build_householder_basis_vector(state,j+1,size(V,1));
end

function q = build_householder_basis_vector(state,j,n)
q = zeros(n,1);
q(j) = 1;
for k = j:-1:1
    q = apply_householder(q,state.u(:,k),state.tau(k));
end
end

function [u,tau,alpha] = make_householder_from_tail(tail,k,n)
[v,tau,alpha] = householder_vector(tail);
u = zeros(n,1);
u(k:n) = v;
end

function [v,tau,alpha] = householder_vector(x)
v = zeros(size(x));
tau = 0;
alpha = 0;

xnorm = norm(x);
if xnorm == 0
    return;
end

alpha = -sign_nonzero(x(1)) * xnorm;
v = x;
v(1) = v(1) - alpha;
tau = 2 / (v' * v);
end

function x = apply_householder(x,u,tau)
if tau ~= 0
    x = x - (tau * (u' * x)) * u;
end
end

function s = sign_nonzero(x)
if x < 0
    s = -1;
else
    s = 1;
end
end

function tf = is_text_scalar(x)
tf = ischar(x) || (isstring(x) && isscalar(x));
end

function orth_method = parse_orth_method(orth_method)
if isstring(orth_method)
    orth_method = char(orth_method);
end

orth_method = lower(strtrim(orth_method));
if ~ismember(orth_method,{'householder','mgs'})
    error('arnoldiSqrt_fast:orth_method', ...
        'orth_method must be ''householder'' or ''mgs''.');
end
end

function self_test()
% Compare both Arnoldi orthogonalizations against the SPD Lanczos fast path.

rng(7);

n = 160;
A = randn(n);
M = A'*A + 0.5*eye(n);
Mfun = @(x) M*x;
z = randn(n,1);

tol = 0;
max_it = 40;

[y_lanczos, err_lanczos, iters_lanczos] = lanczosSqrt_fast(Mfun,z,tol,max_it);
[y_house, err_house, iters_house] = arnoldiSqrt_fast(Mfun,z,tol,max_it);
[y_mgs, err_mgs, iters_mgs] = arnoldiSqrt_fast(Mfun,z,tol,max_it,[],'mgs');

rel_house = norm(y_house-y_lanczos) / norm(y_lanczos);
rel_mgs = norm(y_mgs-y_lanczos) / norm(y_lanczos);

fprintf('arnoldiSqrt_fast self-test (SPD, n=%d):\n', n);
fprintf('  Lanczos steps:     %d\n', numel(err_lanczos));
fprintf('  Householder steps: %d\n', numel(err_house));
fprintf('  MGS steps:         %d\n', numel(err_mgs));
fprintf('  iters Lanczos:     %d\n', iters_lanczos);
fprintf('  iters Householder: %d\n', iters_house);
fprintf('  iters MGS:         %d\n', iters_mgs);
fprintf('  rel diff house:    %.3e\n', rel_house);
fprintf('  rel diff mgs:      %.3e\n', rel_mgs);

if rel_house >= 1e-10 || rel_mgs >= 1e-10
    error('arnoldiSqrt_fast:self_test', ...
        'Arnoldi result did not match lanczosSqrt_fast on the SPD self-test.');
end

fprintf('  self-test passed.\n');

end
