function [y,iter_errv, iters, exact_errv] = lanczosSqrt(M,z,tol,max_it,pre_cond,yexact)
%LANCZOSSQRT Apply Lanczos to approximate sqrt(M)*z with optional preconditioning.
%
%   [y, iter_errv, iters, exact_errv] = lanczosSqrt(M, z, tol, max_it, pre_cond, yexact)
%
%   This routine uses a Lanczos/Krylov subspace method to approximate
%   y ≈ sqrt(M) * z for a symmetric (semi-)definite operator M. It builds the
%   tridiagonal Lanczos matrix H and applies sqrt(H) to the first basis vector.
%   An optional preconditioner is applied to the iterate y at each step.
%
%   INPUTS:
%     M         - matrix or function handle for the linear operator.
%     z         - RHS vector.
%     tol       - stopping tolerance based on relative change in y.
%     max_it    - maximum number of Krylov iterations.
%     pre_cond  - (optional) matrix or function handle applied to y each step.
%     yexact    - (optional) reference solution for exact error tracking.
%
%   OUTPUTS:
%     y         - approximate sqrt(M)*z (optionally preconditioned).
%     iter_errv - relative change per iteration (||y-yold||/||yold||).
%     iters     - number of iterations performed.
%     exact_errv- relative error vs yexact (if provided).
%
%   NOTES:
%     - M can be provided as a matrix or a function handle.
%     - The algorithm assumes M is symmetric (semi-)definite.
%     - pre_cond can be a matrix or a function handle.
%
% Borrowed from FBIM

if nargin < 5
    pre_cond = [];
end

if isa(M,'function_handle')
   afun = M; 
else 
   afun = @(x) M*x;
end

iter_err = 1;
v = z/norm(z);   % first Lanczos basis vector
y = v;

%yexact = sqrtm((M+M')/2)*z;
%yexact = real(yexact);

j=1;
while (iter_err > tol) && (j<max_it)
    
    % Lanczos step: apply operator and orthogonalize
    w = afun(v(:,j));
    h(j,j) = w'*v(:,j);
    
    
    % Three-term Lanczos recurrence; for j=1 there is no v_{j-1} term
    if j>1
       w = w - h(j-1,j)*v(:,j-1); 
    end
    w = w - h(j,j)*v(:,j);
    
    h(j+1,j) = norm(w);
    h(j,j+1) = h(j+1,j);
    v(:,j+1) = w/norm(w);
    
    % Apply sqrt of small tridiagonal matrix to the first basis vector:
    % In the Krylov basis, z/||z|| maps to e1, so we use sqrt(H)*e1 = hsqrt(:,1).
    % Then, project back to the original basis: y = V * hsqrt(:,1) * ||z||.
    hsqrt = sqrtm(h(1:j,1:j));
    yold = y;
    y = norm(z) * (v(:,1:j) * hsqrt(:,1));
    y = real(y);
    % if ~isempty(pre_cond)
    %     if isa(pre_cond,'function_handle')
    %         y = pre_cond(y);
    %     else
    %         y = pre_cond*y;
    %     end
    % end
    
   
    iter_err = norm(y-yold)/norm(yold);
    iter_errv(j) = iter_err;
    
    if nargin == 6
        exact_errv(j) = norm(y-yexact)/norm(yexact);
    end
    
    j = j+1;
%    [j,iter_err]    
end

if ~isempty(pre_cond)
    if isa(pre_cond,'function_handle')
        y = pre_cond(y);
    else
        y = pre_cond*y;
    end
end

if nargin<6
    exact_errv = [];
end

iters = j; 

end
