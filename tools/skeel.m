function c = skeel(A, p)
%SKEEL   Skeel's condition number.
%   SKEEL(A) estimates Skeel's condition number,
%
%      NORM(ABS(INV(A))*ABS(A), INF),
%
%   without computing ABS(INV(A))*ABS(A). SKEEL(A) is always less than or
%   equal to COND(A, INF). In practice, SKEEL(A) can be much less than
%   COND(A, INF). SKEEL(A) is invariant to row scaling.
%
%   SKEEL(A, P) directly computes Skeel's condition number in the P-norm:
%
%      NORM(ABS(INV(A))*ABS(A), P).
%
%   Example:
%
%      A = [1 0; 0 1e9];
%      cond(A)  % = 1e9
%      SKEEL(A) % = 1
%
%   References:
%
%      [1] Robert D. Skeel, "Scaling for numerical stability in Gaussian
%          elimination", J. ACM, 26 (1979), pp. 494-526.
%
%      [2] N. J. Higham, "FORTRAN codes for estimating the one-norm of a
%          real or complex matrix, with applications to condition
%          estimation", ACM Trans. Math. Soft., 14 (1988), pp. 381-396.
%
%   See also COND, CONDEST, NORM, NORMEST.

%   Copyright 2020 Dan Fortunato.

% Turn off warnings in case we encounter a singular matrix
warnState = warning('query', 'all');
temp = onCleanup(@() warning(warnState));
warning('off', 'all')

if ( nargin == 1 )
    % Estimate Skeel's condition number in the manner of LAPACK's DLA_GERCOND
    dA = decomposition(A, 'lu', 'CheckCondition', false);

    % Compute the equilibration matrix R such that inv(R)*A has unit 1-norm
    r = full(sum(abs(A), 2));

    % Estimate the norm of inv(A)
    n = size(A,1);
    v = zeros(n,1);
    x = zeros(n,1);
    isgn = zeros(n,1);
    c = 0;
    kase = 0;
    isave = zeros(3,1);
    [v, x, isgn, c, kase, isave] = dlacn2(n, v, x, isgn, c, kase, isave);
    while ( kase ~= 0 )
        if ( kase == 2 )
            x = dA\(x.*r);
        else
            x = (dA'\x).*r;
        end
        [v, x, isgn, c, kase, isave] = dlacn2(n, v, x, isgn, c, kase, isave);
    end
elseif ( isequal(p, 2) )
    if ( any(isinf(A)) )
        c = nan(class(A));
    else
        Ainv = inv(A);
        if ( any(isinf(Ainv)) )
            c = inf(class(A));
        else
            % Get the largest singular value
            c = svds(abs(Ainv)*abs(A), 1);
            if  ( c == 0 )
                c = inf(class(A));
            elseif ( isempty(c) )
                c = zeros(class(A));
            end
        end
    end
else
    % We'll let NORM pick up any invalid P argument
    c = norm(abs(inv(A))*abs(A), p);
end

end


function [v, x, isgn, est, kase, isave] = dlacn2(n, v, x, isgn, est, kase, isave)
%DLACN2   Estimate the 1-norm of a real square matrix.
%   This code is a MATLAB implementation of LAPACK's DLACN2.

itmax = 5;

if ( kase == 0 )
    x = ones(n,1)./n;
    kase = 1;
    isave(1) = 1;
    return
end

% Handle fall-through
if ( isave(1) < 1 || isave(1) > 5 )
    isave(1) = 1;
end

switch isave(1)
    case 1
        if ( n == 1 )
            v(1) = x(1);
            est = abs(v(1));
            kase = 0;
        else
            est = sum(abs(x));
            x = sign(x);
            isgn = round(x);
            kase = 2;
            isave(1) = 2;
        end

    case 2
        [~, isave(2)] = max(abs(x));
        isave(3) = 2;
        x = zeros(n,1);
        x(isave(2)) = 1;
        kase = 1;
        isave(1) = 3;

    case 3
        v = x;
        estold = est;
        est = sum(abs(v));
        if ( any(sign(x) ~= isgn) && est > estold )
            % Test for cycling
            x = sign(x);
            isgn = round(x);
            kase = 2;
            isave(1) = 4;
        else
            % Repeated sign vector detected, hence algorithm has converged
            x = (-1).^(0:n-1).*(1+(0:n-1)/(n-1)); x = x.';
            kase = 1;
            isave(1) = 5;
        end

    case 4
        jlast = isave(2);
        [~, isave(2)] = max(abs(x));
        kase = 1;
        if ( x(jlast) ~= abs(x(isave(2))) && isave(3) < itmax )
            isave(3) = isave(3)+1;
            x = zeros(n,1);
            x(isave(2)) = 1;
            isave(1) = 3;
        else
            x = (-1).^(0:n-1).*(1+(0:n-1)/(n-1)); x = x.';
            isave(1) = 5;
        end

    case 5
        temp = 2*sum(abs(x))/(3*n);
        if ( temp > est )
            v = x;
            est = temp;
        end
        kase = 0;

    otherwise
        error('Bad case.')
end

end
