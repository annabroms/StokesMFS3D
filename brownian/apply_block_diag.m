function y = apply_block_diag(x, B, P, V, D)
%APPLY_BLOCK_DIAG Apply block-diagonal operator with P identical blocks.
%
%   y = APPLY_BLOCK_DIAG(x, B, P)
%   y = APPLY_BLOCK_DIAG(x, B, P, V, D)
%
%   x is (nb*P) x 1, B is nb x nb. If V and D are provided, applies
%   y_k = V * D * V' * x_k (equivalent to B*x_k when B = V*D*V').
%
%   D may be a vector of diagonal entries or a diagonal matrix.

if nargin < 4
    V = [];
    D = [];
end

nb = size(B,1);
y = zeros(size(x));

use_eig = ~isempty(V) && ~isempty(D);
if use_eig && isvector(D)
    D = diag(D);
end

for k = 1:P
    idx = (k-1)*nb + (1:nb);
    if use_eig
        y(idx) = V * (D * (V' * x(idx)));
    else
        y(idx) = B * x(idx);
    end
end
end
