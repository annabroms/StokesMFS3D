function y = apply_block_diag(x, B, P, varargin)
%APPLY_BLOCK_DIAG Apply a repeated one-body block, optionally rotated per body.
%
%   y = APPLY_BLOCK_DIAG(x, B, P)
%   y = APPLY_BLOCK_DIAG(x, B, P, R)
%
%   x is (nb*P) x 1, B is nb x nb. If R is supplied as a 1 x P cell array
%   of 3 x 3 rotation matrices, then body k uses the rotated block
%   R_k * B * R_k' in the stacked vector basis.

R = [];
if ~isempty(varargin)
    R_in = varargin{1};
    if iscell(R_in)
        R = R_in;
    elseif isnumeric(R_in) && isequal(size(R_in), [3 3]) && P == 1
        R = {R_in};
    end
end

nb = length(x)/P;
y = zeros(size(x));

for k = 1:P
    idx = (k-1)*nb + (1:nb);
    if isempty(R)
        y(idx) = B * x(idx);
    else
        x_local = rotate_vector(x(idx), R{k}.');
        y(idx) = rotate_vector(B * x_local, R{k});
    end
end

end
