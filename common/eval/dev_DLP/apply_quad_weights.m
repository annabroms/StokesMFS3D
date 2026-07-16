function y = apply_quad_weights(x,w)
%APPLY_QUAD_WEIGHTS Apply scalar quadrature weights to stacked 3-vectors.
%
%   y = APPLY_QUAD_WEIGHTS(x,w) multiplies each consecutive three-vector
%   in x by the corresponding entry of w.

y = reshape(x,3,[]);
y = y .* reshape(w(:).',1,[]);
y = y(:);

end
