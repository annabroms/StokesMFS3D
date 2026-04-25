function res = rotate_vector_pages(x, R, npts)
%ROTATE_VECTOR_PAGES Apply one 3x3 rotation to each stacked block of vectors.
%
%   res = ROTATE_VECTOR_PAGES(x, R, npts)
%
%   INPUTS:
%       x    - stacked vector data of size (3*npts*P) x 1 or (3*npts) x P
%       R    - 3 x 3 x P array or 1 x P cell array of rotation matrices
%       npts - number of 3-vectors in each particle block
%
%   OUTPUT:
%       res  - rotated data with the same shape as x

if iscell(R)
    R = cat(3, R{:});
end

if isvector(x)
    x_pages = reshape(x, 3, npts, []);
    rotated = pagemtimes(R, x_pages);
    res = rotated(:);
    return;
end

if size(x, 1) ~= 3*npts
    error('rotate_vector_pages:bad_size', ...
        'Expected size(x,1) == 3*npts for matrix input.');
end

x_pages = reshape(x, 3, npts, size(x, 2));
rotated = pagemtimes(R, x_pages);
res = reshape(rotated, 3*npts, []);

end
