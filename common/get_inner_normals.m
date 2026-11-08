function nin = get_inner_normals(rin, q, opt, ~)
%GET_INNER_NORMALS Return stacked inner proxy normals for rank-one terms.
%
%   nin = GET_INNER_NORMALS(rin, q, opt, R) uses opt.nin or opt.n_in when
%   supplied. For spheres it can infer normals from rin-q. For ellipsoids
%   explicit normals are required because the proxy surface is not radial.

if nargin < 3 || isempty(opt)
    opt = struct();
end

nin = [];
if isfield(opt,'nin') && ~isempty(opt.nin)
    nin = opt.nin;
elseif isfield(opt,'n_in') && ~isempty(opt.n_in)
    nin = opt.n_in;
end

P = size(q,1);
N = size(rin,1) / P;
if abs(N-round(N)) > eps
    error('get_inner_normals:badSize', ...
        'rin must contain the same number of nodes for each particle.');
end
N = round(N);

if ~isempty(nin)
    nin = validate_and_normalize(nin,size(rin,1));
    return;
end

is_ellipsoid = isfield(opt,'ellipsoid') && logical(opt.ellipsoid);
if is_ellipsoid
    error('get_inner_normals:missingNormals', ...
        ['opt.add_rank1=true for ellipsoids requires explicit inner normals ', ...
         'in opt.nin or opt.n_in.']);
end

centers = kron(q,ones(N,1));
nin = validate_and_normalize(rin - centers,size(rin,1));

end

function nin = validate_and_normalize(nin,nrows)
nin = full(nin);
if size(nin,2) ~= 3 || size(nin,1) ~= nrows
    error('get_inner_normals:badNormals', ...
        'Inner normals must be a %d-by-3 array.',nrows);
end
nrms = vecnorm(nin,2,2);
if any(nrms <= eps)
    error('get_inner_normals:badNormals', ...
        'Inner normals contain a zero vector.');
end
nin = nin ./ nrms;
end
