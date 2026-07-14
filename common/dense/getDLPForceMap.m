function [Fmap, Kin, Kout, Kouter, T, W] = getDLPForceMap(rin, rout, nout, wout, q, opt)
%GETDLPFORCEMAP Build inner and outer DLP force/torque maps.
%
%   [Fmap, Kin, Kout, Kouter, T, W] = GETDLPFORCEMAP(rin,rout,nout,wout,q,opt)
%   returns the selected source-to-force map Fmap. If opt.outer_force is
%   true, Fmap = Kouter = T*W*Kout; otherwise Fmap = Kin.

if nargin < 6 || isempty(opt)
    opt = struct();
end
if ~isfield(opt,'outer_force')
    opt.outer_force = false;
end
if nargin < 5 || isempty(q)
    q = [0 0 0];
end

Kin = getKmat(rin,q);
Kout = getKmat(rout,q);

if nargin < 3 || isempty(nout)
    if isfield(opt,'ellipsoid') && logical(opt.ellipsoid)
        error('getDLPForceMap:missingNormals', ...
            'nout is required for ellipsoids.');
    end
    nout = rout - q(1,:);
    nout = nout ./ vecnorm(nout,2,2);
end

M = size(rout,1);
if nargin < 4 || isempty(wout)
    W = eye(3*M);
else
    wout = full(wout(:));
    if numel(wout) ~= M
        error('getDLPForceMap:badWeights', ...
            'wout must contain one quadrature weight for each row of rout.');
    end
    W = diag(repelem(wout,3));
end

T = stokes_DLP_mat(rout,rin,nout);
Kouter = T * W * Kout;

if logical(opt.outer_force)
    Fmap = Kouter;
else
    Fmap = Kin;
end

end
