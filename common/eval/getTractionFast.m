function res = getTractionFast(tau_stokes, rin, rout, nn,vars)
%GETTRACTIONFAST Evaluate traction from Stokeslet sources (fast mex).
%
%   res = GETTRACTIONFAST(tau_stokes, rin, rout, vars)
%
%   Computes the traction due to a collection of Stokeslets located
%   at the source points `rin`, with strengths `tau_stokes`, evaluated at the
%   target points `rout`.
%
%   INPUTS:
%       tau_stokes : 3N x 1 vector of source strengths (Stokeslets), where N is
%                    the number of source locations. It is assumed to be stacked
%                    as [f1; f2; ...; fN] with each fi a 3D vector (x,y,z)
%       rin        : N x 3 matrix of source point locations.
%       rout       : M x 3 matrix of target point locations where the velocity is evaluated.
%       nn         : M x 3 matrix of normals at target points
%       vars       : Struct containing optional settings:
%                    - vars.fmm  : Logical. If true, use FMM3D for fast evaluation.
%                    - vars.fmm_tol : Precision parameter used by FMM3D.
%
%   OUTPUT:
%       res        : 3M x 1 vector of traction at each target point, stacked
%                    as [t1; t2; ...; tM], where each ti is a 3D vector.
%
%   NOTES:
%       - If vars.fmm is true, the function uses the FMM3D interface `stfmm3d`
%         (requires FMM3D to be installed and compiled). Currently, this
%         version is slow as it requires a lot more computation. 
%       - Otherwise, it uses compiled direct summation via the FMM3D routine
%         `st3ddir_opt` with target normals, which evaluates the Stokeslet
%         traction by its own kernel (OpenMP-parallel).
%       - Sign: the direct branch returns -sigma*n (the convention of the
%         former Stokes_Direct traction kernel), while the FMM branch
%         returns sigma*n, so the two branches differ in sign.
%       - Modify this function for solving MFS with other fast summation
%         technique or implementation of direction summation.
%
%   DEPENDENCIES:
%       - FMM3D, branch faster_st3ddir of https://github.com/annabroms/FMM3D,
%         built with 'make matlab': stfmm3d (evaluation via gradu and p)
%         and st3ddir_opt (direct traction, much faster; this branch only).
%
%   See also: getTractionDirectSLP
%
% Anna Broms Feb 3, 2025

if nargin < 1
    self_test;
    return;
end

if vars.fmm 
    % -------- Use FMM3D --------
    warning('too slow!')
    if ~isfield(vars,'fmm_tol')
        vars.fmm_tol = 1e-10;
    end
    ifppreg = 0;      % no eval at sources
    ifppregtarg = 3;  % eval at targets

    srcinfo.sources = rin';                 % 3 x Ns
    srcinfo.stoklet = reshape(tau_stokes,3,[]); % 3 x Ns

    targ = rout';                           % 3 x Nt
    targnor = nn';                          % 3 x Nt
    nt = size(rout,1);
    tic
    U = stfmm3d(vars.fmm_tol,srcinfo,ifppreg,targ,ifppregtarg);
    toc
    p = U.pretarg(:).';    % pressure (1 x Nt)
    gradu = U.gradtarg;    % grad vel (3 x 3 x Nt)
    shearstress = gradu + permute(gradu,[2 1 3]); % gradu + gradu^T
    % stress sigma = -pI + mu*(gradu+gradu^T), mu=1; traction T = sigma*n
    T = -(ones(3,1)*p).*targnor + ...
        squeeze(sum(shearstress .* reshape(kron(ones(3,1),targnor),3,3,nt),1));

    res = T(:);

else
    % -------- Use compiled direct summation (FMM3D) --------
    srcinfo.sources = rin';                       % 3 x Ns
    srcinfo.stoklet = reshape(tau_stokes, 3, []); % 3 x Ns

    % Traction only (flag 0) at targets rout with normals nn: sigma*n
    U = st3ddir_opt(srcinfo, rout', 0, nn');
    res = -U.tractarg(:); % keep the former Stokes_Direct sign, -sigma*n

end

function self_test()

if exist('stfmm3d','file') ~= 2 || exist('st3ddir_opt','file') ~= 2
    fprintf('getTractionFast self-test skipped: stfmm3d or st3ddir_opt not found.\n');
    return;
end

rng(1);
Nsrc = 4000;
Ntar = 3000;

rin = rand(Nsrc,3);
rout = rand(Ntar,3) + 0.3; % avoid near coincidences
nn = rand(Ntar,3);
nn = nn ./ vecnorm(nn,2,2);
tau = rand(3*Nsrc,1);

% Mexed direct traction
vars.fmm = 0;
tic;
res_mex = getTractionFast(tau, rin, rout, nn, vars);
t_mex = toc;

% FMM traction
vars.fmm = 1;
vars.fmm_tol = 1e-10;
tic;
res_fmm = getTractionFast(tau, rin, rout, nn, vars);
t_fmm = toc;

% Explicit traction matrix
tic;
T = getTractionMat(rin, rout, nn);
res_mat = T * tau;
t_mat = toc;

% Direct SLP (matrix-free) traction
tic;
res_slp = getTractionDirectSLP(rin, rout, nn, tau);
t_slp = toc;

rel_fmm_mex = norm(res_fmm + res_mex) / norm(res_mex);
rel_mat_mex = norm(res_mat - res_mex) / norm(res_mex);
rel_slp_mex = norm(res_slp - res_mex) / norm(res_mex);

fprintf('getTractionFast self-test timings (Ns=%d, Nt=%d):\n', Nsrc, Ntar);
fprintf('  mexed:  %.3g s\n', t_mex);
fprintf('  fmm:    %.3g s', t_fmm);
fprintf('  matrix: %.3g s\n', t_mat);
fprintf('  direct: %.3g s\n', t_slp);
fprintf('Relative errors vs mexed: fmm=%.3e, matrix=%.3e, direct=%.3e\n', ...
    rel_fmm_mex, rel_mat_mex, rel_slp_mex);

end


end
