function res = getStokesletFlow(tau_stokes, rin, rout, vars)
%GETSTOKESLETFLOW Evaluate Stokes flow velocity from Stokeslet sources.
%
%   res = GETSTOKESLETFLOW(tau_stokes, rin, rout, vars)
%
%   Computes the fluid velocity induced by a collection of Stokeslets located
%   at the source points `rin`, with strengths `tau_stokes`, evaluated at the
%   target points `rout`.
%
%   INPUTS:
%       tau_stokes : 3N x 1 vector of source strengths (Stokeslets), where N is
%                    the number of source locations. It is assumed to be stacked
%                    as [f1; f2; ...; fN] with each fi a 3D vector (x,y,z)
%       rin        : N x 3 matrix of source point locations.
%       rout       : M x 3 matrix of target point locations where the velocity is evaluated.
%       vars       : Struct containing optional settings:
%                    - vars.fmm  : Logical. If true, use FMM3D for fast evaluation.
%                    - vars.fmm_tol : Precision parameter used by FMM3D.
%
%   OUTPUT:
%       res        : 3M x 1 vector of fluid velocity at each target point, stacked
%                    as [u1; u2; ...; uM], where each ui is a 3D vector.
%
%   NOTES:
%       - If vars.fmm is true, the function uses the FMM3D interface `stfmm3d`
%         (requires FMM3D to be installed and compiled).
%       - Otherwise, it uses compiled direct summation via the FMM3D routine
%         `st3ddir_opt` (velocity only, OpenMP-parallel; the number of
%         threads is set by OMP_NUM_THREADS when MATLAB starts).
%       - Both branches use the standard 1/(8*pi) Stokeslet scaling.
%       - Modify this function for solving MFS with other fast summation
%         technique or implementation of direction summation.
%
%   See also: getTractionFast, getStressletFlow
%
%   DEPENDENCIES:
%       - FMM3D (https://github.com/flatironinstitute/fmm3d), built with
%         'make matlab': stfmm3d (accelerated) and st3ddir_opt (direct).
%
% Anna Broms June 12, 2025, updated Nov 21 2026

if nargin < 1
    self_test;
    return;
end

% Normalize common transposed inputs for mex compatibility
if size(rin,2) ~= 3 && size(rin,1) == 3
    rin = rin.';
end
if size(rout,2) ~= 3 && size(rout,1) == 3
    rout = rout.';
end
if ~isvector(tau_stokes) && size(tau_stokes,2) ~= 3 && size(tau_stokes,1) == 3
    tau_stokes = tau_stokes.';
end


% Source/target layout shared by both FMM3D routines
srcinfo.nd = 1;                                % One density (velocity only)
srcinfo.sources = rin';                        % Transpose: 3 x N
srcinfo.stoklet = reshape(tau_stokes, 3, []);  % Format: 3 x N
targ = rout';                                  % Transpose target: 3 x M
ifppregtarg = 1;                               % Velocity at targets only

if vars.fmm
    % -------- Use FMM3D for fast evaluation --------
    ifppreg = 0;                      % Don't compute potential at source
    fmm_tol = vars.fmm_tol;           % FMM precision (e.g., 1e-6)

    % Call FMM3D driver for Stokeslet potentials
    U = stfmm3d(fmm_tol, srcinfo, ifppreg, targ, ifppregtarg);
else
    % -------- Use compiled direct summation (FMM3D) --------
    U = st3ddir_opt(srcinfo, targ, ifppregtarg);
end

% Extract velocity result, stacked as [u1; u2; ...; uM]
res = U.pottarg(:);

clear U srcinfo;

function self_test()

if exist('stfmm3d','file') ~= 2 || exist('st3ddir_opt','file') ~= 2
    fprintf('getStokesletFlow self-test skipped: stfmm3d or st3ddir_opt not found.\n');
    return;
end

rng(1);
Nsrc = 40;
Ntar = 30;

rin = rand(Nsrc,3);
rout = rand(Ntar,3);
tau = rand(3*Nsrc,1);

vars.fmm = 0;
res_direct = getStokesletFlow(tau, rin, rout, vars);

vars.fmm = 1;
vars.fmm_tol = 1e-10;
res_fmm = getStokesletFlow(tau, rin, rout, vars);

rel_err = norm(res_fmm - res_direct) / norm(res_direct);
fprintf('getStokesletFlow self-test: rel err (fmm vs direct) = %.3e\n', rel_err);

end


end
