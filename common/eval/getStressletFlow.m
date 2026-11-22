function u = getStressletFlow(rvec_in,rvec_out,nn,strslet,Ns,vars)
%GETSTRESSLETFLOW Evaluate Stokes flow velocity from stresslet sources.
%
%   u = GETSTRESSLETFLOW(rvec_in, rvec_out, nn, strslet, Ns, vars)
%
%   Computes the velocity induced by a collection of stresslets located at
%   `rvec_out`, evaluated at the target points `rvec_in`. The stresslet
%   strengths are in `strslet`, and `nn` provides the stresslet orientation
%   vectors (one per source). The layout matches the FMM3D stresslet
%   interface, and the direct branch uses a mexed stresslet kernel.
%
%   INPUTS:
%       rvec_in  : Nt x 3 target points.
%       rvec_out : Ns x 3 source points.
%       nn       : Ns x 3 stresslet orientation vectors.
%       strslet  : 3*Ns x 1 or Ns x 3 stresslet strengths.
%       Ns        : Number of sources 
%       vars     : Struct with options:
%                  - vars.fmm : Logical. If true, use FMM3D for fast evaluation.
%                  - vars.fmm_tol : Precision parameter for FMM3D.
%
%   OUTPUT:
%       u        : 3*Nt x 1 stacked velocity at targets.
%
%   NOTES:
%       - If vars.fmm is true, uses stfmm3d stresslet evaluation.
%       - Otherwise, uses compiled direct summation via the FMM3D routine
%         `st3ddir_opt` (stresslets only, velocity only, OpenMP-parallel).
%       - FMM3D's stresslet has the opposite sign to the one used here, so
%         both branches negate its output.
%
%   See also: getStokesletFlow, getTractionFast

if nargin < 1
    self_test;
    return;
end

% Normalize common transposed inputs for mex compatibility
if size(rvec_out,2) ~= 3 && size(rvec_out,1) == 3
    rvec_out = rvec_out.';
end
if size(rvec_in,2) ~= 3 && size(rvec_in,1) == 3
    rvec_in = rvec_in.';
end
if size(nn,2) ~= 3 && size(nn,1) == 3
    nn = nn.';
end
if ~isvector(strslet) && size(strslet,2) ~= 3 && size(strslet,1) == 3
    strslet = strslet.';
end

if nargin<6
    vars.fmm = 1; 
end


% Source/target layout shared by both FMM3D routines
targ = rvec_out';                          % 3 x Nt
srcinfo.sources = rvec_in';                % 3 x Ns
srcinfo.strslet = reshape(strslet, 3, []); % Format: 3 x Ns
srcinfo.strsvec = nn';                     % 3 x Ns
ifppregtarg = 1;                           % just vel out

if vars.fmm
    ifppreg = 0;      % no eval at sources
    U = stfmm3d(vars.fmm_tol,srcinfo,ifppreg,targ,ifppregtarg);
else
    % -------- Use compiled direct summation (FMM3D) --------
    U = st3ddir_opt(srcinfo, targ, ifppregtarg);
end

u = -U.pottarg(:); % there is a sign difference with FMM3D's stresslet

end

function self_test()

if exist('stfmm3d','file') ~= 2 || exist('st3ddir_opt','file') ~= 2
    fprintf('getStressletFlow self-test skipped: stfmm3d or st3ddir_opt not found.\n');
    return;
end

rng(2);
Nt = 25;
Ns = 20;

rvec_out = rand(Nt,3);
rvec_in = rand(Ns,3);
nn = rand(Ns,3);
nn = nn ./ vecnorm(nn,2,2);
strslet = rand(3*Ns,1);

vars.fmm = 0;
u_direct = getStressletFlow(rvec_in, rvec_out, nn, strslet, Ns, vars);

vars.fmm = 1;
vars.fmm_tol = 1e-10;
u_fmm = getStressletFlow(rvec_in, rvec_out, nn, strslet, Ns, vars);

rel_err = norm(u_fmm - u_direct) / norm(u_direct);
fprintf('getStressletFlow self-test: rel err (fmm vs direct) = %.3e\n', rel_err);

end
