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
%                  - vars.eps : Precision parameter for FMM3D.
%
%   OUTPUT:
%       u        : 3*Nt x 1 stacked velocity at targets.
%
%   NOTES:
%       - If vars.fmm is true, uses stfmm3d stresslet evaluation.
%       - Otherwise, uses SE0P_Stresslet_direct_full_ext_mex (fast direct sum).
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


if vars.fmm
    ifppreg = 0;      % no eval at sources
    ifppregtarg = 1;  % just vel out
    
    targ = rvec_out';
    
    srcinfo.sources = rvec_in';
    srcinfo.stoklet = zeros(3,Ns); 
    
    srcinfo.strslet = reshape(strslet, 3, []); % Format: 3 x Ns   
    srcinfo.strsvec = nn';
    
    U = stfmm3d(vars.eps,srcinfo,ifppreg,targ,ifppregtarg);
    u = U.pottarg;
    u = -u(:); %there is a sign difference!
else
    targ = rvec_out; 
    srcinfo.strslet = reshape(strslet, 3, []); % Format: 3 x N

    % Evaluate stresslet flow at targets via direct sum
    U = SE0P_Stresslet_direct_full_ext_mex(rvec_in, srcinfo.strslet', nn, struct('eval_ext_x', targ));
    
    U = U';              % Transpose to match output stacking
    u = 1/(8*pi) * U(:); % Apply Stokeslet prefactor
end

end

function self_test()

if exist('stfmm3d','file') ~= 2
    fprintf('getStressletFlow self-test skipped: stfmm3d not found.\n');
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
vars.eps = 1e-10;
u_fmm = getStressletFlow(rvec_in, rvec_out, nn, strslet, Ns, vars);

rel_err = norm(u_fmm - u_direct) / norm(u_direct);
fprintf('getStressletFlow self-test: rel err (fmm vs direct) = %.3e\n', rel_err);

end
