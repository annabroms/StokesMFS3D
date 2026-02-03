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
%                    - vars.eps  : Precision parameter used by FMM3D.
%
%   OUTPUT:
%       res        : 3M x 1 vector of traction at each target point, stacked
%                    as [t1; t2; ...; tM], where each ti is a 3D vector.
%
%   NOTES:
%       - If vars.fmm is true, the function uses the FMM3D interface `stfmm3d`
%         (requires FMM3D to be installed and compiled).
%       - Otherwise, it falls back to direct evaluation via a mex interface to SE0P_Stokestraction_direct_full_ext_mex.
%       - Modify this function for solving MFS with other fast summation
%         technique or implementation of direction summation. 
%
%   DEPENDENCIES:
%       - Accelerated evaluation: FMM3D (https://github.com/flatironinstitute/fmm3d)
%       - Direct evaluation: SE0P_Stokestraction_direct_full_ext_mex
%       (https://github.com/annabroms/Stokes_Direct -- precompiled binary exists)
%
%   See also: getTractionDirectSLP
%
% Anna Broms Feb 3, 2025

if nargin < 1
    self_test;
    return;
end

if vars.fmm 
    % -------- Use FMM3D for fast evaluation --------
    error('To be included')

else
    % -------- Use direct evaluation --------
    targ = rout;
    srcinfo.stoklet = reshape(tau_stokes, 3, []); % 3 x N

    % Traction at targets via direct sum (mex)
    
    U = SE0P_Stokestraction_direct_full_ext_mex( ...
        rin, srcinfo.stoklet', nn, struct('eval_ext_x', targ));
 

    U = U';                 % 3 x M
    res = (1/(8*pi)) * U(:); % Apply traction prefactor (mu=1)

end

function self_test()

rng(1);

Nsrc = 30;
Ntar = 25;

rvec_in = rand(Nsrc,3);
rvec_out = rand(Ntar,3);
nn = rand(Ntar,3);
nn = nn ./ vecnorm(nn,2,2); % normalize

mu = rand(3*Nsrc,1);

vars.fmm = 0;

res_fast = getTractionFast(mu, rvec_in, rvec_out, nn, vars);
res_direct = getTractionDirectSLP(rvec_in, rvec_out, nn, mu);

T = getTraction(rvec_in, rvec_out, nn);
res_mat = T * mu;

rel_err_fast = norm(res_fast - res_direct) / norm(res_direct);
rel_err_mat = norm(res_mat - res_direct) / norm(res_direct);

fprintf('Traction fast vs direct: rel err = %.3e\n', rel_err_fast);
fprintf('Traction mat  vs direct: rel err = %.3e\n', rel_err_mat);

end


end
