function res = getFlow(tau_stokes, rin, rout, vars)
%GETFLOW Evaluate Stokes flow velocity from Stokeslet sources 
%
%   res = GETFLOW(tau_stokes, rin, rout, vars)
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
%                    - vars.eps  : Precision parameter used by FMM3D.
%
%   OUTPUT:
%       res        : 3M x 1 vector of fluid velocity at each target point, stacked
%                    as [u1; u2; ...; uM], where each ui is a 3D vector.
%
%   NOTES:
%       - If vars.fmm is true, the function uses the FMM3D interface `stfmm3d`
%         (requires FMM3D to be installed and compiled).
%       - Otherwise, it falls back to direct evaluation via a mex interface to SE0P_Stokeslet_direct_full_ext_mex.
%       - Modify this function for solving MFS with other fast summation
%         technique or implementation of direction summation. 
%
%   DEPENDENCIES:
%       - Accelerated evaluation: FMM3D (https://github.com/flatironinstitute/fmm3d)
%       - Direct evaluation: SE0P_Stokeslet_direct_full_ext_mex
%       (https://github.com/annabroms/Stokes_Direct -- precompiled binary exists)
%
% Anna Broms June 12, 2025

if vars.fmm 
    % -------- Use FMM3D for fast evaluation --------
    nd = 1;                            % Number of densities per target (velocity only)
    srcinfo.nd = nd;                  % Required field for FMM3D
    ifppreg = 0;                      % Don't compute potential at source
    ifppregtarg = 1;                  % Evaluate potential at target
    
    srcinfo.sources = rin';           % Transpose: 3 x N
    srcinfo.stoklet = reshape(tau_stokes, 3, []); % Format: 3 x N
    
    targ = rout';                     % Transpose target: 3 x M
    eps = vars.eps;                   % FMM precision (e.g., 1e-6)
    
    % Call FMM3D driver for Stokeslet potentials
    U = stfmm3d(eps, srcinfo, ifppreg, targ, ifppregtarg); 
    
    % Extract and reshape velocity result
    res = U.pottarg(:); 

    clear U srcinfo;

else
    % -------- Use direct evaluation --------
    targ = rout; 
    srcinfo.stoklet = reshape(tau_stokes, 3, []); % Format: 3 x N

    % Evaluate Stokeslet flow at targets via direct sum
    U = SE0P_Stokeslet_direct_full_ext_mex(rin, srcinfo.stoklet', struct('eval_ext_x', targ));
    
    U = U';              % Transpose to match output stacking
    res = 1/(8*pi) * U(:); % Apply Stokeslet prefactor
end



end