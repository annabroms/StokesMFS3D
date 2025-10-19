function res = applyPmat(vel,rvec_in,rvec_out,Rinv,Zi,Yi,R,opt)
%APPLYPMAT  Apply long-range preconditioning projection matrix to a velocity field.
%
%   res = APPLYPMAT(vel, rvec_in, rvec_out, Rinv, Zi, Yi, opt)
%
%   Applies the projector 
%
%       P = I - G * (Z * (R \ I) * Y')
%
%   to the input velocity vector. This corresponds to projecting out
%   contributions from the coarse subspace used in long-range preconditioning.
%
%   INPUTS:
%     vel       - 3PM×1 velocity vector at the M boundary points of P bodies.
%     rvec_in   - PNx3 array of source points for the projection flow.
%     rvec_out  - PMx3 array of target points where the projection flow is evaluated.
%     Zi         - (3N×k) matrix mapping coarse coefficients to proxy forces for a single body (diagonal block in Z). 
%                 k is set by the number of coarse basis functions per body
%     Yi         - (3M×k) matrix such that Y' maps surface flow to coarse space for a single body.
%     Rinv      - (Pk×Pk) inverse of coarse interaction matrix R.  
%     opt       - Struct with fmm and lr flags etc
%
%   OUTPUT:
%     res       - Projected velocity vector res = P vel
%                 where lambda are the coarse source strengths.
%
%   NOTES:
%     - The projector removes the coarse flow contribution defined by the span
%       of Z, and used in the long-range preconditioner.
%     - The matrices Z and Y are sorted as [x y z x y z ... ]
%     - Currently implemented without image system contributions.
%
%   See also: applyQmat, getFlow
%
% Anna Broms Oct 14, 2025

lambda = getCoarseSource(vel,Rinv,Zi,Yi,R,opt);

%compute velocities using these source strengths

proj_vel = getFlow(lambda, rvec_in, rvec_out, opt); 

res = vel-proj_vel;


end