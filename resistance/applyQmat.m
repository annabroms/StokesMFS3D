function lambda_fine = applyQmat(lambda, rvec_in, rvec_out, Rinv, Zi,Yi,R, opt)
%APPLYQMAT  Apply long-range preconditioning projection matrix on sources
%
%   lambda_fine = APPLYQMAT(vel, rvec_in, rvec_out, Rinv, Zi, Yi, opt)
%
%   Applies the projection operator 
%
%       Q = I - (Z * (R \ I) * Y') * G
%
%   to the solution vector lambda. 
%
%   INPUTS:
%     lambda    - 3PN×1 vector of proxy source strengths 
%     rvec_in   - PN×3 array of source points for evaluating the proxy flow.
%     rvec_out  - PM×3 array of target surface points where the flow is evaluated.
%     Rinv      - (3Pk×3Pk) inverse of the coarse interaction matrix R,
%                 where k is the number of coarse basis functions per body per coordinate.
%     Yi         - (3M×k) matrix mapping from surface flow to coarse space
%                  for a single body (a diagonal block of Y) 
%     Zi         - (3N×k) matrix mapping from proxy source space to coarse space 
%                   for a single body (a diagonal block of Z).
%     R         - Cell array of P rotation matrices (needed for ellipsoids)
%     opt       - Struct with flow options (e.g., FMM flags, kernel type).
%
%   OUTPUT:
%     lamda_fine       - Projected vector lambda_fine = Q * lambda, with coarse flow contributions removed.
%
%   NOTES:
%     - Unlike APPLYPMAT, this function applies the projection directly on 
%       (i.e., it modifies test functions or right-hand side vectors).
%     - The flow is computed via getStokesletFlow using proxy points and Stokes kernel.
%     - Image systems are not included in this implementation.
%
%   See also: applyPmat, getStokesletFlow
%
% Anna Broms, Oct 14, 2025

proj_vel = getStokesletFlow(lambda, rvec_in, rvec_out, opt); 
lambda_proj = getCoarseSource(proj_vel,Rinv,Zi,Yi,R,opt);

lambda_fine = lambda-lambda_proj;


end