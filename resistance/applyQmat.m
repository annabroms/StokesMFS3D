function lambda_fine = applyQmat(lambda, rvec_in, rvec_out, Rinv, Z, Y, opt)
%APPLYQMAT  Apply long-range preconditioning projection matrix on sources
%
%   lambda_fine = APPLYQMAT(vel, rvec_in, rvec_out, Rinv, Z, Y, opt)
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
%     Y         - (3PM×3Pk) matrix mapping from surface flow to coarse space.
%     Z         - (3PN×3Pk) matrix mapping from proxy source space to coarse space.
%     opt       - Struct with flow options (e.g., FMM flags, kernel type).
%
%   OUTPUT:
%     lamda_fine       - Projected vector lambda_fine = Q * lambda, with coarse flow contributions removed.
%
%   NOTES:
%     - Unlike APPLYPMAT, this function applies the projection directly on 
%       (i.e., it modifies test functions or right-hand side vectors).
%     - The flow is computed via getFlow using proxy points and Stokes kernel.
%     - Image systems are not included in this implementation.
%
%   See also: applyPmat, getFlow
%
% Anna Broms, Oct 14, 2025

proj_vel = getFlow(lambda, rvec_in, rvec_out, opt); 
lambda_proj = Z*Rinv*(Y'*proj_vel);

lambda_fine = lambda-lambda_proj;


end