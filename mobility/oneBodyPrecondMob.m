function [Y, UU, LL, Kin, Kout] = oneBodyPrecondMob(rin, rout, q)
%oneBodyPrecondMob Determines pseudoinverse factors for the single body MFS
% Stokes mobility system matrix.
%
%   [Y, UU, LL, Kin, Kout] = oneBodyPrecondMob(rin, rout, q)
%  
%   The function builds and factorizes a one-body matrix that does not
%   contribute to net force and torque and that utilizes an unused subspace (the range of Kin) to
%   express rigid body motions in terms of the unknown source vector,
%
%   INPUTS:
%       rin  - N x 3 matrix of proxy source points.
%       rout - M x 3 matrix of collocation points.
%       q    - (Optional) 1 x 3 vector for the particle center. Default: [0 0 0].
%
%   OUTPUTS:
%       U - Matrix of left singular vectors of the 1-body system matrix 
%           corresponding to retained singular values
%       Y - Product VS⁺, where:
%           - S⁺ is a diagonal matrix with entries 1/σ for retained
%           singular values of the 1-body system matrix
%           - V contains the corresponding right singular vectors
%       LL    - Force/torque projection matrix: LL = Kin * inv(Kin^T Kin) * Kin^T.
%       Kin   - Mapping from source strength to net force and torque.
%       Kout  - Mapping from RBM velocity to velocity on boundary.
%
%   DEPENDENCIES:
%       - generate_stokes_mat, getKmat, getPseudoFactors
%
% Anna Broms, June 12, 2025

if nargin < 3
    q = [0 0 0];
end

% Build MFS and rigid-body matrices
S    = generate_stokes_mat(rin, rout);
Kin  = getKmat(rin, q);
Kout = getKmat(rout, q);

% Construct the force/torque projection matrix
LL = Kin * ((Kin' * Kin) \ Kin');

% Build one-body operator that does not contribute to force/torque. The
% last term takes care of velocity constraint.
A = S * (eye(size(LL)) - LL) + Kout * Kin';

% Compute pseudoinverse factors 
tol = 1e-15;
visualise = 0; 
[Y, UU] = getPseudoFactors(A,tol,visualise);
UU = UU';

end

