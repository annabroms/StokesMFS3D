function [Y, UU, LL, Kin, Kout] = getSL(rin, rout, q)
%GETSL Preconditioner factors for the Stokes mobility matrix (single body).
%
%   [Y, UU, LL, Kin, Kout] = GETSL(rin, rout, q)
%
%   Constructs preconditioner and pseudoinverse factors for a single particle
%   in the Stokes mobility problem using the Method of Fundamental Solutions (MFS).
%   The function builds and factorizes a one-body matrix that enforces force/torque
%   constraints while ensuring minimal-norm source densities.
%
%   INPUTS:
%       rin  - N x 3 matrix of proxy source points.
%       rout - M x 3 matrix of collocation points.
%       q    - (Optional) 1 x 3 vector for the particle center. Default: [0 0 0].
%
%   OUTPUTS:
%       Y     - Right-side pseudoinverse factor from getPseudoFactors.
%       UU    - Left-side projection factor from getPseudoFactors.
%       LL    - Force/torque projection matrix: LL = Kin * inv(Kinᵗ Kin) * Kinᵗ.
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

% Build one-body operator
A = S * (eye(size(LL)) - LL) + Kout * Kin';

% Compute pseudoinverse factors 
tol = 1e-15;
visualise = 0; 
[Y, UU] = getPseudoFactors(A,tol,visualise);

end

