function lambda0 = getRecompletionSource(F, T, Kin)
%GETRECOMPLETIONSOURCE Construct recompletion source lambda0 for a particle in Stokes flow.
%
%   lambda0 = GETRECOMPLETIONSOURCE(F, T, Kin)
%
%   Computes the Stokeslet source density lambda0 on a proxy surface (inner
%   surface) such that it reproduces the specified net force F and
%   torque T on a single rigid particle in Stokes flow. 
%
%   INPUTS:
%       F     - 3 x 1 force vector applied to the particle.
%       T     - 3 x 1 torque vector applied to the particle.
%       Kin   - 6 x 3N matrix that maps the source density lambda to the
%               net force and torque: [F; T] = Kin' * lambda.
%
%   OUTPUT:
%       lambda0 - 3N x 1 vector of source strengths on the proxy surface
%                 that yields the specified F and T.
%
%   METHOD:
%       - Solve the least-squares problem:
%           min || lambda ||²   subject to   Kin' * lambda = [F; T]
%       - The minimal-norm solution is given by:
%           lambda0 = Kin * ((Kin' * Kin) \ [F; T])
%
%   USE CASE:
%       - Used to "complete" the right-hand side in mobility problems to
%         ensure prescribed forces/torques are matched by the flow.
%
%   SEE ALSO:
%       solve_mobility
%
% Anna Broms, June 12, 2025

D = Kin'*Kin;
ab = D\[F;T];
lambda0 = Kin*ab;

end