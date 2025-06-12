function [Y, U] = oneBodyPrecondRes(rin, rout)
%ONEBODYPRECONDRES Construct preconditioner for one-body Stokes resistance problem.
%
%   [Y, UU] = ONEBODYPRECONDRES(rin, rout)
%
%   Constructs pseudoinverse factors for a single body's Stokes resistance problem,
%   using the Method of Fundamental Solutions (MFS). The function builds the 
%   dense single-layer Stokeslet matrix from sources to collocation points, then 
%   computes a truncated SVD to extract pseudoinverse factors, to later
%   enable a backward-stable apply.
%
%   INPUTS:
%       rin  - N x 3 matrix of source (proxy) point locations.
%       rout - M x 3 matrix of collocation (boundary) point locations.
%
%   OUTPUTS:
%   U - Matrix of left singular vectors corresponding to retained singular values
%   Y - Product VS⁺, where:
%         - S⁺ is a diagonal matrix with entries 1/σ for retained singular values
%         - V contains the corresponding right singular vectors
%
%   USAGE:
%       - This function is intended for resistance problems involving a single 
%         particle.
%       - The outputs (Y, UU) can be reused for all identical particles in a 
%         many-body resistance solve.
%
%   NOTES:
%       - Singular values below tol * max(Σ) are truncated, with a
%         tol=1e-15 as a default.
%       - The condition number of S and its spectrum are not returned but could
%         be visualized for diagnostics with the flag visualise = 1. 
%
%   DEPENDENCIES:
%       - generate_stokes_mat.m, getPseudoFactors.m
%
%  Anna Broms, June 12, 2025

tol = 1e-15;
visualise = 0; 

S = generate_stokes_mat(rin,rout);

[Y, U] = getPseudoFactors(S, tol, visualise);
U = U'; 

end
