function y = apply_sym_Amat(x,rin,rout,nout,wout,vars)
%APPLY_SYM_AMAT Apply the matrix-free symmetrized DLP/SLP operator.
%
%   y = APPLY_SYM_AMAT(x,rin,rout,nout,wout,vars)
%
%   Applies
%
%       A_sym = (A + A')/2,       A = T*W*S,
%
%   without assembling A or its transpose. The forward action A*x is
%   evaluated by APPLY_AMAT. The transpose action uses the Stokeslet
%   traction operator for T', applies the outer quadrature weights, and
%   then evaluates the reciprocal Stokeslet map from rout to rin.
%
%   INPUTS:
%       x     - 3*N*P-by-1 stacked vector density at the N inner nodes.
%       rin   - NP-by-3 inner proxy/source-node coordinates.
%       rout  - MP-by-3 outer quadrature-node coordinates.
%       nout  - MP-by-3 directions used by the stresslet and traction
%               operators, normally the outward normals at ROUT.
%       wout  - MP-vector of scalar quadrature weights at ROUT.
%       vars  - Evaluation options passed to the Stokes kernel routines.
%               Relevant fields include fmm and eps.
%
%   OUTPUT:
%       y     - 3*N*P-by-1 stacked vector A_sym*x at the inner nodes.
%
%   The traction evaluation used in the transpose action is forced to the
%   direct implementation because the current FMM traction path is slower.
%   Other kernel evaluations continue to follow vars.fmm.
%
%   See also APPLY_AMAT, APPLY_QUAD_WEIGHTS, GETTRACTIONFAST,
%   GETSTOKESLETFLOW.

y = apply_Amat(x,rin,rout,nout,wout,vars);

traction_vars = vars;
traction_vars.fmm = 0;
u = getTractionFast(x,rin,rout,nout,traction_vars);
u = apply_quad_weights(u,wout);
y_transpose = getStokesletFlow(u,rout,rin,vars);

y = 0.5 * (y + y_transpose);

end
