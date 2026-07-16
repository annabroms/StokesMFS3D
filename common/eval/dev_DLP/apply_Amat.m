function y = apply_Amat(x,rin,rout,nout,wout,vars)
%APPLY_AMAT Apply the matrix-free weighted DLP/SLP operator A = T*W*S.
%
%   y = APPLY_AMAT(x,rin,rout,nout,wout,vars)
%
%   Applies the composed operator
%
%       A = T(rout -> rin) * W * S(rin -> rout),
%
%   without assembling its dense matrix. Here S is the Stokes single-layer
%   (Stokeslet) operator from the inner proxy nodes rin to the outer
%   quadrature nodes rout, W applies the scalar quadrature weights wout to
%   all three vector components, and T is the double-layer (stresslet)
%   operator from rout back to rin with directions nout.
%
%   INPUTS:
%       x     - 3*N*P-by-1 stacked vector density at the N inner proxy nodes for the P bodies.
%       rin   - NP-by-3 inner proxy/source-node coordinates.
%       rout  - MP-by-3 outer quadrature-node coordinates.
%       nout  - MP-by-3 stresslet directions, normally the outward normals
%               at ROUT.
%       wout  - MP-vector of scalar quadrature weights at ROUT.
%       vars  - Evaluation options passed to GETSTOKESLETFLOW and
%               GETSTRESSLETFLOW. Relevant fields include fmm and eps.
%
%   OUTPUT:
%       y     - 3*N*P-by-1 stacked vector A*x at the inner nodes.
%
%   The vector layout is [v1_x; v1_y; v1_z; v2_x; ...].
%
%   See also APPLY_SYM_AMAT, APPLY_QUAD_WEIGHTS, GETSTOKESLETFLOW,
%   GETSTRESSLETFLOW.


u = getStokesletFlow(x,rin,rout,vars);
u = apply_quad_weights(u,wout);
% The magnitudes fed in to the stresslet are smaller. Better results with
% more strict tolerance here
vars_stokes = vars;
vars_stokes.fmm_tol = vars.fmm_tol*1e-2;
y = getStressletFlow(rout,rin,nout,u,numel(wout),vars_stokes);

end
