%AUG05_DEFLATION_PROJECTION_CHECK Verify particlewise near-null projection.
%
% This deterministic experiment checks the shared particlewise projectors
% against explicit dense matrices.  All unit-sphere grids are spherical
% designs; no ellipsoid quadrature is used.

clear;
clc;
rng(1);

%% Two identical proxy spheres, with a nontrivial rotation on particle two
Rp = 0.30;
[source_unit,~] = get_sphdesign(100);
N = size(source_unit,1);
rin_reference = Rp * source_unit;

q = [0 0 0; 0 0 3];
alpha = 0.47;
beta = -0.61;
Rz = [cos(alpha) -sin(alpha) 0; ...
      sin(alpha)  cos(alpha) 0; ...
      0           0          1];
Ry = [cos(beta) 0 sin(beta); ...
      0         1 0; ...
     -sin(beta) 0 cos(beta)];
R = {eye(3); Rz*Ry};
P = numel(R);

rvec_in = zeros(P*N,3);
for k = 1:P
    rows = (k-1)*N + (1:N);
    rvec_in(rows,:) = q(k,:) + (R{k} * rin_reference')';
end

opt.deflation_check_nmax = 150;
opt.deflation_verify_nmax = 180;
opt.null_subspace_tol = 1e-6;
opt.null_field_tol = 1e-6;
opt.constraint_tol = 1e-12;
opt.normal_correlation_min = 0.9;

projection = build_deflation_projection(q,rvec_in,opt,R);
assert(projection.info.accepted, ...
    'Projection mode rejected: %s',projection.info.rejection_reason);
assert(projection.info.source_nodes < projection.info.actual_check_nodes && ...
       projection.info.actual_check_nodes < projection.info.actual_verify_nodes, ...
    'The available spherical designs do not give N < Mcheck < Mverify.');

%% Explicit matrices are references for this experiment only
ndof = projection.ndof;
L = zeros(ndof);
mode_matrix = zeros(ndof,P);
for k = 1:P
    rows = (k-1)*3*N + (1:3*N);
    Kk = projection.Kblocks{k};
    L(rows,rows) = Kk * ((Kk' * Kk) \ Kk');
    mode_matrix(rows,k) = projection.modes{k};
end
Nprojector = mode_matrix * mode_matrix';
Q = eye(ndof) - L - Nprojector;

x = randn(ndof,4);
parts = {'L','N','P','Q'};
explicit_actions = {L*x, Nprojector*x, (eye(ndof)-L)*x, Q*x};
action_errors = zeros(numel(parts),1);
for j = 1:numel(parts)
    y = projection.apply(x,parts{j});
    action_errors(j) = norm(y-explicit_actions{j},'fro') / ...
        max(norm(explicit_actions{j},'fro'),eps);
end

Kglobal = getKmat(rvec_in,q);
symmetry_errors = [relative_matrix_error(L,L'); ...
                   relative_matrix_error(Nprojector,Nprojector'); ...
                   relative_matrix_error(Q,Q')];
idempotence_errors = [relative_matrix_error(L*L,L); ...
                      relative_matrix_error(Nprojector*Nprojector,Nprojector); ...
                      relative_matrix_error(Q*Q,Q)];
completeness_error = norm(L+Nprojector+Q-eye(ndof),'fro') / sqrt(ndof);
KQ_error = norm(Kglobal' * Q,'fro') / ...
    max(norm(Kglobal,'fro')*norm(Q,'fro'),eps);
NQ_error = norm(mode_matrix' * Q,'fro') / ...
    max(norm(mode_matrix,'fro')*norm(Q,'fro'),eps);
LQ_error = norm(L * Q,'fro') / max(norm(L,'fro')*norm(Q,'fro'),eps);
mode_orthogonality = norm(mode_matrix' * mode_matrix-eye(P),'fro');
force_torque_error = norm(Kglobal' * mode_matrix,'fro') / ...
    max(norm(mode_matrix,'fro'),eps);

%% Report and enforce the projection contract
info = projection.info;
fprintf('Aug 05 deflation projection check\n');
fprintf('  Nsource / Ncheck / Nverify = %d / %d / %d\n', ...
    info.source_nodes,info.actual_check_nodes,info.actual_verify_nodes);
fprintf('  preliminary rank = %d, used smallest fallback = %d\n', ...
    info.preliminary_rank,info.used_smallest_fallback);
fprintf('  normal correlation = %.16f\n',info.normal_correlation);
fprintf('  weighted field: check %.3e, verify %.3e, rotated max %.3e\n', ...
    info.check_trace_relative,info.verify_trace_relative, ...
    max(info.particle_verify_trace_relative));
fprintf('  reference K''N residual = %.3e\n',info.force_torque_relative);
fprintf('  projector symmetry errors   [L N Q] = [%s]\n', ...
    sprintf(' %.3e',symmetry_errors));
fprintf('  projector idempotence errors[L N Q] = [%s]\n', ...
    sprintf(' %.3e',idempotence_errors));
fprintf('  projector action errors [L N P Q] = [%s]\n', ...
    sprintf(' %.3e',action_errors));
fprintf('  completeness / K''Q / N''Q / LQ = %.3e / %.3e / %.3e / %.3e\n', ...
    completeness_error,KQ_error,NQ_error,LQ_error);
fprintf('  mode orthogonality / force-torque = %.3e / %.3e\n', ...
    mode_orthogonality,force_torque_error);

assert(info.normal_correlation > 0.99, ...
    'The selected source mode is not sufficiently normal-like.');
assert(info.verify_trace_relative < 1e-6 && ...
       max(info.particle_verify_trace_relative) < 1e-6, ...
    'The accepted source mode is not negligible on the verification grid.');
assert(info.force_torque_relative < 1e-12 && force_torque_error < 1e-12, ...
    'The accepted source mode changes force or torque.');
assert(mode_orthogonality < 1e-12, ...
    'Particle modes are not orthonormal.');
assert(max(symmetry_errors) < 1e-12, ...
    'A projection matrix is not symmetric.');
assert(max(idempotence_errors) < 1e-12, ...
    'A projection matrix is not idempotent.');
assert(completeness_error < 1e-12, ...
    'L + N*N'' + Q does not equal the identity.');
assert(max(action_errors) < 1e-12, ...
    'Particlewise projector application disagrees with the explicit matrix.');
assert(max([KQ_error,NQ_error,LQ_error]) < 1e-12, ...
    'The retained projector is not orthogonal to a removed subspace.');

fprintf('  projection check passed.\n');

function err = relative_matrix_error(A,B)
% Return a scale-aware Frobenius error for two matrices.
err = norm(A-B,'fro') / max(norm(B,'fro'),eps);
end
