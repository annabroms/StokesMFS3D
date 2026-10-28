function res = wall_images_matvec(x,problem)
%WALL_IMAGES_MATVEC Apply the square one-body-preconditioned wall operator.
% x and res have 3*sum(M) interleaved collocation entries. The factors
% Y{k} (3*N_total(k)-by-rank) and UU{k} (rank-by-3*M(k)) map x to source
% coefficients. U is never used for singular vectors: U denotes rigid motion.
%
% Resistance: evaluate the flow from these source coefficients.
% Mobility: first remove their force/torque using the six-column completion
% map. This is the main mobility solver's LL projection, applied without
% constructing a full source-by-source projector.
%
% Subtract the factored self response UU'*(UU*x) and add x, so each diagonal
% block becomes identity. For mobility, also restore Kout*Kin'*lambda_i,
% which belongs to the factored self matrix but not the physical flow.
% Thus all cross-body interactions retain the complete wall kernels.

lambda = zeros(problem.sourceOffsets(end),1);
res = x;
mobility = strcmp(problem.mode,'mobility');
for k = 1:problem.P
    rows = problem.targetOffsets(k)+(1:3*problem.M(k));
    cols = problem.sourceOffsets(k)+(1:3*problem.N_total(k));
    body = problem.body{k};
    % Apply the stored pseudoinverse in two stable factor multiplications.
    step1 = problem.UU{k}*x(rows);
    lambda_i = problem.Y{k}*step1;
    res(rows) = res(rows)-problem.UU{k}.'*step1;
    if mobility
        moments = body.Kin.'*lambda_i;
        res(rows) = res(rows)+body.Kout*moments;
        lambda_i = lambda_i-body.completion*moments;
    end
    lambda(cols) = lambda_i;
end
res = res+apply_wall_images_velocity(lambda,problem);
end
