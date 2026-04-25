function res = matvecStokesMFS(mu, rin, rout, q, Uii, Yii, vars, resistance_flag,R,L)
%MATVECSTOKESMFS Matrix-vector product for basic Stokes MFS (without image enhancement) 
%
%   res = MATVECSTOKESMFS(mu, rin, rout, q, Uii, Yii, vars, R, resistance_flag,R,L)
%
%   Computes the matrix-vector product A*mu for the linear system arising
%   in a 1-body precomputed Stokes problem solved via the Method of Fundamental Solutions (MFS).
%   Used as a callback in GMRES, both for the resistance and mobility
%   problems
%
%   INPUTS:
%       mu    - 3*M*P x 1 vector of boundary data at collocation points on all particles.
%       rin    - N*P x 3 matrix of all source (proxy) point positions.
%       rout   - M*P x 3 matrix of all target (collocation) point positions.
%       q      - P x 3 array of particle centers.
%       Uii    - Cell array {U} containing left preconditioner matrix from
%               one-body SVD for body i in cell i
%       Yii    - Cell array {Y} containing right preconditioner matrix from one-body SVD.
%       vars   - Struct with solver settings and flags:
%                - vars.fmm: if true, use FMM3D to evaluate flow.
%                - vars.profile: if true, calls memorygraph profiling tool.
%                - vars.ellipsoid: if true, need rotation matrices in R.
%                - vars.Sblock: cached one-body self block.
%       resistance_flag - boolean to determine whether a resistance or mobility
%                   problem is solved
%       R      - (Optional) Cell array of 3x3 rotation matrices for ellipsoidal particles.
%       L      - Single body projection matrix (only needed for mobility)
%      
%
%   OUTPUT:
%       res    - Resulting 3*M*P x 1 velocity vector corresponding to A*mu.
%
%   METHOD:
%       1. Applies preconditioner mapping: source density from boundary
%       data: lambda <- mu. If mobility, also computes lambda <- (I-L)lambda
%       2. Uses FMM or direct sum to evaluate flow at all collocation
%       points due to all sources.
%       3. Correct identity blocks: Adds mu and removes diagonal
%       self-interaction contributions to improve stability and conditioning
%
%   ASSUMPTIONS:
%       - All particles have the same geometry and source setup, up to rotations (however
%       easy to modify). 
%
%   DEPENDENCIES:
%       - getStokesletFlow, rotate_vector, 
%
%  Anna Broms, June 12, 2025

P = size(q,1); %number of particles

M = vars.M; %points per particle on outer grids
N = vars.N; %points per particle on proxy surface

%For now, we assume everyone has the same shape
U = Uii{1};
Y = Yii{1};


%% First, map density for all particles lambda <- mu, using blocks for
%pseudoinverse

if vars.profile
    memorygraph('label','apply precond in matvec');
end

if vars.ellipsoid
    Rpages = cat(3, R{:});
    mu_body = rotate_vector_pages(mu, permute(Rpages, [2 1 3]), M);
    step1 = U * reshape(mu_body, 3*M, P);
    tau_body = Y * step1;
    if resistance_flag
        lambda_stokes = rotate_vector_pages(tau_body, Rpages, N);
    else
        lambda_stokes = rotate_vector_pages(tau_body - L * tau_body, Rpages, N);
    end
    lambda_stokes = lambda_stokes(:);
else
    step1 = U * reshape(mu, 3*M, P);
    tau_body = Y * step1;
    if resistance_flag
        lambda_stokes = tau_body(:);
    else
        lambda_stokes = reshape(tau_body - L * tau_body, [], 1);
    end
end

% Reference loop kept for comparison:
% lambda_stokes = zeros(3*P*N,1);
% for i = 1:P
%     rows_m = (i-1)*3*M + (1:3*M);
%     rows_n = (i-1)*3*N + (1:3*N);
%     if vars.ellipsoid
%         Ri = R{i};
%         step0 = rotate_vector(mu(rows_m), Ri');
%         tau_i = Y * (U * step0);
%         if resistance_flag
%             lambda_i = rotate_vector(tau_i, Ri);
%         else
%             lambda_i = rotate_vector(tau_i - L*tau_i, Ri);
%         end
%     else
%         tau_i = Y * (U * mu(rows_m));
%         if resistance_flag
%             lambda_i = tau_i;
%         else
%             lambda_i = tau_i - L*tau_i;
%         end
%     end
%     lambda_stokes(rows_n) = lambda_i(1:3*N);
% end


if vars.profile
    memorygraph('label','compute FMM');
end

%% Do one call to FMM (or direct evaluation) with all sources and targets
 %points
res = getStokesletFlow(lambda_stokes,rin,rout,vars);

%% Adjust to obtain identity blocks on diagonal of system matrix
res = res+mu; 

vars.fmm = 0; %a small block, fmm not needed. Maybe better compute full matrix vector product?
%Correct self evaluation: subtract self-interaction 

if vars.ellipsoid
    Rpages = cat(3, R{:});
    lambda_body = rotate_vector_pages(lambda_stokes, permute(Rpages, [2 1 3]), N);
    u_body = vars.Sblock * reshape(lambda_body, 3*N, P);
    u_self = rotate_vector_pages(u_body, Rpages, M);
    res = res - u_self(:);
else
    u_self = vars.Sblock * reshape(lambda_stokes, 3*N, P);
    res = res - u_self(:);
end

% Reference loop kept for comparison:
% for i = 1:P
%     rin_i = rin(N*(i-1)+1:N*i,:);
%     rows_i = (i-1)*M+1:i*M;
%     targ = rout(rows_i,:);
%     u_self_i = getStokesletFlow(lambda_stokes(3*(i-1)*N+1:i*3*N),rin_i,targ,vars);
%     res((i-1)*3*M+1:i*3*M) = res((i-1)*3*M+1:i*3*M)-u_self_i;
% end

end
