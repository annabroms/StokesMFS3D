function res = matvecStokesMFS_DLP(mu, rin, rout, q, Uii, Yii, TGblock, nn, vars, resistance_flag,R,L)
%MATVECSTOKESMFS_DLP Matrix-vector product for basic but non-standard
%Stokes MFS. This is without image enhancement, but using a combined
%ansatz, TG, with T a DLP (stresslet) mapping and G are the standard
%stokeslets. 
%
%   res = MATVECSTOKESMFS_DLP(mu, rin, rout, q, Uii, Yii, TSblock, nn, vars, resistance_flag, R, L)
%
%   Computes the matrix-vector product A*mu for the linear system arising
%   in a 1-body precomputed Stokes problem solved via the Method of Fundamental Solutions 
%   (MFS) with the TG ansatz (composed DLP and SLP).
%   Used as a callback in GMRES, both for the resistance and mobility problems
%
%   INPUTS:
%       mu    - 3*N*P x 1 vector of boundary data at proxy points on all particles.
%       rin    - N*P x 3 matrix of all source (proxy) point positions.
%       rout   - M*P x 3 matrix of all target (collocation) point positions.
%       q      - P x 3 array of particle centers.
%       Uii    - Cell array {U} containing left preconditioner matrix from
%               one-body SVD for body i in cell i
%       Yii    - Cell array {Y} containing right preconditioner matrix from one-body SVD.
%       TSblock - Precomputed body-frame T*S matrix (DLP times SLP) for a
%               single body. For ellipsoids this must be formed in the body
%               frame: TSblock = stokes_DLP_mat(rout_body,...) * stokes_SLP_mat(rin_body,rout_body).
%       vars   - Struct with solver settings and flags:
%                - vars.fmm: if true, use FMM3D to evaluate flow.
%                - vars.profile: if true, calls memorygraph profiling tool.
%                - vars.ellipsoid: if true, need rotation matrices in R.
%       resistance_flag - boolean to determine whether a resistance or mobility
%                   problem is solved
%       R      - (Optional) Cell array of 3x3 rotation matrices for ellipsoidal particles.
%       L      - Single body projection matrix (only needed for mobility)
%      
%
%   OUTPUT:
%       res    - Resulting 3*N*P x 1 velocity vector corresponding to A*mu.
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
%       - getStokesletFlow, getStressletFlow, rotate_vector, 
%
%  Anna Broms, November 12, 2025, updated April 25, 2026

P = size(q,1); %number of particles 
N = vars.N; %points per particle on proxy surface

%For now, we assume everyone has the same shape
U = Uii{1};
Y = Yii{1};


%% First, map density for all particles lambda <- mu, using blocks for
%pseudoinverse
lambda_stokes = zeros(3*P*N,1);

if vars.profile
    memorygraph('label','apply precond in matvec');
end

if ~vars.ellipsoid
    step1 = U * reshape(mu, 3*N, P);
    tau_body = Y * step1;

    if resistance_flag
        lambda_stokes = tau_body(:);
    else
        lambda_stokes = reshape(tau_body - L * tau_body, [], 1);
    end

    % Reference loop kept for comparison:
    % for i = 1:P
    %     if resistance_flag
    %         step1_i = U*mu((i-1)*3*N+1:i*3*N);
    %         lambda_i = Y*step1_i;
    %     else
    %         step1_i = U*mu((i-1)*3*N+1:i*3*N);
    %         tau_mapped = Y*step1_i;
    %         lambda_i = tau_mapped-L*tau_mapped;
    %     end
    %     lambda_stokes(3*(i-1)*N+1:3*i*N) = lambda_i(1:3*N);
    % end
else
    for i = 1:P
    %Precomputation is done only for a single particle (all are assumed to have
    %the same shape). Otherwise we would need to retrieve self evaluation
    %blocks U{i} and Y{i} here.
        Ri = R{i};
        rows_n = (i-1)*3*N+1:i*3*N;

        if resistance_flag %solving a resistance problem
            step0 = rotate_vector(mu(rows_n),Ri');
            step1 = U*step0;
            lambda_i = rotate_vector(Y*step1,Ri);
        else %solving a mobility problem
            step0 = rotate_vector(mu(rows_n),Ri');
            step1 = U*step0;
            tau_mapped1 = Y*step1;
            tau_mapped2 = rotate_vector(tau_mapped1,Ri);

            lambda_i = tau_mapped2-rotate_vector(L*tau_mapped1,Ri);
        end

        lambda_stokes(rows_n) = lambda_i(1:3*N);
    end
end


if vars.profile
    memorygraph('label','compute FMM');
end

%% Do one call to FMM (or direct evaluation) with all sources and targets
 %points
res_stokes = getStokesletFlow(lambda_stokes,rin,rout,vars);

% apply the matrix T from the left, i.e. evaluate a stresslet flow
res = getStressletFlow(rout,rin,nn,res_stokes,vars.M*P,vars);

%% Adjust to obtain identity blocks on diagonal of system matrix
res = res+mu; 

vars.fmm = 0; %a small block, fmm not needed. Maybe better compute full matrix vector product?
%Correct self evaluation: subtract self-interaction 

if vars.ellipsoid
    Rpages = cat(3, R{:});
    lambda_body = rotate_vector_pages(lambda_stokes, permute(Rpages, [2 1 3]), N);
    u_body = TGblock * reshape(lambda_body, 3*N, P);
    u_self = rotate_vector_pages(u_body, Rpages, N);
    res = res - u_self(:);
else
    u_self = TGblock * reshape(lambda_stokes, 3*N, P);
    res = res - u_self(:);
end

% Reference loop kept for comparison:
% for i = 1:P
%     vec_n = (i-1)*3*N + (1:3*N);
%     if vars.ellipsoid
%         Ri = R{i};
%         lambda_body_i = rotate_vector(lambda_stokes(vec_n), Ri');
%         u_self_i = rotate_vector(TGblock * lambda_body_i, Ri);
%     else
%         u_self_i = TGblock * lambda_stokes(vec_n);
%     end
%     res(vec_n) = res(vec_n) - u_self_i;
% end


end
