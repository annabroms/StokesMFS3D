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
%       rin    - 3*N*P x 1 matrix of all source (proxy) point positions.
%       rout   - 3*M*P x 1 matrix of all target (collocation) point positions.
%       q      - P x 3 array of particle centers.
%       Uii    - Cell array {U} containing left preconditioner matrix from
%               one-body SVD for body i in cell i
%       Yii    - Cell array {Y} containing right preconditioner matrix from one-body SVD.
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
%       - getFlow, rotate_vector, 
%
%  Anna Broms, June 12, 2025

P = size(q,1); %number of particles

M = size(rout,1)/P; %points per particle on outer grids
N = size(rin,1)/P; %points per particle on proxy surface

%For now, we assume everyone has the same shape
U = Uii{1};
Y = Yii{1};


%% First, map density for all particles lambda <- mu, using blocks for
%pseudoinverse
lambda_stokes = zeros(3*P*N,1);

if vars.profile
    memorygraph('label','apply precond in matvec');
end

for i = 1:P
%Precomputation is done only for a single paricle (all are assumed to have
%the same shape). Otherwise we would need to retrieve self evaluation
%blocks U{i} and Y{i} here.
    if resistance_flag %solving a resistance problem
        if vars.ellipsoid %rotations needed
            Ri = R{i};
        
            step0 = rotate_vector(mu((i-1)*3*M+1:M*i*3),Ri');
            step1 = U*step0;
            lambda_i = rotate_vector(Y*step1,Ri);
        else
            %spheres: no rotations needed
            step1 = U*mu((i-1)*3*M+1:M*i*3);
            lambda_i = Y*step1;
        end
    else %solving a mobility problem
        if vars.ellipsoid 
       
            Ri = R{i};
        
            step0 = rotate_vector(mu((i-1)*3*M+1:M*i*3),Ri');
            step1 = U*step0; 
            tau_mapped1 = Y*step1;
            tau_mapped2 = rotate_vector(tau_mapped1,Ri);
    
            lambda_i = tau_mapped2-rotate_vector(L*tau_mapped1,Ri);
            
        else
            step1 = U*mu((i-1)*3*M+1:M*i*3);
            tau_mapped = Y*step1; 
            lambda_i = tau_mapped-L*tau_mapped; 
        end
    end

    lambda_stokes(3*(i-1)*N+1:3*i*N) = lambda_i(1:3*N);

end


if vars.profile
    memorygraph('label','compute FMM');
end

%% Do one call to FMM (or direct evaluation) with all sources and targets
 %points
res = getFlow(lambda_stokes,rin,rout,vars);

%% Adjust to obtain identity blocks on diagonal of system matrix
res = res+mu; 

vars.fmm = 0; %a small block, fmm not needed. Maybe better compute full matrix vector product?
%Correct self evaluation: subtract self-interaction 
for i = 1:P

    rin_i = rin(N*(i-1)+1:N*i,:); %sources on body i    
    rows_i = (i-1)*M+1:i*M;
    targ = rout(rows_i,:);  %collocation points on body i
    u_self = getFlow(lambda_stokes(3*(i-1)*N+1:3*i*N),rin_i,targ,vars); %self-interaction
    res((i-1)*3*M+1:i*3*M) = res((i-1)*3*M+1:i*3*M)-u_self;

end

end
