function [U, iters, lambda_norm, uerr] = solve_mobility(q, Fvec, fmm, Rp, N)
%SOLVE_MOBILITY Solve a Stokes mobility problem for spherical particles.
%
%   [U, iters, lambda_norm, uerr] = SOLVE_MOBILITY(q, Fvec, fmm, Rp, N)
%
%   Solves a mobility problem in Stokes flow for spherical particles using
%   the Method of Fundamental Solutions (MFS). Given the applied forces and
%   torques on each particle, the function computes the resulting translational
%   and angular velocities.
%
%   INPUTS:
%       q       - P x 3 matrix of particle centers.
%       Fvec    - 6P x 1 vector of forces and torques on particles.
%                 Format: [F1; T1; F2; T2; ...] with Fi, Ti in R^3.
%       fmm     - Logical flag indicating whether to use FMM3D for fast evaluation.
%       Rp      - (Optional) Radius of the proxy surface around each particle.
%                 Default: 0.68.
%       N       - (Optional) Number of proxy source points per particle.
%                 Default: 700.
%
%   OUTPUTS:
%       U           - 6P x 1 vector of rigid body velocities of the particles.
%                     Format: [u1; omega1; u2; omega2; ...] with each in R^3.
%       iters       - Number of GMRES iterations used to reach the tolerance.
%       lambda_norm - Infinity norm of the final source density vector.
%                     (Large values may signal unfavourable representation
%                     for the problem, e.g. more sources needed or proxy surface 
%                     needs to be closer to the particle surface).
%       uerr        - Maximum relative velocity residual on the particle surfaces.
%                     
%
%   METHOD OVERVIEW:
%       - Initializes proxy and collocation grids on a unit sphere.
%       - Constructs MFS RHS using "completion sources" based on input forces/torques.
%       - Solves the MFS linear system using GMRES with one-body preconditioning.
%       - Maps the solution source vector to rigid body velocity for each particle.
%       - Validates accuracy by checking surface residuals at points different from the collocation nodes.
%
%   DEPENDENCIES:
%       - init_MFS, getDesignGrid, oneBodyPrecondMob, getCompletionSource, getFlow, helsing_gmres,
%         matvecStokesMFS, getKmat
%
%   EXAMPLE USAGE:
%       q = [0 0 0; 2 0 0];
%       Fvec = [1; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0];
%       [U, iters, lambda_norm, uerr] = solve_mobility(q, Fvec, true);
%
% Anna Broms June 13, 2025


P = size(q,1);

% Set default values if Rg and N are not provided
if nargin < 4
    Rp = 0.68; %proxy radius
    N = 700; % approximate number of proxy sources on every particle
elseif nargin < 5
    N = 700;
end

% initialize a bunch of parameters. Do not change if you don't really want.
opt = init_MFS(N);
%opt = init_MFS(1000);
opt.Rp = Rp;
opt.fmm = fmm; 
opt.maxit = 200; %max number of gmres iterations
gmres_tol = 1e-7;
opt.plot = 0; %visualise?

%% Discretize one body
% Inner proxy surface, outer collocation grid for single sphere using spherical design nodes
[rin,rout] =  getDesignGrid(opt.Rp,opt); 

%% One-body preconditioning
%Create pseudoinverse of self-interaction matrix,
[Y,UU,LL,Kin,~] = oneBodyPrecondMob(rin,rout);

%The format is used to prepare for the case when different shapes are is
%use
Yii{1} = Y;
UUii{1} = UU; 

N = size(rin,1); %number of sources per particle
M = size(rout,1); %numer of collocation points per particle

%% Assemble grid and completion source
%Create grid on every particle. Also create completion source on every
%particle, given force and torque.

rvec_in = [];
rvec_out = [];
lambda_vec = []; 

for k = 1:P

    rvec_in = [rvec_in; rin+q(k,:)];
    rvec_out = [rvec_out; rout+q(k,:)];
    
    %Create right hand side, given forces and torques on the particles
    F = Fvec(6*(k-1)+1:6*(k-1)+3);
    T = Fvec(6*(k-1)+4:6*k);
    lambda_k = getCompletionSource(F,T,Kin);
    lambda_vec = [lambda_vec; lambda_k];
 
end

%% Get flow field due to completion source.
uvec = getFlow(lambda_vec,rvec_in,rvec_out,opt); 
uvec = -uvec;


%% Solve for source strengths
[x_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvecStokesMFS(x,rvec_in,rvec_out,q,UUii,Yii,opt,0,[],LL),uvec,3*size(rvec_out,1),opt.maxit,gmres_tol);

%check residual
abs_res = norm(matvecStokesMFS(x_gmres,rvec_in,rvec_out,q,UUii,Yii,opt,0,[],LL)-uvec);


%% Map back to the sought density in source points, determine rigid body velocities 
U = zeros(6*P,1);
for i = 1:P
    lambda_gmres((i-1)*3*N+1:i*3*N) = Y*(UU*x_gmres((i-1)*3*M+1:i*3*M));
    lambda_i = Y*(UU*x_gmres((i-1)*3*M+1:i*3*M));
    U(6*(i-1)+1:6*i) = -Kin'*lambda_i;  
end
lambda_norm = norm(lambda_gmres,inf);


%% Check residual
% Get new nodes for evaluating velocity residuals
opt.des_n = 2000; %fine grid
[~, rout_check] = getDesignGrid(opt.Rp, opt);
n_check = size(rout_check,1); 

% Assemble all check points
rcheck = [];
for k = 1:P
    x = q(k,:) + rout_check;
    rcheck = [rcheck; x];
end


%Assign velocities at checkpoints
ucheck = zeros(n_check*3*size(q,1),1); 
for k = 1:size(q,1)
    Kcheck = getKmat(rcheck(n_check*(k-1)+1:k*n_check,:),q(k,:));
    ucheck((k-1)*3*n_check+1:3*k*n_check) = Kcheck*U((k-1)*6+1:k*6);   
end

for i =1:P
    densityK_particle = (eye(3*N)-LL)*lambda_gmres(3*(i-1)*N+1:i*3*N)'+lambda_vec(3*(i-1)*N+1:i*3*N);
    densityK(3*(i-1)*N+1:i*3*N) = densityK_particle;
end

%get flow and compare RHS and LHS of representation
ubdry = getFlow(densityK,rvec_in,rcheck,opt);
uerr_vec = vecnorm(reshape(ucheck-ubdry,3,[]),2,1)/max(vecnorm(reshape(ucheck,3,[]),2,1));
uerr = max(uerr_vec);



   

end