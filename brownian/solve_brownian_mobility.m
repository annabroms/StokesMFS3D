function [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec,randvel, opt,sqrtA,dW,R,E0)
%SOLVE_BROWNIAN_MOBILITY Solve a Stokes mobility problem with fluctuating velocity field for a configuration of ellipsoidal particles using MFS.
%
%   [U, iters, lambda_norm, uerr] = SOLVE_BROWNIAN_MOBILITY(q, rvec_in, rvec_out, Fvec, opt, sqrtA,W)
%
%   Given the forces and torques applied to each particle, the function computes the resulting translational and rotational
%   velocities.
%
%   INPUTS:
%       q         - P × 3 matrix of particle center positions.
%       rvec_in   - 3NP × 1 vector of collocation points on particle surfaces (stacked).
%       rvec_out  - 3MP × 1 vector of proxy source points (stacked).
%       Y         - Factor in one-body pseudoinverse
%       UU        - Factor in one-body pseudoinverse, such that Y*U' is
%                   TS_L^{(ii}).
%       LL        - Projection matrix for sources not to contribute to net
%                   force or torque.
%       Tblock    - DLP matrix for single body
%       Fvec      - 6P × 1 vector of applied forces and torques, format: [F1; T1; F2; T2; ...].
%       opt       - Struct containing solver options (e.g., gmres tolerance, fmm flag).
%       sqrtA     - Matrix of size 3NP x 3NP needed to set the correlation of the fluctuating
%                   velocity field
%       dW         - Vector of size 3NP x 1 containing scaled Wiener increment (drawn from standard normal
%                   distribution)
%       R         - P x 1 cell array of rotation matrices for the P
%                   ellipsoids [NOT YET SUPPORTED]
%       E0        - 1 × 3 vector of semiaxes [a, b, c] of the ellipsoidal particles.
%
%   OUTPUTS:
%       U           - 6P × 1 vector of resulting rigid body velocities: [u1; omega1; u2; omega2; ...].
%       iters       - Number of GMRES iterations until convergence.
%       lambda_norm - Infinity norm of the final density vector (for diagnostic use).
%       uerr        - Maximum relative residual of the velocity field on the surface.
%
%   METHOD OVERVIEW:
%       - Builds MFS representation from collocation and proxy surfaces.
%       - Uses a "completion source" to represent force and torque.
%       - Computes rigid body motion from the computed source density.
%       - Determines the surface residuals in new points.
%
%   DEPENDENCIES:
%       init_MFS, getDesignGrid, getCompletionSource, matvecStokesMFS, 
%       oneBodyPrecondMob, helsing_gmres, getKmat
%
%   See also: LARGE_ELLIPSOID_EX, SOLVE_RESISTANCE
%
%   Anna Broms, Aug 4, 2025

P = size(q,1);
N = size(rvec_in,1)/P; %number of sources per particle
M = size(rvec_out,1)/P; %numer of collocation points per particle

opt.M = M; 
opt.N = N;

getSqrt = 0; 


if nargin < 12
    getSqrt = 1;
    R = eye(3);
    E0 = [1 1 1];
elseif nargin<13
    dW =  randn(3*N*P,1);
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 14
    getSqrt = 0;
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 15
    getSqrt = 0; 
    E0 = [1 1 1];
end



%% One-body preconditioning already computed 
%The format is used to prepare for the case when different shapes are is
%use
Yii{1} = Y;
UUii{1} = UU; 

%% Assemble completion source, given force and torque

lambda_vec = []; 

for k = 1:P

    %Create right hand side, given forces and torques on the particles
    F = Fvec(6*(k-1)+1:6*(k-1)+3);
    T = Fvec(6*(k-1)+4:6*k);
    
    if opt.ellipsoid
        Rk = R{k};
        if any([F;T])
            lambda_k = getCompletionSource(Rk'*F,Rk'*T,Kin);
        else
            lambda_k = zeros(3*N,1); 
        end
        lambda_vec = [lambda_vec; rotate_vector(lambda_k,Rk)];
    else
        if any([F;T])
            lambda_k = getCompletionSource(F,T,Kin);
        else
            lambda_k = zeros(3*N,1); 
        end
        lambda_vec = [lambda_vec; lambda_k];
    end

    
 
end

%% Get flow field due to completion source.
uvec = getFlow(lambda_vec,rvec_in,rvec_out,opt); 

%multiply through with blockdiag T
uvec_T = zeros(3*N*P,1);
for k = 1:P
    uvec_T((k-1)*3*N+1:k*3*N) = Tblock*uvec((k-1)*3*M+1:k*3*M);
end

%% Get flow due to fluctuating velocity field
if randvel
    if getSqrt
        %urand = getRandomField(q,rvec_in,rvec_out,P);
    else
        urand = sqrtA*dW;
    end
    uvec_T = uvec_T - sqrt(2/opt.dt)*urand;
   
end



%% Solve for source strengths
[x_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvecStokesMFS_DLP(x,rvec_in,rvec_out,q,UUii,Yii,Tblock,opt,0,R,LL),uvec_T,3*size(rvec_in,1),opt.maxit,opt.gmres_tol);

%check residual
%abs_res = norm(matvecStokesMFS(x_gmres,rvec_in,rvec_out,q,UUii,Yii,opt,0,R,LL)-uvec);


%% Map back to the sought density in source points, determine rigid body velocities 
U = zeros(6*P,1);
lambda_gmres = zeros(N*P*3,1);
for i = 1:P
    if opt.ellipsoid
        temp_i = Y*(UU*(rotate_vector(x_gmres((i-1)*3*N+1:i*3*N),R{i}')));
        lambda_i = rotate_vector(temp_i,R{i});
        Kin_i = getKmat(rvec_in(N*(i-1)+1:N*i,:),q(i,:));
        U(6*(i-1)+1:6*i) = -Kin_i'*lambda_i; 
    else
        lambda_i = Y*(UU*x_gmres((i-1)*3*N+1:i*3*N));
        U(6*(i-1)+1:6*i) = -Kin'*lambda_i; 
    end
    
    lambda_gmres((i-1)*3*N+1:i*3*N) = lambda_i;    
     
end
lambda_norm = norm(lambda_gmres,inf);


   

end