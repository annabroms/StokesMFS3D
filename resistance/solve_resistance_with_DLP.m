function [Fvec, iters, lambda_norm, err_res] = solve_resistance_with_DLP(q, rvec_in, rvec_out, U, opt, R, E0)
%SOLVE_RESISTANCE_WITH_DLP Solve a Stokes resistance problem using an alternative
%formulation, with DS\lambda = KU
%
%   [Fvec, iters, lambda_norm, err_res] = SOLVE_RESISTANCE_WITH_DLP(q, rvec_in, rvec_out, U, opt, R, E0)
%
%   Given prescribed translational and angular velocities, the function computes the resulting hydrodynamic forces 
%   and torques acting on each particle.
%
%   INPUTS:
%       q         - P × 3 matrix of particle center positions, x,y,z
%       rvec_in   - 3NP × 1 vector of proxy source points (stacked).
%       rvec_out  - 3MP × 1 vector of collocation points on particle surfaces (stacked).
%       U         - 6P × 1 vector of prescribed rigid body velocities: [u1; omega1; u2; omega2; ...].
%       opt       - Struct containing solver options (e.g., gmres tolerance, fmm flag).
%       R         - P × 1 cell array of rotation matrices for the P particles.
%       E0        - 1 × 3 vector of semiaxes [a, b, c] of the ellipsoidal particles.
%
%   OUTPUTS:
%       Fvec        - 6P × 1 vector of computed forces and torques: [F1; T1; F2; T2; ...].
%       iters       - Number of GMRES iterations until convergence.
%       lambda_norm - Infinity norm of the final source density vector (for diagnostic use).
%       err_res     - Maximum relative residual in velocity on the particle surfaces.
%
%   METHOD OVERVIEW:
%       - Builds MFS representation from collocation and proxy surfaces.
%       - Enforces rigid body velocity boundary conditions on the true particle surfaces.
%       - Computes surface force densities using a one-body preconditioned GMRES solve.
%       - Extracts the resulting net forces and torques.
%       - Validates accuracy by evaluating velocity residuals in new check points.
%
%   DEPENDENCIES:
%       init_MFS, getDesignGrid, getVelocityData, matvecStokesMFS,
%       oneBodyPrecondRes, helsing_gmres, getKmat, getFlow
%
%   See also: SOLVE_RESISTANCE
%
%   Anna Broms, Nov 19, 2025

if nargin==0, test_solve; 
    return; end

P = size(q,1); %number of spheres

if nargin < 6
    R = eye(3);
    E0 = [1 1 1];
elseif nargin < 7
    E0 = [1 1 1];
end


N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P;

opt.N = N;
opt.M = M; 
opt.P = P; 


%% Visualize geometry
% Optional block for displaying the particle configuration 

if opt.plot
    [XX,YY,ZZ] = sphere(12);
    r = 1; % assumed unit sphere
    figure()
    for k = 1:size(q,1)
        hSurface = surf(r*XX+q(k,1), r*YY+q(k,2), r*ZZ+q(k,3));
        set(hSurface, 'FaceColor', [1 0 0], ...
                      'FaceAlpha', 0.9, ...
                      'FaceLighting', 'gouraud', ...
                      'EdgeColor', 'none');
        hold on
        axis equal
    end
    camlight
end


%% Assign RHS in resistance problem    
Kin = getKmat(rvec_in(1:N,:),[0,0,0]);

%For comparison with Kin
Kout = getKmat(rvec_out(1:M,:),[0,0,0]);
Tblock = getTraction(rvec_in(1:N,:),rvec_out(1:M,:),rvec_out(1:M,:)-q(1));
%For each particle, get data at surface, given rigid body motion
for k = 1:P
    if opt.ellipsoid
        Kin = getKmat(rvec_in(1:N,:),[0,0,0]); %If spheres, this will be the same matrix for every one.
    end
    %u_bndry((k-1)*3*N+1:3*k*N) = M/4/pi*Kin*U((k-1)*6+1:k*6);

    %should give the exact same thing, but turns out to give slightly
    %smaller residual.
    u_bndry((k-1)*3*N+1:3*k*N) = -Tblock*Kout*U((k-1)*6+1:k*6);


end

%% Compute preconditioning. It's enough to do this for a single particle 
%as it's assumed that everyone has the same shape.
if opt.ellipsoid
    [Yk,UUk,Si] = oneBodyPrecondResDLP((R{1}'*rvec_in(1:N,:)')',...
    (R{1}'*rvec_out(1:M,:)')');
else
    nvec = repmat(rvec_out(1:M,:),P,1);
    [Yk,UUk,Si] = oneBodyPrecondResDLP(rvec_in(1:N,:),rvec_out(1:M,:),nvec);
end

%The format is used to prepare for the case when different shapes are is
%use
Y{1} = Yk;
UU{1} = UUk; 
disp('Done preconditioning')

if opt.profile
    memorygraph('label','start matvec resistance');
end

%% What do eigvals look like for the system matrix?

debug = 0; 
if debug
    x = zeros(3*N*P,1);
    tic
    for k = 1:3*N*P
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvecStokesMFS_DLP(x,rvec_in,rvec_out,q,UU,Y,Tblock,opt,1,R);
        CC(:,k) = uu;
    end
    toc
    figure(13);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    skeel(CC)

    figure(5)
    [V,D] = eig(CC);
    D = diag(D); 
    plot(real(D),imag(D),'ro')
end

%% Prepare long-range preconditioning with deflation

u_bndry = u_bndry';
if opt.lr
    error('not yet implemented with DLP extension')
    if opt.lr > 2 %use svd for single body
        [Rinv,Zi,Yi,db] = get_long_range_precond(q,rvec_in,rvec_out,R,opt,Yk*diag(Si),UUk');
    else
        [Rinv,Zi,Yi,db] = get_long_range_precond(q,rvec_in,rvec_out,R,opt);
    end
    opt.db = db; 
    tau_coarse = getCoarseSource(u_bndry,Rinv,Zi,Yi,R,opt);
    disp('Done coarse projection')
else
    disp('No deflation')
end


%% Solve problem
verbose = 1; 
if opt.lr
    %not yet implemented
    Pf = applyPmat(u_bndry,rvec_in,rvec_out,Rinv,Zi,Yi,R,opt); 
    [mu_gmres,iters,resvec,real_res] = helsing_gmres(@(x) lr_matvecStokesMFS(x,rvec_in,rvec_out,q,UU,Y,opt,R,Rinv,Zi,Yi),Pf,3*size(rvec_out,1),opt.maxit,opt.gmres_tol,verbose,0);
else
    [mu_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvecStokesMFS_DLP(x,rvec_in,rvec_out,q,UU,Y,Tblock,opt,1,R),u_bndry,3*size(rvec_in,1),opt.maxit,opt.gmres_tol,verbose,0);
end

if opt.profile
    memorygraph('label','done matvec resistance, remap and determine force')
end

%% Determine source strengths on proxy sources from the solution at the boundary,
% lambda <- mu. Then, determine forces and torques on particles, given lambda

Fvec = zeros(6*P,1);
Kin = getKmat(rvec_in(1:N,:),[0,0,0]);
for i = 1:P
    if opt.ellipsoid
        temp_i = Y{1}*(UU{1}*(rotate_vector(mu_gmres((i-1)*3*N+1:i*3*N),R{i}')));
        lambda_i = rotate_vector(temp_i,R{i});
        if ~opt.lr
            Kin = getKmat(rvec_in(N*(i-1)+1:N*i,:),q(i,:));
        end
    else
        lambda_i = Y{1}*(UU{1}*(mu_gmres((i-1)*3*N+1:i*3*N)));
    end
    lambda_gmres(3*(i-1)*N+1:i*3*N) = lambda_i;

    if ~opt.lr
        Fvec(6*(i-1)+1:6*i) = Kin'*lambda_i; 
    end
end

%lambda_gmres = M*Y{1}*(UU{1}*Kin*U)/4/pi; %should be the same thing! 

if opt.lr
    %not yet implemented
    lambda_gmres_old = lambda_gmres';
    tau_stokes = applyQmat(lambda_gmres_old,rvec_in,rvec_out,Rinv,Zi,Yi,R,opt);
    lambda_gmres = tau_stokes+tau_coarse; 

    for i = 1:P
        if opt.ellipsoid
            Kin = getKmat(rvec_in(N*(i-1)+1:N*i,:),q(i,:));
        end

        lambda_i = lambda_gmres(3*(i-1)*N+1:i*3*N);
        Fvec(6*(i-1)+1:6*i) = Kin'*lambda_i; 

    end

end


% lambda_norm gives a sanity check on the source distribution. If large,
% the representation for the problem is not optimal.
lambda_norm = norm(lambda_gmres, inf);

%% Estimate surface residual
% Evaluate how well the flow generated by the sources matches the prescribed boundary velocity (from bc).
if opt.profile
    memorygraph('label','done solving resistance, compute velocities')
end

% Get new nodes for evaluating velocity residuals
b = ellipsoid_param(E0(1),E0(2),E0(3));   % baseline object at the oridin, aligned
b = setupsurfquad(b,[46,55]);

rcheck = []; 
for k = 1:P
    if opt.ellipsoid
        x = q(k,:)+(R{k}*b.x)';
    else
        x = q(k,:) + b.x'; 
    end

    rcheck = [rcheck; x];
 
end
n_check = size(b.x,2);

% Assign expected velocity at check points using prescribed U
ucheck = zeros(n_check*3*P,1);
for k = 1:P
    Kcheck = getKmat(rcheck(n_check*(k-1)+1:k*n_check,:), q(k,:));
    ucheck((k-1)*3*n_check+1:3*k*n_check) = Kcheck * U((k-1)*6+1:k*6);
end

% Evaluate flow from solution to resistance problem
ubdry = getFlow(lambda_gmres, rvec_in, rcheck, opt);

% Compute relative residual
uerr_vec = vecnorm(reshape(ucheck - ubdry, 3, []), 2, 1) ./ ...
           max(vecnorm(reshape(ucheck, 3, []), 2, 1));
err_res = max(uerr_vec)



end

function test_solve
rng(5); %reproducable

P = 10; %number of bodies
delta = 7; %smallest particle particle distance 
q = [0 0 0; 2+delta 0 0]; %center coordiante matrix for P particles, x,y,z: size P x 3
%q = [0 0 0]; 

%random configurations
[q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
  
fmm = 0; %only activate if many particles (say, more than 40)

%% Solve resistance problem (given velocities)
disp('Start with resistance: ')
Uref = rand(6*P,1); 

%Test first with very low resolution
Rp = 0.30;
N = 100; 

Rp = 0.15; %Proxy sphere radius -- very coarse resolution
N = 50;  % Number of proxy sources
a = 2; 
%a = 2; %or play with SVD truncation level

% Rp = 0.68; %proxy radius
% N = 700; % approximate number of proxy sources on every particle


[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a);
%[rvec_in,rvec_out,opt] = init_spheres(q);
opt.fmm = fmm;
opt.gmres_tol = 1e-3;
opt.lr = 0;
%opt.gmres_tol = 1e-10; %does this matter for the accuracy? 
[Fvec,it_res,lambda_norm_res,err_res] = solve_resistance_with_DLP(q,rvec_in,rvec_out,Uref, opt);
opt.gmres_tol = 1e-10; 
[Fvec2,it_res2,lambda_norm_res2,err_res2] = solve_resistance(q,rvec_in,rvec_out,Uref, opt); 

norm(Fvec-Fvec2,inf)/norm(Fvec2,inf)

end

