function [Fvec, iters, lambda_norm, err_res] = solve_resistance(q, rvec_in, rvec_out, U, opt, R, E0)
%SOLVE_RESISTANCE Solve a Stokes resistance problem for a configuration of ellipsoidal particles using MFS.
%
%   [Fvec, iters, lambda_norm, err_res] = SOLVE_RESISTANCE(q, rvec_in, rvec_out, U, opt, R, E0)
%
%   Given prescribed translational and angular velocities, the function computes the resulting hydrodynamic forces 
%   and torques acting on each particle.
%
%   INPUTS:
%       q         - P × 3 matrix of particle center positions, x,y,z
%       rvec_in   - (N*P) × 3 array of proxy source points (stacked by particle).
%       rvec_out  - (M*P) × 3 array of collocation points on particle surfaces (stacked by particle).
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
%       oneBodyPrecondRes, helsing_gmres, getKmat, getStokesletFlow,
%       ellipsoid_param, setupsurfquad, getTractionFast
%
%   See also: LARGE_ELLIPSOID_EX, SOLVE_MOBILITY
%
%   NOTES:
%       - If opt.plot is true, a visualization of the geometry is produced.
%       - Boundary velocity and traction visualization is computed only when opt.plot is true.
%
%   Anna Broms, June 13, 2025

if nargin==0, test_solve; 
    return; end


P = size(q,1); %number of spheres

if ~isfield(opt,'gmres_verbose')
    opt.gmres_verbose = 1;
end
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

if opt.ellipsoid
    rin_block = (R{1}' * (rvec_in(1:N,:) - q(1,:))')';
    rout_block = (R{1}' * (rvec_out(1:M,:) - q(1,:))')';
    opt.Sblock = stokes_SLP_mat(rin_block, rout_block);
else
    opt.Sblock = stokes_SLP_mat(rvec_in(1:N,:), rvec_out(1:M,:));
end


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
    title('Particle configuration')
end


%% Assign RHS in resistance problem    

Kout = getKmat(rvec_out(1:M,:),[0,0,0]);
%For each particle, get data at surface, given rigid body motion
for k = 1:P
    if opt.ellipsoid
        Kout = getKmat(rvec_out(M*(k-1)+1:k*M,:),q(k,:)); %If spheres, this will be the same matrix for every one.
    end
    u_bndry((k-1)*3*M+1:3*k*M) = Kout*U((k-1)*6+1:k*6);
end

%% Compute preconditioning. It's enough to do this for a single particle 
%as it's assumed that everyone has the same shape.
if opt.ellipsoid
    [Yk,UUk,Si] = oneBodyPrecondRes((R{1}'*rvec_in(1:N,:)')',...
    (R{1}'*rvec_out(1:M,:)')');
else
    [Yk,UUk,Si] = oneBodyPrecondRes(rvec_in(1:N,:),rvec_out(1:M,:));
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
    x = zeros(3*M*P,1);
    tic
    for k = 1:3*M*P
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvecStokesMFS(x,rvec_in,rvec_out,q,UU,Y,opt,1,R);
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
    if opt.lr > 2 %use svd for single body
        [Rinv,Zi,Yi,db] = get_long_range_precond(q,rvec_in,rvec_out,R,opt,Yk*diag(Si),UUk');
    else
        [Rinv,Zi,Yi,db] = get_long_range_precond(q,rvec_in,rvec_out,R,opt);
    end
    opt.db = db; 
    tau_coarse = getCoarseSource(u_bndry,Rinv,Zi,Yi,R,opt);
    disp('Done coarse projection')
else
    disp('No deflation preconditioning')
end


%% Solve problem
verbose = opt.gmres_verbose; 
if opt.lr
    Pf = applyPmat(u_bndry,rvec_in,rvec_out,Rinv,Zi,Yi,R,opt); 
    [mu_gmres,iters,resvec,real_res] = helsing_gmres(@(x) lr_matvecStokesMFS(x,rvec_in,rvec_out,q,UU,Y,opt,R,Rinv,Zi,Yi),Pf,3*size(rvec_out,1),opt.maxit,opt.gmres_tol,verbose,0);
else
    [mu_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvecStokesMFS(x,rvec_in,rvec_out,q,UU,Y,opt,1,R),u_bndry,3*size(rvec_out,1),opt.maxit,opt.gmres_tol,verbose,0);
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
        temp_i = Y{1}*(UU{1}*(rotate_vector(mu_gmres((i-1)*3*M+1:i*3*M),R{i}')));
        lambda_i = rotate_vector(temp_i,R{i});
        if ~opt.lr
            Kin = getKmat(rvec_in(N*(i-1)+1:N*i,:),q(i,:));
        end
    else
        lambda_i = Y{1}*(UU{1}*(mu_gmres((i-1)*3*M+1:i*3*M)));
    end
    lambda_gmres(3*(i-1)*N+1:i*3*N) = lambda_i;

    if ~opt.lr
        Fvec(6*(i-1)+1:6*i) = Kin'*lambda_i; 
    end
end

if opt.lr
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
ubdry_check = getStokesletFlow(lambda_gmres, rvec_in, rcheck, opt);

% Compute relative residual
uerr_vec = vecnorm(reshape(ucheck - ubdry_check, 3, []), 2, 1) ./ ...
           max(vecnorm(reshape(ucheck, 3, []), 2, 1));
err_res = max(uerr_vec);

%% Optional visualization: color by surface velocity and traction magnitudes
if opt.plot
    ubdry_surf = getStokesletFlow(lambda_gmres, rvec_in, rvec_out, opt);
    umag = vecnorm(reshape(ubdry_surf,3,[]),2,1).';
    nvec = rvec_out - kron(q, ones(M,1));
    nvec = nvec ./ vecnorm(nvec,2,2);
    traction = getTractionFast(lambda_gmres, rvec_in, rvec_out, nvec, opt);
    tmag = vecnorm(reshape(traction,3,[]),2,1).';

    plot_surface_scalar(rvec_out, M, P, umag, ...
        'Resistance: surface velocity magnitude', 'parula');
    plot_surface_scalar(rvec_out, M, P, tmag, ...
        'Resistance: traction magnitude on surface', 'hot');
end



end

function test_solve
% Self-test: mobility -> resistance -> compare forces/torques.
rng(5); %reproducable
close all; 

P = 8; %number of bodies
delta = 2; %smallest particle particle distance
[q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta

fmm = 0; %only activate if many particles (say, more than 40)

Fref = rand(6*P,1);

[rvec_in,rvec_out,opt] = init_spheres(q);
opt.fmm = fmm;
opt.lr = 0;
opt.gmres_tol = 1e-10;
opt.plot = 1; 

[U,it_mob,lambda_norm_mob,err_mob]  = solve_mobility( ...
    q,rvec_in,rvec_out,Fref,[],opt);
[Fvec,it_res,lambda_norm_res,err_res] = solve_resistance(q,rvec_in,rvec_out,U, opt);

rel_err = norm(Fvec-Fref,inf)/norm(Fref,inf);

fprintf('\nSelf-test summary (solve_resistance)\n');
fprintf('  P = %d, delta = %.3g, N = %d, M = %d\n', P, delta, size(rvec_in,1)/P, size(rvec_out,1)/P);
fprintf('  gmres_tol = %.1e\n', opt.gmres_tol);
fprintf('  Mobility:   iters = %d, lambda_norm = %.3e, rel surf residual = %.3e\n', it_mob, lambda_norm_mob, err_mob);
fprintf('  Resistance: iters = %d, lambda_norm = %.3e, rel surf residual = %.3e\n', it_res, lambda_norm_res, err_res);
fprintf('  Relative force/torque error = %.3e\n', rel_err);

alignfigs;

end
