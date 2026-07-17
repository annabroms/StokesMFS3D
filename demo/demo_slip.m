clear;
close all
%%%%%% Demo for Stokes MFS in 3D with prescribed slip on P unit spheres.
%
% The particles are force- and torque-free. Their tangential slip uses the
% B1 and C1 modes from Pak and Lauga, J. Eng. Math. 88 (2014). With the
% convention in this repository, the boundary condition on particle k is
%
%   u_slip = -B1*(p-(n dot p)*n) - C1*(c x n),
%   u = V_k + Omega_k x (x-q_k) + u_slip.
%
% Thus u_slip is the physical surface-distortion velocity u_s used by Pak
% and Lauga.
%
% The surface residual is evaluated on a different spherical design from
% the one used by the solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Particle configuration
P = 4;                    % number of particles
delta = 1;                % surface-to-surface separation
radius = 1;
q = [(0:P-1)'*(2+delta), zeros(P,2)];

%% Pak--Lauga slip parameters
B1 = ones(P,1);           % translational squirming strength
C1 = 0.25*ones(P,1);      % rotational squirming strength
p_dir = repmat([0 0 1],P,1); % translational axes
c_dir = repmat([0 1 0],P,1); % rotational axes

%% MFS discretization
Rp = 0.7;                 % proxy-sphere radius
N_request = 400;          % approximate proxy points per particle
a_glob = 1.2;             % collocation oversampling
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N_request,a_glob);

opt.fmm = 0;              % use direct summation for this small example
opt.gmres_tol = 1e-10;
opt.gmres_verbose = 0;
opt.compute_residual = false; % slip is not known on the solver check grid
opt.plot = 0;

N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P;

%% Prescribed slip on the collocation points
normals = rvec_out - repelem(q,M,1);
normals = normals ./ vecnorm(normals,2,2);
u_slip = zeros(P*M,3);

for k = 1:P
    rows = (k-1)*M + (1:M);
    u_slip(rows,:) = pak_lauga_slip( ...
        normals(rows,:),B1(k),C1(k),p_dir(k,:),c_dir(k,:));
end

Fvec = zeros(6*P,1);      % force- and torque-free particles
[U,iters,lambda_norm,~,lambda] = solve_mobility( ...
    q,rvec_in,rvec_out,Fvec,u_slip,opt);

%% Surface residual on a new set of points
N_check_request = 1000;
[normal_check,w_check_one] = get_sphdesign(N_check_request);
N_check = size(normal_check,1);

rcheck = zeros(P*N_check,3);
u_slip_check = zeros(P*N_check,3);
u_rigid_check = zeros(3*P*N_check,1);

for k = 1:P
    point_rows = (k-1)*N_check + (1:N_check);
    velocity_rows = (k-1)*3*N_check + (1:3*N_check);
    body_velocity = U((k-1)*6 + (1:6));

    rcheck(point_rows,:) = q(k,:) + radius*normal_check;
    u_slip_check(point_rows,:) = pak_lauga_slip( ...
        normal_check,B1(k),C1(k),p_dir(k,:),c_dir(k,:));
    u_rigid_check(velocity_rows) = ...
        getKmat(rcheck(point_rows,:),q(k,:))*body_velocity;
end

u_mfs_check = getStokesletFlow(lambda,rvec_in,rcheck,opt);
u_slip_check_vec = reshape(u_slip_check.',[],1);
boundary_residual = ...
    u_mfs_check - u_rigid_check - u_slip_check_vec;

w_check = repmat(radius^2*w_check_one,P,1);
residual_at_nodes = reshape(boundary_residual,3,[]).';
residual_norm = sqrt(sum(w_check.*sum(residual_at_nodes.^2,2)));
slip_norm = sqrt(sum(w_check.*sum(u_slip_check.^2,2)));
relative_surface_residual = residual_norm/slip_norm;

%% Report
fprintf('\nPrescribed-slip mobility demo\n');
fprintf('  P = %d, N = %d, M = %d, check points/body = %d\n', ...
    P,N,M,N_check);
fprintf('  GMRES iterations = %d, source norm = %.3e\n', ...
    iters,lambda_norm);
fprintf('  Relative surface residual = %.3e\n\n', ...
    relative_surface_residual);

fprintf(' particle       Vx       Vy       Vz  Omega_x  Omega_y  Omega_z\n');
for k = 1:P
    body_velocity = U((k-1)*6 + (1:6));
    fprintf('%9d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
        k,body_velocity);
end

function u_slip = pak_lauga_slip(normals,B1,C1,p_dir,c_dir)
% Evaluate the tangential B1 and C1 slip modes on a unit sphere.
p_values = repmat(p_dir,size(normals,1),1);
c_values = repmat(c_dir,size(normals,1),1);
u_slip = -B1*(p_values-sum(normals.*p_values,2).*normals) - ...
    C1*cross(c_values,normals,2);
end
