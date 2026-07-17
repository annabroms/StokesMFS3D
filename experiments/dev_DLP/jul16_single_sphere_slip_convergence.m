%JUL16_SINGLE_SPHERE_SLIP_CONVERGENCE
%
% Validate nonzero prescribed slip for a force- and torque-free unit sphere
% using solve_mobility and solve_mobility_with_DLP independently. The exact
% swimming velocities follow from the reciprocal-theorem identities in
% Section 4.2 of:
%
%   O. S. Pak and E. Lauga, "Generalized squirming motion of a sphere,"
%   Journal of Engineering Mathematics 88 (2014), 1--28.
%   https://www.damtp.cam.ac.uk/user/lauga/papers/90.pdf
%
% Pak and Lauga denote the physical surface-distortion velocity by u_s and
% obtain velocities proportional to -u_s. The boundary convention here is
%
%   u = V + Omega x (x-c) + u_slip,
%
% so u_slip = u_s. The slip modes below include the signs needed for the
% exact velocities to be positive along p_dir and q_dir.
%
% Full experiment:
%   run('experiments/dev_DLP/jul16_single_sphere_slip_convergence.m')
%
% Created July 16, 2026.

clear;
close all;
clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
repo_root = find_repo_root(script_dir);
run(fullfile(repo_root,'setpath.m'));

%% Physical configuration
centre = [0 0 0];
radius = 1;
B1 = 1;
C1 = 1;
p_dir = [1; 2; 3] / sqrt(14);
q_dir = [-2; 1; 1] / sqrt(6);

V_exact = (2*B1/3) * p_dir;
Omega_exact = (C1/radius) * q_dir;
Fvec = zeros(6,1);

%% Discretization and solver configuration
N_request_vec = [100 200 400 700 1000 1500];
Rp = 0.7;
a_glob = 1.2;
gmres_tol = 1e-12;
maxit = 600;
use_fmm = false;
fmm_tol = 1e-12;

dlp_inner_only = false;
dlp_symmetrize_weighted = false;
dlp_add_rank1 = false;
dlp_outer_force = true;

%% Independent rotated check grid
N_check_request = 5000;
[n_check_base,~] = get_sphdesign(N_check_request);
rotation_axis = [2; -3; 1] / sqrt(14);
rotation_angle = pi/7;
R_check = axis_angle_rotation(rotation_axis,rotation_angle);

n_check = (R_check * n_check_base.').';
r_check = centre + radius*n_check;
u_slip_check = evaluate_slip(n_check,B1,C1,p_dir,q_dir);
u_slip_check_vec = reshape(u_slip_check.',[],1);
slip_check_norm = norm(u_slip_check_vec,inf);

%% Storage
n_cases = numel(N_request_vec);
results = struct();
results.created = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));
results.description = [ ...
    'Single force- and torque-free sphere with prescribed translational ', ...
    'and rotational squirming modes.'];
results.reference = [ ...
    'Pak and Lauga, Generalized squirming motion of a sphere, ', ...
    'J. Eng. Math. 88 (2014), Section 4.2.'];
results.config.centre = centre;
results.config.radius = radius;
results.config.B1 = B1;
results.config.C1 = C1;
results.config.p_dir = p_dir;
results.config.q_dir = q_dir;
results.config.N_request = N_request_vec;
results.config.N_check_request = N_check_request;
results.config.N_check_actual = size(r_check,1);
results.config.check_rotation_axis = rotation_axis;
results.config.check_rotation_angle = rotation_angle;
results.config.Rp = Rp;
results.config.a_glob = a_glob;
results.config.gmres_tol = gmres_tol;
results.config.maxit = maxit;
results.config.use_fmm = use_fmm;
results.config.fmm_tol = fmm_tol;
results.config.inner_only = dlp_inner_only;
results.config.symmetrize_weighted = dlp_symmetrize_weighted;
results.config.add_rank1 = dlp_add_rank1;
results.config.outer_force = dlp_outer_force;
results.exact.V = V_exact;
results.exact.Omega = Omega_exact;
results.actual.N = nan(1,n_cases);
results.actual.M = nan(1,n_cases);
results.tangential_defect = nan(1,n_cases);
method_results.V = nan(3,n_cases);
method_results.Omega = nan(3,n_cases);
method_results.error.trans = nan(1,n_cases);
method_results.error.rot = nan(1,n_cases);
method_results.error.surface = nan(1,n_cases);
method_results.iters = nan(1,n_cases);
method_results.lambda_total_norm = nan(1,n_cases);
method_results.solve_time = nan(1,n_cases);
results.standard = method_results;
results.dlp = method_results;
N_visualization = -inf;
u_total_visualization.standard = [];
u_total_visualization.dlp = [];

fprintf('Single-sphere prescribed-slip convergence\n');
fprintf('  B1 = %.3g, C1 = %.3g, radius = %.3g, Rp = %.3g\n', ...
    B1,C1,radius,Rp);
fprintf('  requested N = [%s]\n',num2str(N_request_vec));
fprintf('  independent rotated check nodes = %d\n',size(r_check,1));
fprintf('  gmres_tol = %.1e, maxit = %d\n\n',gmres_tol,maxit);
fprintf('%9s %6s %6s %6s %7s %11s %11s %11s %11s %9s\n', ...
    'method','N','M','iters','time','e_trans','e_rot', ...
    'e_surface','|lambda|inf','n dot us');

%% Convergence sweep
for i = 1:n_cases
    N_request = N_request_vec(i);
    [rvec_in,rvec_out,opt] = init_spheres( ...
        centre,Rp,N_request,a_glob);

    N_actual = size(rvec_in,1);
    M_actual = size(rvec_out,1);
    [nout,wout] = sphere_normals_and_weights( ...
        centre,rvec_out,radius);
    u_slip = evaluate_slip(nout,B1,C1,p_dir,q_dir);

    tangential_defect = max(abs(sum(nout.*u_slip,2)));
    if tangential_defect > 1e-12
        error('jul16_single_sphere_slip:tangentiality', ...
            'Slip tangentiality defect is %.3e.',tangential_defect);
    end

    opt = configure_options( ...
        opt,gmres_tol,maxit,use_fmm,fmm_tol, ...
        dlp_inner_only,dlp_symmetrize_weighted, ...
        dlp_add_rank1,dlp_outer_force);

    standard_start = tic;
    [U_standard,iters_standard,~,~,lambda_standard] = solve_mobility( ...
        centre,rvec_in,rvec_out,Fvec,u_slip,opt);
    solve_time_standard = toc(standard_start);

    dlp_start = tic;
    [U_dlp,iters_dlp,~,~,lambda_dlp] = solve_mobility_with_DLP( ...
        centre,rvec_in,rvec_out,nout,wout,Fvec,u_slip,opt);
    solve_time_dlp = toc(dlp_start);

    [V_standard,Omega_standard,e_trans_standard,e_rot_standard, ...
        e_surface_standard,u_mfs_standard] = evaluate_solution( ...
        U_standard,lambda_standard,rvec_in,r_check,centre,opt, ...
        V_exact,Omega_exact,u_slip_check_vec,slip_check_norm);
    [V_dlp,Omega_dlp,e_trans_dlp,e_rot_dlp, ...
        e_surface_dlp,u_mfs_dlp] = evaluate_solution( ...
        U_dlp,lambda_dlp,rvec_in,r_check,centre,opt, ...
        V_exact,Omega_exact,u_slip_check_vec,slip_check_norm);

    results.actual.N(i) = N_actual;
    results.actual.M(i) = M_actual;
    results.tangential_defect(i) = tangential_defect;
    results.standard.V(:,i) = V_standard;
    results.standard.Omega(:,i) = Omega_standard;
    results.standard.error.trans(i) = e_trans_standard;
    results.standard.error.rot(i) = e_rot_standard;
    results.standard.error.surface(i) = e_surface_standard;
    results.standard.iters(i) = iters_standard;
    results.standard.lambda_total_norm(i) = norm(lambda_standard,inf);
    results.standard.solve_time(i) = solve_time_standard;
    results.dlp.V(:,i) = V_dlp;
    results.dlp.Omega(:,i) = Omega_dlp;
    results.dlp.error.trans(i) = e_trans_dlp;
    results.dlp.error.rot(i) = e_rot_dlp;
    results.dlp.error.surface(i) = e_surface_dlp;
    results.dlp.iters(i) = iters_dlp;
    results.dlp.lambda_total_norm(i) = norm(lambda_dlp,inf);
    results.dlp.solve_time(i) = solve_time_dlp;

    if N_actual > N_visualization
        N_visualization = N_actual;
        u_total_visualization.standard = u_mfs_standard;
        u_total_visualization.dlp = u_mfs_dlp;
    end

    fprintf('%9s %6d %6d %6d %7.2f %11.3e %11.3e %11.3e %11.3e %9.2e\n', ...
        'standard',N_actual,M_actual,iters_standard,solve_time_standard, ...
        e_trans_standard,e_rot_standard,e_surface_standard, ...
        results.standard.lambda_total_norm(i),tangential_defect);
    fprintf('%9s %6d %6d %6d %7.2f %11.3e %11.3e %11.3e %11.3e %9.2e\n', ...
        'DLP',N_actual,M_actual,iters_dlp,solve_time_dlp, ...
        e_trans_dlp,e_rot_dlp,e_surface_dlp, ...
        results.dlp.lambda_total_norm(i),tangential_defect);
end
results.visualization.N = N_visualization;

validate_convergence(results.standard,'standard',maxit);
validate_convergence(results.dlp,'DLP',maxit);

%% Convergence plot
figure('Color','w','Name','Single-sphere prescribed-slip convergence');
layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
metric_names = {'trans','rot','surface'};
metric_titles = { ...
    '$e_{\mathrm{trans}}$', ...
    '$e_{\mathrm{rot}}$', ...
    '$e_{\partial\Omega}$'};
for i_metric = 1:numel(metric_names)
    ax = nexttile(layout);
    metric_name = metric_names{i_metric};
    semilogy(ax,results.actual.N, ...
        results.standard.error.(metric_name),'o-', ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName','standard');
    hold(ax,'on');
    semilogy(ax,results.actual.N, ...
        results.dlp.error.(metric_name),'s--', ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName','DLP');
    axis(ax,'tight');
    grid(ax,'on');
    box(ax,'on');
    set(ax,'TickLabelInterpreter','latex','FontSize',12);
    xlabel(ax,'$N$','Interpreter','latex');
    ylabel(ax,'$\mathrm{relative\ error}$','Interpreter','latex');
    title(ax,metric_titles{i_metric},'Interpreter','latex');
    legend(ax,'Location','best','Interpreter','latex');
end
title(layout, ...
    '$\mathrm{Single\ spherical\ swimmer\ with\ prescribed\ slip}$', ...
    'Interpreter','latex');

%% Prescribed slip velocity on the particle surface
arrow_stride = 10;
plot_surface_vector_field( ...
    r_check,u_slip_check_vec,arrow_stride, ...
    'Prescribed slip velocity', ...
    '$\mathrm{Prescribed\ slip\ velocity}\ \widehat{\mathbf{u}}_{\mathrm{slip}}$', ...
    '$||u_{slip}||_2$');

%% Total velocities on the particle surface at the finest resolution
plot_surface_vector_field( ...
    r_check,u_total_visualization.standard,arrow_stride, ...
    'Standard total surface velocity', ...
    sprintf('$\\mathrm{Standard\\ total\\ surface\\ velocity},\\ N=%d$', ...
    N_visualization), ...
    '$||u_N||_2$');
plot_surface_vector_field( ...
    r_check,u_total_visualization.dlp,arrow_stride, ...
    'DLP total surface velocity', ...
    sprintf('$\\mathrm{DLP\\ total\\ surface\\ velocity},\\ N=%d$', ...
    N_visualization), ...
    '$||u_N||_2$');

%% Local functions

function repo_root = find_repo_root(start_dir)
% Locate the repository root so the experiment is independent of pwd.
repo_root = start_dir;
while true
    if isfile(fullfile(repo_root,'setpath.m'))
        return;
    end

    parent_dir = fileparts(repo_root);
    if isempty(parent_dir) || strcmp(parent_dir,repo_root)
        error('jul16_single_sphere_slip:repoRootNotFound', ...
            'Could not find setpath.m starting from %s.',start_dir);
    end
    repo_root = parent_dir;
end
end

function opt = configure_options( ...
        opt,gmres_tol,maxit,use_fmm,fmm_tol, ...
        inner_only,symmetrize_weighted,add_rank1,outer_force)
% Set the deterministic DLP solver options used throughout the sweep.
opt.fmm = use_fmm;
opt.fmm_tol = fmm_tol;
opt.lr = 0;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = 0;
opt.maxit = maxit;
opt.plot = 0;
opt.debug = 0;
opt.profile = 0;
opt.compute_residual = false;
opt.inner_only = inner_only;
opt.symmetrize_weighted = symmetrize_weighted;
opt.add_rank1 = add_rank1;
opt.outer_force = outer_force;
opt.transform_slip = true;
end

function [nout,wout] = sphere_normals_and_weights(centre,rvec_out,radius)
% Construct outward unit normals and equal surface quadrature weights.
M = size(rvec_out,1);
nout = (rvec_out-centre) / radius;
nout = nout ./ vecnorm(nout,2,2);
wout = (4*pi*radius^2/M) * ones(M,1);
end

function u_slip = evaluate_slip(n,B1,C1,p_dir,q_dir)
% Evaluate the tangential translational and rotational slip modes.
n_nodes = size(n,1);
p_values = repmat(p_dir.',n_nodes,1);
q_values = repmat(q_dir.',n_nodes,1);
u_slip = -B1 * (p_values - (n*p_dir).*n) - ...
    C1 * cross(q_values,n,2);
end

function [V,Omega,e_trans,e_rot,e_surface,u_mfs_check] = ...
        evaluate_solution(U,lambda,rvec_in,r_check,centre,opt, ...
        V_exact,Omega_exact,u_slip_check_vec,slip_check_norm)
% Evaluate one numerical solution only against the exact swimming data and boundary condition.
V = U(1:3);
Omega = U(4:6);
u_mfs_check = getStokesletFlow(lambda,rvec_in,r_check,opt);
u_rigid_check = getKmat(r_check,centre) * U;
boundary_residual = ...
    u_mfs_check - u_rigid_check - u_slip_check_vec;

e_trans = norm(V-V_exact,2) / norm(V_exact,2);
e_rot = norm(Omega-Omega_exact,2) / norm(Omega_exact,2);
e_surface = norm(boundary_residual,inf) / slip_check_norm;
end

function validate_convergence(method_results,method_name,maxit)
% Enforce independent exact-solution and boundary-residual acceptance criteria.
all_values = [ ...
    method_results.V(:); ...
    method_results.Omega(:); ...
    method_results.error.trans(:); ...
    method_results.error.rot(:); ...
    method_results.error.surface(:); ...
    method_results.lambda_total_norm(:); ...
    method_results.solve_time(:)];
if any(~isfinite(all_values))
    error('jul16_single_sphere_slip:nonfiniteResult', ...
        '%s produced a nonfinite result.',method_name);
end
if any(method_results.iters >= maxit)
    error('jul16_single_sphere_slip:maxIterations', ...
        '%s reached the maximum GMRES iteration count.',method_name);
end

error_fields = {'trans','rot','surface'};
for i_error = 1:numel(error_fields)
    error_name = error_fields{i_error};
    error_values = method_results.error.(error_name);
    if error_values(end) >= error_values(1)
        error('jul16_single_sphere_slip:noConvergence', ...
            '%s %s error did not decrease from the coarsest grid.', ...
            method_name,error_name);
    end
    if error_values(end) > 1e-6
        error('jul16_single_sphere_slip:accuracyTarget', ...
            '%s %s error %.3e exceeds 1e-6 on the finest grid.', ...
            method_name,error_name,error_values(end));
    end
end
end

function R = axis_angle_rotation(axis,angle)
% Build a proper rotation matrix using Rodrigues' formula.
axis = axis(:) / norm(axis);
K = [0 -axis(3) axis(2); ...
     axis(3) 0 -axis(1); ...
     -axis(2) axis(1) 0];
R = eye(3) + sin(angle)*K + (1-cos(angle))*(K*K);
end

function plot_surface_vector_field( ...
        points,u_vec,arrow_stride,figure_name,title_text,colorbar_label)
% Plot a surface vector field using coloured nodes and subsampled arrows.
u = reshape(u_vec,3,[]).';
speed = vecnorm(u,2,2);
arrow_indices = 1:arrow_stride:size(points,1);

figure('Color','w','Name',figure_name);
ax = axes();
scatter3(ax,0.95*points(:,1),0.95*points(:,2),0.95*points(:,3),30,speed,'filled');
hold(ax,'on');
quiver3(ax, ...
    points(arrow_indices,1), ...
    points(arrow_indices,2), ...
    points(arrow_indices,3), ...
    u(arrow_indices,1), ...
    u(arrow_indices,2), ...
    u(arrow_indices,3), ...
    0.8,'k','LineWidth',1);

axis(ax,'equal');
axis(ax,'tight');
grid(ax,'on');
box(ax,'on');
view(ax,3);
colormap(ax,parula);
set(ax,'TickLabelInterpreter','latex','FontSize',12);
xlabel(ax,'$x$','Interpreter','latex');
ylabel(ax,'$y$','Interpreter','latex');
zlabel(ax,'$z$','Interpreter','latex');
title(ax,title_text,'Interpreter','latex');

cb = colorbar(ax);
cb.TickLabelInterpreter = 'latex';
cb.Label.String = colorbar_label;
cb.Label.Interpreter = 'latex';
end
