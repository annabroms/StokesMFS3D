%AUG05_SINGLE_SPHERE_SLIP_COLLOCATION_CONVERGENCE
%
% This experiment is based on the same physical and numerical setting as
% jul16_single_sphere_slip_convergence.m: a force- and torque-free unit
% sphere has prescribed non-axisymmetric tangential slip, and both the
% standard and DLP mobility solvers are compared with exact swimming
% velocities and an independently evaluated surface residual.
%
% Unlike the Jul 16 source-convergence test, this script fixes each of two
% proxy source grids (requested N = 200 and 400) in turn and varies only the
% number M of physical surface collocation/quadrature nodes. Only designs
% with actual M >= actual N are solved. Spherical designs are not available
% at every integer, so both requested and actual N and M are reported.
% A separate DLP-only comparison repeats the admissible cases for proxy
% radii Rp = 0.3 and 0.6.
%
% The boundary convention is
%
%   u = V + Omega x (x-c) + u_slip,
%
% and the signs of the slip modes below give positive exact velocities
% along p_dir and q_dir. The surface residual is measured on the same type
% of fixed, finer, rotated check grid used by the Jul 16 experiment.
%
% Run from any working directory with
%
%   run('experiments/dev_DLP/aug05_single_sphere_slip_collocation_convergence.m')
%
% Created August 5, 2026.

clear;
close all;
clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
repo_root = find_repo_root(script_dir);
run(fullfile(repo_root,'setpath.m'));

%% Physical configuration: identical to the Jul 16 test
centre = [0 0 0];
radius = 1;
B1 = 1;
C1 = 1;
p_dir = [1; 2; 3]/sqrt(14);
q_dir = [-2; 1; 1]/sqrt(6);

V_exact = (2*B1/3)*p_dir;
Omega_exact = (C1/radius)*q_dir;
Fvec = zeros(6,1);

%% Two fixed source grids and variable collocation grids
N_request_vec = [200 300];
M_request_vec = 100:50:1000;
Rp = 0.6;
Rp_dlp_vec = [0.3 0.6];

n_N = numel(N_request_vec);
n_M = numel(M_request_vec);
source_directions_by_N = cell(1,n_N);
N_actual_vec = nan(1,n_N);
for iN = 1:n_N
    [source_directions_by_N{iN},~] = ...
        get_sphdesign(N_request_vec(iN));
    N_actual_vec(iN) = size(source_directions_by_N{iN},1);
end

surface_directions_by_M = cell(1,n_M);
surface_weights_by_M = cell(1,n_M);
M_actual_vec = nan(1,n_M);
for iM = 1:n_M
    [surface_directions_by_M{iM},surface_weights_by_M{iM}] = ...
        get_sphdesign(M_request_vec(iM));
    M_actual_vec(iM) = size(surface_directions_by_M{iM},1);
end
t_actual_vec = get_spherical_design_degrees(M_actual_vec);
if any(diff(N_actual_vec) <= 0) || any(diff(M_actual_vec) <= 0)
    error('aug05_single_sphere_slip:duplicateDesign', ...
        ['Requested N and M values must produce strictly increasing ', ...
         'actual spherical-design sizes.']);
end
n_Rp_dlp = numel(Rp_dlp_vec);
baseline_Rp_index = find(abs(Rp_dlp_vec-Rp) <= ...
    32*eps(max(1,abs(Rp))),1);
if isempty(baseline_Rp_index)
    error('aug05_single_sphere_slip:missingBaselineProxyRadius', ...
        'Rp_dlp_vec must contain the main proxy radius Rp.');
end

%% Solver configuration: identical to the Jul 16 test
gmres_tol = 1e-14;
maxit = 600;
use_fmm = false;
fmm_tol = 1e-12;

dlp_inner_only = false;
dlp_symmetrize_weighted = false;
dlp_add_rank1 = false;
dlp_outer_force = false;

%% Independent rotated surface-check grid
N_check_request = 5000;
[n_check_base,~] = get_sphdesign(N_check_request);
rotation_axis = [2; -3; 1]/sqrt(14);
rotation_angle = pi/7;
R_check = axis_angle_rotation(rotation_axis,rotation_angle);

n_check = (R_check*n_check_base.').';
r_check = centre + radius*n_check;
u_slip_check = evaluate_slip(n_check,B1,C1,p_dir,q_dir);
u_slip_check_vec = reshape(u_slip_check.',[],1);
slip_check_norm = norm(u_slip_check_vec,inf);

%% Storage
results = struct();
results.created = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));
results.description = [ ...
    'Jul 16 single-sphere prescribed-slip setting with fixed source ', ...
    'count and varying collocation/quadrature count.'];
results.reference_experiment = ...
    'experiments/dev_DLP/jul16_single_sphere_slip_convergence.m';
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
results.config.M_request = M_request_vec;
results.config.N_check_request = N_check_request;
results.config.N_check_actual = size(r_check,1);
results.config.check_rotation_axis = rotation_axis;
results.config.check_rotation_angle = rotation_angle;
results.config.Rp = Rp;
results.config.Rp_dlp = Rp_dlp_vec;
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
results.actual.N = N_actual_vec;
results.actual.M = M_actual_vec;
results.actual.t = t_actual_vec;
results.actual.valid_M = M_actual_vec >= N_actual_vec(:);
results.actual.a = nan(n_N,n_M);
results.tangential_defect = nan(n_N,n_M);

method_results.V = nan(3,n_N,n_M);
method_results.Omega = nan(3,n_N,n_M);
method_results.error.trans = nan(n_N,n_M);
method_results.error.rot = nan(n_N,n_M);
method_results.error.surface = nan(n_N,n_M);
method_results.iters = nan(n_N,n_M);
method_results.lambda_total_norm = nan(n_N,n_M);
method_results.solve_time = nan(n_N,n_M);
results.standard = method_results;
results.dlp = method_results;

results.dlp_proxy.Rp = Rp_dlp_vec;
results.dlp_proxy.V = nan(3,n_Rp_dlp,n_N,n_M);
results.dlp_proxy.Omega = nan(3,n_Rp_dlp,n_N,n_M);
results.dlp_proxy.error.trans = nan(n_Rp_dlp,n_N,n_M);
results.dlp_proxy.error.rot = nan(n_Rp_dlp,n_N,n_M);
results.dlp_proxy.error.surface = nan(n_Rp_dlp,n_N,n_M);
results.dlp_proxy.iters = nan(n_Rp_dlp,n_N,n_M);
results.dlp_proxy.lambda_total_norm = nan(n_Rp_dlp,n_N,n_M);
results.dlp_proxy.solve_time = nan(n_Rp_dlp,n_N,n_M);

fprintf('Single-sphere prescribed-slip collocation sweep\n');
fprintf('  based on jul16_single_sphere_slip_convergence.m\n');
fprintf('  requested N = [%s]\n',num2str(N_request_vec));
fprintf('  actual N = [%s], Rp = %.3g\n',num2str(N_actual_vec),Rp);
fprintf('  DLP-only proxy-radius comparison Rp = [%s]\n', ...
    num2str(Rp_dlp_vec));
fprintf('  requested M = [%s]\n',num2str(M_request_vec));
fprintf('  actual M = [%s]\n',num2str(M_actual_vec));
fprintf('  only cases satisfying actual M >= actual N are solved\n');
fprintf('  independent rotated check nodes = %d\n',size(r_check,1));
fprintf('  gmres_tol = %.1e, maxit = %d\n\n',gmres_tol,maxit);
fprintf(['%9s %6s %6s %8s %6s %7s %11s %11s %11s ', ...
    '%11s %9s\n'], ...
    'method','N','M','M/N','iters','time','e_trans','e_rot', ...
    'e_res_max','|lambda|inf','n dot us');

%% Collocation-count sweeps for both fixed source grids
for iN = 1:n_N
    N_request = N_request_vec(iN);
    N_actual = N_actual_vec(iN);
    rvec_in = centre + radius*Rp*source_directions_by_N{iN};
    valid_M_indices = find(results.actual.valid_M(iN,:));

    fprintf('\nFixed requested N = %d, actual N = %d\n', ...
        N_request,N_actual);
    for iM = valid_M_indices
        surface_directions = surface_directions_by_M{iM};
        surface_weights = surface_weights_by_M{iM};
        rvec_out = centre + radius*surface_directions;
        nout = surface_directions;
        wout = radius^2*surface_weights;
        M_actual = M_actual_vec(iM);

        u_slip = evaluate_slip(nout,B1,C1,p_dir,q_dir);
        tangential_defect = max(abs(sum(nout.*u_slip,2)));
        if tangential_defect > 1e-12
            error('aug05_single_sphere_slip:tangentiality', ...
                'Slip tangentiality defect is %.3e at N = %d, M = %d.', ...
                tangential_defect,N_actual,M_actual);
        end

        opt = configure_options( ...
            N_request,Rp,gmres_tol,maxit,use_fmm,fmm_tol, ...
            dlp_inner_only,dlp_symmetrize_weighted, ...
            dlp_add_rank1,dlp_outer_force);

        standard_start = tic;
        [U_standard,iters_standard,~,~,lambda_standard] = ...
            solve_mobility( ...
            centre,rvec_in,rvec_out,Fvec,u_slip,opt);
        solve_time_standard = toc(standard_start);

        dlp_start = tic;
        [U_dlp,iters_dlp,~,~,lambda_dlp] = ...
            solve_mobility_with_DLP( ...
            centre,rvec_in,rvec_out,nout,wout,Fvec,u_slip,opt);
        solve_time_dlp = toc(dlp_start);

        [V_standard,Omega_standard,e_trans_standard,e_rot_standard, ...
            e_surface_standard] = evaluate_solution( ...
            U_standard,lambda_standard,rvec_in,r_check,centre,opt, ...
            V_exact,Omega_exact,u_slip_check_vec,slip_check_norm);
        [V_dlp,Omega_dlp,e_trans_dlp,e_rot_dlp,e_surface_dlp] = ...
            evaluate_solution( ...
            U_dlp,lambda_dlp,rvec_in,r_check,centre,opt, ...
            V_exact,Omega_exact,u_slip_check_vec,slip_check_norm);

        results.actual.a(iN,iM) = M_actual/N_actual;
        results.tangential_defect(iN,iM) = tangential_defect;
        results.standard.V(:,iN,iM) = V_standard;
        results.standard.Omega(:,iN,iM) = Omega_standard;
        results.standard.error.trans(iN,iM) = e_trans_standard;
        results.standard.error.rot(iN,iM) = e_rot_standard;
        results.standard.error.surface(iN,iM) = e_surface_standard;
        results.standard.iters(iN,iM) = iters_standard;
        results.standard.lambda_total_norm(iN,iM) = ...
            norm(lambda_standard,inf);
        results.standard.solve_time(iN,iM) = solve_time_standard;
        results.dlp.V(:,iN,iM) = V_dlp;
        results.dlp.Omega(:,iN,iM) = Omega_dlp;
        results.dlp.error.trans(iN,iM) = e_trans_dlp;
        results.dlp.error.rot(iN,iM) = e_rot_dlp;
        results.dlp.error.surface(iN,iM) = e_surface_dlp;
        results.dlp.iters(iN,iM) = iters_dlp;
        results.dlp.lambda_total_norm(iN,iM) = norm(lambda_dlp,inf);
        results.dlp.solve_time(iN,iM) = solve_time_dlp;

        results.dlp_proxy.V(:,baseline_Rp_index,iN,iM) = V_dlp;
        results.dlp_proxy.Omega(:,baseline_Rp_index,iN,iM) = Omega_dlp;
        results.dlp_proxy.error.trans(baseline_Rp_index,iN,iM) = ...
            e_trans_dlp;
        results.dlp_proxy.error.rot(baseline_Rp_index,iN,iM) = ...
            e_rot_dlp;
        results.dlp_proxy.error.surface(baseline_Rp_index,iN,iM) = ...
            e_surface_dlp;
        results.dlp_proxy.iters(baseline_Rp_index,iN,iM) = iters_dlp;
        results.dlp_proxy.lambda_total_norm( ...
            baseline_Rp_index,iN,iM) = norm(lambda_dlp,inf);
        results.dlp_proxy.solve_time( ...
            baseline_Rp_index,iN,iM) = solve_time_dlp;

        validate_case( ...
            results.standard,iN,iM,'standard',N_actual,M_actual,maxit);
        validate_case(results.dlp,iN,iM,'DLP',N_actual,M_actual,maxit);

        for iRp = 1:n_Rp_dlp
            if iRp == baseline_Rp_index
                continue;
            end
            Rp_proxy = Rp_dlp_vec(iRp);
            rvec_in_proxy = centre + ...
                radius*Rp_proxy*source_directions_by_N{iN};
            opt_proxy = configure_options( ...
                N_request,Rp_proxy,gmres_tol,maxit,use_fmm,fmm_tol, ...
                dlp_inner_only,dlp_symmetrize_weighted, ...
                dlp_add_rank1,dlp_outer_force);

            proxy_start = tic;
            [U_proxy,iters_proxy,~,~,lambda_proxy] = ...
                solve_mobility_with_DLP( ...
                centre,rvec_in_proxy,rvec_out,nout,wout,Fvec, ...
                u_slip,opt_proxy);
            solve_time_proxy = toc(proxy_start);
            [V_proxy,Omega_proxy,e_trans_proxy,e_rot_proxy, ...
                e_surface_proxy] = evaluate_solution( ...
                U_proxy,lambda_proxy,rvec_in_proxy,r_check,centre, ...
                opt_proxy,V_exact,Omega_exact,u_slip_check_vec, ...
                slip_check_norm);

            results.dlp_proxy.V(:,iRp,iN,iM) = V_proxy;
            results.dlp_proxy.Omega(:,iRp,iN,iM) = Omega_proxy;
            results.dlp_proxy.error.trans(iRp,iN,iM) = e_trans_proxy;
            results.dlp_proxy.error.rot(iRp,iN,iM) = e_rot_proxy;
            results.dlp_proxy.error.surface(iRp,iN,iM) = ...
                e_surface_proxy;
            results.dlp_proxy.iters(iRp,iN,iM) = iters_proxy;
            results.dlp_proxy.lambda_total_norm(iRp,iN,iM) = ...
                norm(lambda_proxy,inf);
            results.dlp_proxy.solve_time(iRp,iN,iM) = ...
                solve_time_proxy;

            validate_proxy_case( ...
                results.dlp_proxy,iRp,iN,iM,Rp_proxy, ...
                N_actual,M_actual,maxit);
            fprintf(['DLP proxy Rp=%3.1f, N=%d, M=%d: ', ...
                'e_trans=%.3e, e_rot=%.3e, e_res=%.3e\n'], ...
                Rp_proxy,N_actual,M_actual,e_trans_proxy,e_rot_proxy, ...
                e_surface_proxy);
        end

        fprintf(['%9s %6d %6d %8.3f %6d %7.2f %11.3e %11.3e ', ...
            '%11.3e %11.3e %9.2e\n'], ...
            'standard',N_actual,M_actual,M_actual/N_actual, ...
            iters_standard,solve_time_standard,e_trans_standard, ...
            e_rot_standard,e_surface_standard, ...
            results.standard.lambda_total_norm(iN,iM), ...
            tangential_defect);
        fprintf(['%9s %6d %6d %8.3f %6d %7.2f %11.3e %11.3e ', ...
            '%11.3e %11.3e %9.2e\n'], ...
            'DLP',N_actual,M_actual,M_actual/N_actual, ...
            iters_dlp,solve_time_dlp,e_trans_dlp,e_rot_dlp, ...
            e_surface_dlp,results.dlp.lambda_total_norm(iN,iM), ...
            tangential_defect);
    end
end

%% Figures: the same metrics versus M and versus a=M/N
metric_names = {'trans','rot','surface'};
metric_titles = { ...
    '$\epsilon_{\mathrm{trans}}$', ...
    '$\epsilon_{\mathrm{rot}}$', ...
    '$\epsilon_{\mathrm{res}}^{\max}$'};
method_colors = [ ...
    0.0000 0.4470 0.7410; ...
    0.8500 0.3250 0.0980];
method_markers = {'o','s'};
N_line_styles = {'-','--'};

plot_error_panels( ...
    results,metric_names,metric_titles,method_colors,method_markers, ...
    N_line_styles,'M');
plot_error_panels( ...
    results,metric_names,metric_titles,method_colors,method_markers, ...
    N_line_styles,'a');

%% Figure: DLP-only proxy-radius comparison versus a=M/N
Rp_colors = [ ...
    0.0000 0.4470 0.7410; ...
    0.8500 0.3250 0.0980];
Rp_markers = {'o','s'};
plot_dlp_proxy_panels( ...
    results,metric_names,metric_titles,Rp_colors,Rp_markers, ...
    N_line_styles);

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
        error('aug05_single_sphere_slip:repoRootNotFound', ...
            'Could not find setpath.m starting from %s.',start_dir);
    end
    repo_root = parent_dir;
end
end

function t_values = get_spherical_design_degrees(M_values)
% Recover the exact polynomial degree associated with each loaded design.
design_dir = fullfile(fileparts(which('get_sphdesign')), ...
    'spherical_designs');
file_id = fopen(fullfile(design_dir,'filelist.txt'),'r');
if file_id < 0
    error('aug05_single_sphere_slip:designCatalogueNotFound', ...
        'Could not open the spherical-design catalogue.');
end
file_cleanup = onCleanup(@() fclose(file_id));
design_data = textscan(file_id,'%3d.%5d');
available_t = design_data{1};
available_M = design_data{2};
[found,indices] = ismember(M_values(:),available_M);
if ~all(found)
    error('aug05_single_sphere_slip:designDegreeNotFound', ...
        'An actual M value is absent from the spherical-design catalogue.');
end
t_values = reshape(available_t(indices),size(M_values));
clear file_cleanup;
end

function opt = configure_options( ...
        N_request,Rp,gmres_tol,maxit,use_fmm,fmm_tol, ...
        inner_only,symmetrize_weighted,add_rank1,outer_force)
% Build the common solver options without regenerating either node set.
opt = init_MFS(N_request);
opt.Rp = Rp;
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

function u_slip = evaluate_slip(n,B1,C1,p_dir,q_dir)
% Evaluate the Jul 16 tangential translational and rotational slip modes.
n_nodes = size(n,1);
p_values = repmat(p_dir.',n_nodes,1);
q_values = repmat(q_dir.',n_nodes,1);
u_slip = -B1*(p_values-(n*p_dir).*n) - ...
    C1*cross(q_values,n,2);
end

function [V,Omega,e_trans,e_rot,e_surface] = evaluate_solution( ...
        U,lambda,rvec_in,r_check,centre,opt,V_exact,Omega_exact, ...
        u_slip_check_vec,slip_check_norm)
% Evaluate exact rigid velocities and the independent maximum-norm residual.
V = U(1:3);
Omega = U(4:6);
u_mfs_check = getStokesletFlow(lambda,rvec_in,r_check,opt);
u_rigid_check = getKmat(r_check,centre)*U;
boundary_residual = ...
    u_mfs_check-u_rigid_check-u_slip_check_vec;

e_trans = norm(V-V_exact,2)/norm(V_exact,2);
e_rot = norm(Omega-Omega_exact,2)/norm(Omega_exact,2);
e_surface = norm(boundary_residual,inf)/slip_check_norm;
end

function plot_error_panels( ...
        results,metric_names,metric_titles,method_colors,method_markers, ...
        N_line_styles,x_mode)
% Plot the three errors with method and source-count encodings kept separate.
n_N = numel(results.actual.N);
if numel(N_line_styles) ~= n_N
    error('aug05_single_sphere_slip:badLineStyles', ...
        'Supply one line style for each fixed source count.');
end

switch x_mode
    case 'M'
        figure_name = 'Single-sphere slip errors versus M';
        x_label = '$M$';
    case 'a'
        figure_name = 'Single-sphere slip errors versus M over N';
        x_label = '$a=M/N$';
    otherwise
        error('aug05_single_sphere_slip:badPlotAbscissa', ...
            'x_mode must be ''M'' or ''a''.');
end

figure('Color','w','Name',figure_name);
layout = tiledlayout( ...
    1,numel(metric_names),'TileSpacing','compact','Padding','compact');
axes_handles = gobjects(1,numel(metric_names));
for i_metric = 1:numel(metric_names)
    ax = nexttile(layout);
    axes_handles(i_metric) = ax;

    metric_name = metric_names{i_metric};
    for iN = 1:n_N
        valid_M = results.actual.valid_M(iN,:);
        if strcmp(x_mode,'M')
            x_values = results.actual.M(valid_M);
        else
            x_values = results.actual.a(iN,valid_M);
        end
        semilogy(ax,x_values, ...
            results.standard.error.(metric_name)(iN,valid_M), ...
            'Color',method_colors(1,:), ...
            'Marker',method_markers{1}, ...
            'LineStyle',N_line_styles{iN}, ...
            'LineWidth',1.5,'MarkerSize',7, ...
            'HandleVisibility','off');
        hold(ax,'on');
        semilogy(ax,x_values, ...
            results.dlp.error.(metric_name)(iN,valid_M), ...
            'Color',method_colors(2,:), ...
            'Marker',method_markers{2}, ...
            'LineStyle',N_line_styles{iN}, ...
            'LineWidth',1.5,'MarkerSize',7, ...
            'HandleVisibility','off');
    end
    axis(ax,'tight');
    format_axes(ax);
    xlabel(ax,x_label,'Interpreter','latex');
    ylabel(ax,'$\mathrm{relative\ error}$','Interpreter','latex');
    title(ax,metric_titles{i_metric},'Interpreter','latex');
end

legend_handles = gobjects(1,2+n_N);
legend_handles(1) = semilogy(axes_handles(1),nan,nan, ...
    'Color',method_colors(1,:), ...
    'Marker',method_markers{1},'LineStyle','none', ...
    'MarkerSize',7,'DisplayName','$\mathrm{standard}$');
legend_handles(2) = semilogy(axes_handles(1),nan,nan, ...
    'Color',method_colors(2,:), ...
    'Marker',method_markers{2},'LineStyle','none', ...
    'MarkerSize',7,'DisplayName','$\mathrm{DLP}$');
for iN = 1:n_N
    legend_handles(2+iN) = semilogy(axes_handles(1),nan,nan, ...
        'Color','k','LineStyle',N_line_styles{iN}, ...
        'LineWidth',1.5,'DisplayName',sprintf( ...
        '$N=%d\;(N_{\mathrm{req}}=%d)$', ...
        results.actual.N(iN),results.config.N_request(iN)));
end
plot_legend = legend( ...
    axes_handles(1),legend_handles,'Interpreter','latex', ...
    'NumColumns',2);
plot_legend.Layout.Tile = 'south';
end

function plot_dlp_proxy_panels( ...
        results,metric_names,metric_titles,Rp_colors,Rp_markers, ...
        N_line_styles)
% Plot DLP errors with proxy radius encoded by color and marker.
n_Rp = numel(results.dlp_proxy.Rp);
n_N = numel(results.actual.N);
if size(Rp_colors,1) ~= n_Rp || numel(Rp_markers) ~= n_Rp || ...
        numel(N_line_styles) ~= n_N
    error('aug05_single_sphere_slip:badProxyPlotStyles', ...
        ['Supply one color and marker per proxy radius and one line ', ...
         'style per N.']);
end

figure('Color','w','Name','DLP proxy-radius comparison versus M');
layout = tiledlayout( ...
    1,numel(metric_names),'TileSpacing','compact','Padding','compact');
axes_handles = gobjects(1,numel(metric_names));
for i_metric = 1:numel(metric_names)
    ax = nexttile(layout);
    axes_handles(i_metric) = ax;
    metric_name = metric_names{i_metric};
    for iRp = 1:n_Rp
        for iN = 1:n_N
            valid_M = results.actual.valid_M(iN,:);
            x_values = results.actual.M(valid_M);
            error_values = reshape( ...
                results.dlp_proxy.error.(metric_name)( ...
                iRp,iN,valid_M),1,[]);
            semilogy(ax,x_values,error_values, ...
                'Color',Rp_colors(iRp,:), ...
                'Marker',Rp_markers{iRp}, ...
                'LineStyle',N_line_styles{iN}, ...
                'LineWidth',1.5,'MarkerSize',7, ...
                'HandleVisibility','off');
             hold(ax,'on');
        end

        % Anchor the expected C*Rp^t quadrature scaling to the first
        % positive point of the higher-N curve in this panel.
        guide_N_index = n_N;
        guide_valid_M = results.actual.valid_M(1,:);
        guide_M = results.actual.M(guide_valid_M);
        %guide_t = results.actual.t(guide_valid_M);

        guide_t = double(results.actual.t(guide_valid_M));
       
        %help_line = 1e7*results.dlp_proxy.Rp(iRp).^( ...
            %    double(guide_t));
        % semilogy(ax,guide_M,help_line,':', ...
        %         'Color',Rp_colors(iRp,:), ...
        %         'LineWidth',1.5,'HandleVisibility','off');



        % guide_error = reshape( ...
        %     results.dlp_proxy.error.(metric_name)( ...
        %     iRp,guide_N_index,guide_valid_M),1,[]);
        % anchor_index = find(isfinite(guide_error) & guide_error > 0,1);
        % if ~isempty(anchor_index)
        %     help_line = guide_error(anchor_index) * ...
        %         results.dlp_proxy.Rp(iRp).^( ...
        %         guide_t-guide_t(anchor_index));
        %     semilogy(ax,guide_M,help_line,':', ...
        %         'Color','k', ...
        %         'LineWidth',1.5,'HandleVisibility','off');
        % end
    end
    axis(ax,'tight');
    format_axes(ax);
    xlabel(ax,'$M$','Interpreter','latex');
    ylabel(ax,metric_titles{i_metric},'Interpreter','latex');
end

legend_handles = gobjects(1,n_Rp+n_N+1);
for iRp = 1:n_Rp
    legend_handles(iRp) = semilogy(axes_handles(1),nan,nan, ...
        'Color',Rp_colors(iRp,:), ...
        'Marker',Rp_markers{iRp},'LineStyle','none','MarkerSize',7, ...
        'DisplayName',sprintf('$R_p=%g$',results.dlp_proxy.Rp(iRp)));
end
for iN = 1:n_N
    legend_handles(n_Rp+iN) = semilogy(axes_handles(1),nan,nan, ...
        'Color','k','LineStyle',N_line_styles{iN}, ...
        'LineWidth',1.5,'DisplayName',sprintf( ...
        '$N=%d\;(N_{\mathrm{req}}=%d)$', ...
        results.actual.N(iN),results.config.N_request(iN)));
end
legend_handles(end) = semilogy(axes_handles(1),nan,nan, ...
    'Color','k','LineStyle',':','LineWidth',1.5, ...
    'DisplayName','$C R_p^t\sim C\exp(-c\sqrt{M})$');
proxy_legend = legend( ...
    axes_handles(1),legend_handles,'Interpreter','latex', ...
    'NumColumns',2);
proxy_legend.Layout.Tile = 'south';
end

function validate_case( ...
        method_results,iN,iM,method_name,N_actual,M_actual,maxit)
% Reject nonfinite solves while allowing a nonmonotone oversampling sweep.
case_values = [ ...
    method_results.V(:,iN,iM); ...
    method_results.Omega(:,iN,iM); ...
    method_results.error.trans(iN,iM); ...
    method_results.error.rot(iN,iM); ...
    method_results.error.surface(iN,iM); ...
    method_results.lambda_total_norm(iN,iM); ...
    method_results.solve_time(iN,iM)];
if any(~isfinite(case_values))
    error('aug05_single_sphere_slip:nonfiniteResult', ...
        '%s produced a nonfinite result at N = %d, M = %d.', ...
        method_name,N_actual,M_actual);
end
if method_results.iters(iN,iM) >= maxit
    warning('aug05_single_sphere_slip:maxIterations', ...
        '%s reached maxit at N = %d, M = %d.', ...
        method_name,N_actual,M_actual);
end
end

function validate_proxy_case( ...
        proxy_results,iRp,iN,iM,Rp_proxy,N_actual,M_actual,maxit)
% Reject nonfinite DLP proxy-radius cases without imposing monotonicity.
case_values = [ ...
    proxy_results.V(:,iRp,iN,iM); ...
    proxy_results.Omega(:,iRp,iN,iM); ...
    proxy_results.error.trans(iRp,iN,iM); ...
    proxy_results.error.rot(iRp,iN,iM); ...
    proxy_results.error.surface(iRp,iN,iM); ...
    proxy_results.lambda_total_norm(iRp,iN,iM); ...
    proxy_results.solve_time(iRp,iN,iM)];
if any(~isfinite(case_values))
    error('aug05_single_sphere_slip:nonfiniteProxyResult', ...
        ['DLP produced a nonfinite result at Rp = %.3g, ', ...
         'N = %d, M = %d.'],Rp_proxy,N_actual,M_actual);
end
if proxy_results.iters(iRp,iN,iM) >= maxit
    warning('aug05_single_sphere_slip:proxyMaxIterations', ...
        'DLP reached maxit at Rp = %.3g, N = %d, M = %d.', ...
        Rp_proxy,N_actual,M_actual);
end
end

function R = axis_angle_rotation(axis,angle)
% Build the fixed proper rotation used for the independent check grid.
axis = axis(:)/norm(axis);
K = [0 -axis(3) axis(2); ...
     axis(3) 0 -axis(1); ...
     -axis(2) axis(1) 0];
R = eye(3)+sin(angle)*K+(1-cos(angle))*(K*K);
end

function format_axes(ax)
% Apply the experiment's common LaTeX and grid formatting.
grid(ax,'on');
box(ax,'on');
set(ax,'TickLabelInterpreter','latex','FontSize',12);
end
