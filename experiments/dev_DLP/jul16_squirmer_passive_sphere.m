%JUL16_SQUIRMER_PASSIVE_SPHERE
%
% Validate solve_mobility_with_DLP for a coaxial unit spherical squirmer
% interacting with a passive unit sphere. The exact reference is the
% bispherical solution of:
%
%   D. Papavassiliou and G. P. Alexander,
%   "Exact solutions for hydrodynamic interactions of two squirming
%   spheres", J. Fluid Mech. 813 (2017), 618--646.
%   https://wrap.warwick.ac.uk/id/eprint/87388/1/
%       WRAP-exact-solutions-hydrodynamic-squirming-papavassiliou-2017.pdf
%
% The exact map uses equations (3.12)--(3.13), (3.16), (4.1), and
% Appendix B equations (B3)--(B4). The repository slip convention is
%
%   u = V + Omega x (x-c) + u_slip,
%
% where u_slip is the physical surface-distortion velocity denoted u_s in
% the paper.
%
%
% This script performs an N-convergence sweep at h_convergence and a gap sweep at
% N_request_gap. Change those top-level parameters below as needed.
%
% Created July 16, 2026.

clear;
close all;
clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
repo_root = find_repo_root(script_dir);
run(fullfile(repo_root,'setpath.m'));

%% Editable experiment parameters
h_convergence = 0.2;
N_request_convergence = [100 200 400 700 1000 1500];
gap_vec = [4 2 1 0.5 0.2 0.1 0.05];
N_request_gap = 1000;

beta_vec = [-1 0 1];
a = 1;
B1 = 3/2;
U0 = 2*B1/3;

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

series_opts = struct( ...
    'tol',1e-13, ...
    'min_terms',20, ...
    'consecutive_terms',8, ...
    'max_terms',2000);

%% Independent rotated surface-check grid
N_check_request = 5000;
[n_check_base,~] = get_sphdesign(N_check_request);
rotation_axis = [2; -1; 3]/sqrt(14);
rotation_angle = pi/11;
R_check = axis_angle_rotation(rotation_axis,rotation_angle);
n_check = (R_check*n_check_base.').';

%% Exact-reference validation and storage
[mode_map_convergence,info_convergence] = ...
    validated_exact_map(h_convergence,a,series_opts);
exact_convergence = mode_map_convergence * ...
    [B1*ones(1,numel(beta_vec)); B1*beta_vec];

far_h = 1000;
[far_map,far_info] = papavassiliou_squirmer_passive_map( ...
    far_h,a,series_opts);
far_response = far_map*[B1;0];
far_field_defect = norm(far_response-[U0;0],inf)/U0;
if far_field_defect > 1e-6
    error('jul16_squirmer_passive_sphere:farFieldValidation', ...
        'Exact far-field validation failed with defect %.3e.', ...
        far_field_defect);
end

linearity_coefficients_a = [B1; -0.3*B1];
linearity_coefficients_b = [-0.2*B1; 0.7*B1];
linearity_defect = norm( ...
    mode_map_convergence * ...
        (linearity_coefficients_a+linearity_coefficients_b) - ...
    mode_map_convergence*linearity_coefficients_a - ...
    mode_map_convergence*linearity_coefficients_b,inf);
if linearity_defect > 100*eps*max(1,norm(mode_map_convergence,inf))
    error('jul16_squirmer_passive_sphere:linearityValidation', ...
        'Exact-map linearity defect is %.3e.',linearity_defect);
end

%% Results metadata
results = struct();
results.created = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));
results.description = [ ...
    'Coaxial unit spherical squirmer interacting with a passive unit ', ...
    'sphere: DLP-MFS versus the exact bispherical series.'];
results.reference = [ ...
    'Papavassiliou and Alexander, J. Fluid Mech. 813 (2017), ', ...
    'equations (3.12)-(3.13), (3.16), (4.1), and (B3)-(B4).'];
results.config = struct();
results.config.h_convergence = h_convergence;
results.config.N_request_convergence = N_request_convergence;
results.config.gap_vec = gap_vec;
results.config.N_request_gap = N_request_gap;
results.config.beta_vec = beta_vec;
results.config.a = a;
results.config.B1 = B1;
results.config.U0 = U0;
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
results.config.series_opts = series_opts;
results.config.N_check_request = N_check_request;
results.config.N_check_actual = size(n_check,1);
results.config.check_rotation_axis = rotation_axis;
results.config.check_rotation_angle = rotation_angle;
results.validation.exact.linearity_defect = linearity_defect;
results.validation.exact.far_h = far_h;
results.validation.exact.far_response = far_response;
results.validation.exact.far_field_defect = far_field_defect;
results.validation.exact.far_info = far_info;

%% Convergence sweep storage
n_beta = numel(beta_vec);
n_N = numel(N_request_convergence);
results.convergence.h = h_convergence;
results.convergence.beta = beta_vec;
results.convergence.B2 = B1*beta_vec;
results.convergence.N_request = N_request_convergence;
results.convergence.actual.N = nan(1,n_N);
results.convergence.actual.M = nan(1,n_N);
results.convergence.exact.velocity = exact_convergence;
results.convergence.exact.mode_map = mode_map_convergence;
results.convergence.exact.info = info_convergence;
results.convergence.U = nan(12,n_beta,n_N);
results.convergence.velocity = nan(2,n_beta,n_N);
results.convergence.error.squirmer = nan(n_beta,n_N);
results.convergence.error.passive = nan(n_beta,n_N);
results.convergence.error.symmetry = nan(n_beta,n_N);
results.convergence.error.surface = nan(n_beta,n_N);
results.convergence.iters = nan(n_beta,n_N);
results.convergence.lambda_total_norm = nan(n_beta,n_N);
results.convergence.solve_time = nan(n_beta,n_N);
results.convergence.precond_build_time = nan(1,n_N);

Fvec = zeros(12,1);
shared_gap_precond = [];
shared_gap_NM = [];
[gap_N_actual,gap_M_actual] = actual_design_sizes( ...
    N_request_gap,a_glob);
predicted_convergence_NM = nan(n_N,2);
for iN = 1:n_N
    [predicted_convergence_NM(iN,1), ...
     predicted_convergence_NM(iN,2)] = actual_design_sizes( ...
        N_request_convergence(iN),a_glob);
end
gap_cache_index = find( ...
    predicted_convergence_NM(:,1) == gap_N_actual & ...
    predicted_convergence_NM(:,2) == gap_M_actual,1,'last');
convergence_execution_order = 1:n_N;
if ~isempty(gap_cache_index)
    convergence_execution_order = [ ...
        convergence_execution_order( ...
            convergence_execution_order ~= gap_cache_index), ...
        gap_cache_index];
end
results.config.convergence_execution_order = ...
    N_request_convergence(convergence_execution_order);
results.config.gap_cache_reused_from_convergence = ...
    ~isempty(gap_cache_index);

fprintf('Unit-sphere squirmer--passive sphere validation\n');
fprintf('  B1 = %.6g, beta = [%s], U0 = %.6g\n', ...
    B1,num2str(beta_vec),U0);
fprintf('  Rp = %.3g, a_glob = %.3g, gmres_tol = %.1e\n', ...
    Rp,a_glob,gmres_tol);
fprintf('  independent rotated check nodes per sphere = %d\n\n', ...
    size(n_check,1));

fprintf('N-convergence sweep at h = %.6g\n',h_convergence);
fprintf(['%6s %6s %6s %7s %11s %11s %11s %11s ', ...
    '%10s %10s %10s %10s %11s %8s %8s\n'], ...
    'beta','N','M','iters','V_sq','V_sq ex','V_ps','V_ps ex', ...
    'e_sq','e_ps','e_sym','e_bdry','|lambda|','time','pc time');

%% N-convergence sweep
for iN = convergence_execution_order
    N_request = N_request_convergence(iN);
    q = sphere_centres(h_convergence,a);
    [rvec_in,rvec_out,opt] = init_spheres( ...
        q,Rp,N_request,a_glob);
    N_actual = size(rvec_in,1)/2;
    M_actual = size(rvec_out,1)/2;
    [nout,wout] = sphere_normals_and_weights(q,rvec_out,a);
    opt = configure_options( ...
        opt,gmres_tol,maxit,use_fmm,fmm_tol, ...
        dlp_inner_only,dlp_symmetrize_weighted, ...
        dlp_add_rank1,dlp_outer_force);

    precond_start = tic;
    precond = build_complete_preconditioner( ...
        rvec_in(1:N_actual,:)-q(1,:), ...
        rvec_out(1:M_actual,:)-q(1,:), ...
        nout(1:M_actual,:),wout(1:M_actual),opt);
    precond_build_time = toc(precond_start);
    opt.precond = precond;

    results.convergence.actual.N(iN) = N_actual;
    results.convergence.actual.M(iN) = M_actual;
    results.convergence.precond_build_time(iN) = precond_build_time;

    if iN == gap_cache_index
        shared_gap_precond = precond;
        shared_gap_NM = [N_actual M_actual];
    end

    for ib = 1:n_beta
        beta = beta_vec(ib);
        B2 = beta*B1;
        exact_velocity = exact_convergence(:,ib);

        case_result = solve_squirmer_case( ...
            q,rvec_in,rvec_out,nout,wout,Fvec,B1,B2,U0, ...
            exact_velocity,opt,n_check,a);

        results.convergence.U(:,ib,iN) = case_result.U;
        results.convergence.velocity(:,ib,iN) = ...
            [case_result.V_sq; case_result.V_ps];
        results.convergence.error.squirmer(ib,iN) = ...
            case_result.e_sq;
        results.convergence.error.passive(ib,iN) = ...
            case_result.e_ps;
        results.convergence.error.symmetry(ib,iN) = ...
            case_result.e_sym;
        results.convergence.error.surface(ib,iN) = ...
            case_result.e_surface;
        results.convergence.iters(ib,iN) = case_result.iters;
        results.convergence.lambda_total_norm(ib,iN) = ...
            case_result.lambda_total_norm;
        results.convergence.solve_time(ib,iN) = ...
            case_result.solve_time;

        fprintf(['%6.1f %6d %6d %7d %11.4e %11.4e %11.4e ', ...
            '%11.4e %10.3e %10.3e %10.3e %10.3e %11.3e ', ...
            '%8.2f %8.2f\n'], ...
            beta,N_actual,M_actual,case_result.iters, ...
            case_result.V_sq,exact_velocity(1), ...
            case_result.V_ps,exact_velocity(2), ...
            case_result.e_sq,case_result.e_ps,case_result.e_sym, ...
            case_result.e_surface,case_result.lambda_total_norm, ...
            case_result.solve_time,precond_build_time);
    end
end
clear precond opt rvec_in rvec_out nout wout;

%% Gap-sweep exact references and storage
n_h = numel(gap_vec);
results.gap_sweep.h = gap_vec;
results.gap_sweep.beta = beta_vec;
results.gap_sweep.B2 = B1*beta_vec;
results.gap_sweep.N_request = N_request_gap;
results.gap_sweep.actual.N = nan(1,n_h);
results.gap_sweep.actual.M = nan(1,n_h);
results.gap_sweep.exact.velocity = nan(2,n_beta,n_h);
results.gap_sweep.exact.mode_map = nan(2,2,n_h);
results.gap_sweep.exact.info = cell(1,n_h);
results.gap_sweep.exact.tolerance_stability = nan(1,n_h);
results.gap_sweep.U = nan(12,n_beta,n_h);
results.gap_sweep.velocity = nan(2,n_beta,n_h);
results.gap_sweep.error.squirmer = nan(n_beta,n_h);
results.gap_sweep.error.passive = nan(n_beta,n_h);
results.gap_sweep.error.symmetry = nan(n_beta,n_h);
results.gap_sweep.error.surface = nan(n_beta,n_h);
results.gap_sweep.iters = nan(n_beta,n_h);
results.gap_sweep.lambda_total_norm = nan(n_beta,n_h);
results.gap_sweep.solve_time = nan(n_beta,n_h);
results.gap_sweep.precond_build_time = 0;

fprintf('\nGap sweep at requested N = %d\n',N_request_gap);
fprintf(['%6s %7s %6s %6s %7s %11s %11s %11s %11s ', ...
    '%10s %10s %10s %10s %11s %8s %8s\n'], ...
    'beta','h','N','M','iters','V_sq','V_sq ex','V_ps','V_ps ex', ...
    'e_sq','e_ps','e_sym','e_bdry','|lambda|','time','pc time');

%% Gap sweep
for ih = 1:n_h
    h = gap_vec(ih);
    [mode_map,series_info,tolerance_stability] = ...
        validated_exact_map(h,a,series_opts);
    exact_velocity_all = mode_map * ...
        [B1*ones(1,n_beta); B1*beta_vec];

    q = sphere_centres(h,a);
    [rvec_in,rvec_out,opt] = init_spheres( ...
        q,Rp,N_request_gap,a_glob);
    N_actual = size(rvec_in,1)/2;
    M_actual = size(rvec_out,1)/2;
    [nout,wout] = sphere_normals_and_weights(q,rvec_out,a);
    opt = configure_options( ...
        opt,gmres_tol,maxit,use_fmm,fmm_tol, ...
        dlp_inner_only,dlp_symmetrize_weighted, ...
        dlp_add_rank1,dlp_outer_force);

    if isempty(shared_gap_precond)
        precond_start = tic;
        shared_gap_precond = build_complete_preconditioner( ...
            rvec_in(1:N_actual,:)-q(1,:), ...
            rvec_out(1:M_actual,:)-q(1,:), ...
            nout(1:M_actual,:),wout(1:M_actual),opt);
        results.gap_sweep.precond_build_time = toc(precond_start);
        shared_gap_NM = [N_actual M_actual];
    elseif any(shared_gap_NM ~= [N_actual M_actual])
        error('jul16_squirmer_passive_sphere:cachedGridMismatch', ...
            ['The cached one-body preconditioner has actual [N M] = ', ...
             '[%d %d], but the gap grid has [%d %d].'], ...
            shared_gap_NM(1),shared_gap_NM(2),N_actual,M_actual);
    end
    opt.precond = shared_gap_precond;

    results.gap_sweep.actual.N(ih) = N_actual;
    results.gap_sweep.actual.M(ih) = M_actual;
    results.gap_sweep.exact.velocity(:,:,ih) = exact_velocity_all;
    results.gap_sweep.exact.mode_map(:,:,ih) = mode_map;
    results.gap_sweep.exact.info{ih} = series_info;
    results.gap_sweep.exact.tolerance_stability(ih) = ...
        tolerance_stability;

    for ib = 1:n_beta
        beta = beta_vec(ib);
        B2 = beta*B1;
        exact_velocity = exact_velocity_all(:,ib);

        case_result = solve_squirmer_case( ...
            q,rvec_in,rvec_out,nout,wout,Fvec,B1,B2,U0, ...
            exact_velocity,opt,n_check,a);

        results.gap_sweep.U(:,ib,ih) = case_result.U;
        results.gap_sweep.velocity(:,ib,ih) = ...
            [case_result.V_sq; case_result.V_ps];
        results.gap_sweep.error.squirmer(ib,ih) = ...
            case_result.e_sq;
        results.gap_sweep.error.passive(ib,ih) = ...
            case_result.e_ps;
        results.gap_sweep.error.symmetry(ib,ih) = ...
            case_result.e_sym;
        results.gap_sweep.error.surface(ib,ih) = ...
            case_result.e_surface;
        results.gap_sweep.iters(ib,ih) = case_result.iters;
        results.gap_sweep.lambda_total_norm(ib,ih) = ...
            case_result.lambda_total_norm;
        results.gap_sweep.solve_time(ib,ih) = ...
            case_result.solve_time;

        fprintf(['%6.1f %7.3g %6d %6d %7d %11.4e %11.4e ', ...
            '%11.4e %11.4e %10.3e %10.3e %10.3e %10.3e ', ...
            '%11.3e %8.2f %8.2f\n'], ...
            beta,h,N_actual,M_actual,case_result.iters, ...
            case_result.V_sq,exact_velocity(1), ...
            case_result.V_ps,exact_velocity(2), ...
            case_result.e_sq,case_result.e_ps,case_result.e_sym, ...
            case_result.e_surface,case_result.lambda_total_norm, ...
            case_result.solve_time, ...
            results.gap_sweep.precond_build_time);
    end
end

results.validation.exact.max_gap_reciprocity_defect = max( ...
    cellfun(@(x) x.reciprocity_defect,results.gap_sweep.exact.info));
results.validation.exact.max_gap_tolerance_stability = max( ...
    results.gap_sweep.exact.tolerance_stability);

%% Figure 1: N-convergence errors
figure('Color','w','Name','Squirmer--passive sphere N convergence');
tiledlayout(1,n_beta,'TileSpacing','compact','Padding','compact');
for ib = 1:n_beta
    ax = nexttile();
    N_plot = results.convergence.actual.N;
    loglog(ax,N_plot,results.convergence.error.squirmer(ib,:),'o-', ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName','$e_{\mathrm{sq}}$');
    hold(ax,'on');
    loglog(ax,N_plot,results.convergence.error.passive(ib,:),'s-', ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName','$e_{\mathrm{ps}}$');
    loglog(ax,N_plot,results.convergence.error.symmetry(ib,:),'d-', ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName','$e_{\mathrm{sym}}$');
    loglog(ax,N_plot,results.convergence.error.surface(ib,:),'^-', ...
        'LineWidth',1.4,'MarkerSize',6, ...
        'DisplayName','$e_{\partial\Omega}$');
    format_axes(ax);
    xlabel(ax,'$N$','Interpreter','latex');
    ylabel(ax,'$\mathrm{error}$','Interpreter','latex');
    title(ax,sprintf('$\\beta=%g$',beta_vec(ib)),'Interpreter','latex');
    legend(ax,'Location','best','Interpreter','latex');
end
sgtitle(sprintf( ...
    '$\\mathrm{Unit\\ squirmer--passive\\ sphere},\\ h=%g$', ...
    h_convergence),'Interpreter','latex');

%% Figure 2: exact and MFS velocities versus gap
[h_plot,h_order] = sort(gap_vec);
figure('Color','w','Name','Squirmer--passive sphere velocities');
tiledlayout(1,n_beta,'TileSpacing','compact','Padding','compact');
for ib = 1:n_beta
    ax = nexttile();
    exact_sq = squeeze(results.gap_sweep.exact.velocity(1,ib,h_order));
    exact_ps = squeeze(results.gap_sweep.exact.velocity(2,ib,h_order));
    mfs_sq = squeeze(results.gap_sweep.velocity(1,ib,h_order));
    mfs_ps = squeeze(results.gap_sweep.velocity(2,ib,h_order));

    semilogx(ax,h_plot,exact_sq,'k-', ...
        'LineWidth',1.7,'DisplayName','$V_{\mathrm{sq}}^{\mathrm{exact}}$');
    hold(ax,'on');
    semilogx(ax,h_plot,mfs_sq,'ko', ...
        'LineWidth',1.2,'MarkerSize',6, ...
        'DisplayName','$V_{\mathrm{sq}}^{\mathrm{MFS}}$');
    semilogx(ax,h_plot,exact_ps,'b-', ...
        'LineWidth',1.7,'DisplayName','$V_{\mathrm{ps}}^{\mathrm{exact}}$');
    semilogx(ax,h_plot,mfs_ps,'bs', ...
        'LineWidth',1.2,'MarkerSize',6, ...
        'DisplayName','$V_{\mathrm{ps}}^{\mathrm{MFS}}$');
    format_axes(ax);
    xlabel(ax,'$h$','Interpreter','latex');
    ylabel(ax,'$\mathrm{axial\ velocity}$','Interpreter','latex');
    title(ax,sprintf('$\\beta=%g$',beta_vec(ib)),'Interpreter','latex');
    legend(ax,'Location','best','Interpreter','latex');
end
sgtitle(sprintf( ...
    '$\\mathrm{Gap\\ sweep},\\ N_{\\mathrm{requested}}=%d$', ...
    N_request_gap),'Interpreter','latex');

%% Figure 3: errors versus gap
figure('Color','w','Name','Squirmer--passive sphere gap errors');
tiledlayout(1,n_beta,'TileSpacing','compact','Padding','compact');
for ib = 1:n_beta
    ax = nexttile();
    loglog(ax,h_plot,results.gap_sweep.error.squirmer(ib,h_order),'o-', ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName','$e_{\mathrm{sq}}$');
    hold(ax,'on');
    loglog(ax,h_plot,results.gap_sweep.error.passive(ib,h_order),'s-', ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName','$e_{\mathrm{ps}}$');
    loglog(ax,h_plot,results.gap_sweep.error.symmetry(ib,h_order),'d-', ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName','$e_{\mathrm{sym}}$');
    loglog(ax,h_plot,results.gap_sweep.error.surface(ib,h_order),'^-', ...
        'LineWidth',1.4,'MarkerSize',6, ...
        'DisplayName','$e_{\partial\Omega}$');
    format_axes(ax);
    xlabel(ax,'$h$','Interpreter','latex');
    ylabel(ax,'$\mathrm{error}$','Interpreter','latex');
    title(ax,sprintf('$\\beta=%g$',beta_vec(ib)),'Interpreter','latex');
    legend(ax,'Location','best','Interpreter','latex');
end
sgtitle('$\mathrm{Gap-dependent\ errors}$','Interpreter','latex');

%% Local functions

function repo_root = find_repo_root(start_dir)
repo_root = start_dir;
while true
    if isfile(fullfile(repo_root,'setpath.m'))
        return;
    end
    parent_dir = fileparts(repo_root);
    if isempty(parent_dir) || strcmp(parent_dir,repo_root)
        error('jul16_squirmer_passive_sphere:repoRootNotFound', ...
            'Could not find setpath.m starting from %s.',start_dir);
    end
    repo_root = parent_dir;
end
end

function q = sphere_centres(h,a)
q = [0 0 -(a+h/2); 0 0 (a+h/2)];
end

function [N_actual,M_actual] = actual_design_sizes(N_request,a_glob)
[rin,~] = get_sphdesign(N_request);
[rout,~] = get_sphdesign(ceil(a_glob*N_request));
N_actual = size(rin,1);
M_actual = size(rout,1);
end

function [nout,wout] = sphere_normals_and_weights(q,rvec_out,a)
P = size(q,1);
M = size(rvec_out,1)/P;
nout = (rvec_out-repelem(q,M,1))/a;
nout = nout./vecnorm(nout,2,2);
wout = (4*pi*a^2/M)*ones(P*M,1);
end

function opt = configure_options( ...
        opt,gmres_tol,maxit,use_fmm,fmm_tol, ...
        inner_only,symmetrize_weighted,add_rank1,outer_force)
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

function precond = build_complete_preconditioner( ...
        rin_body,rout_body,nout_body,wout_body,opt)
q_body = [0 0 0];
nin_body = [];
if opt.add_rank1
    nin_body = get_inner_normals(rin_body,q_body,opt,[]);
end

[Y,UU,LL,Kin,~,Fmap] = oneBodyPrecondMobDLP( ...
    rin_body,rout_body,opt,q_body,wout_body,nout_body,nin_body);
Sblock = stokes_SLP_mat(rin_body,rout_body);
Tself = stokes_DLP_mat(rout_body,rin_body,nout_body);
Wself = diag(repelem(wout_body(:),3));
TWSblock = Tself*Wself*Sblock;
if opt.symmetrize_weighted
    TWSblock = 0.5*(TWSblock+TWSblock.');
end
if opt.add_rank1
    TWSblock = TWSblock + ...
        normal_rank1_block(nin_body,opt.rank1_scale);
end

precond = struct( ...
    'Y',Y, ...
    'UU',UU, ...
    'LL',LL, ...
    'Kin',Kin, ...
    'Fmap',Fmap, ...
    'Sblock',Sblock, ...
    'TWSblock',TWSblock);
end

function result = solve_squirmer_case( ...
        q,rvec_in,rvec_out,nout,wout,Fvec,B1,B2,U0, ...
        exact_velocity,opt,n_check,a)
P = 2;
M = size(rvec_out,1)/P;
u_slip = zeros(P*M,3);
u_slip(1:M,:) = evaluate_squirmer_slip(nout(1:M,:),B1,B2);

tangential_defect = max(abs(sum(nout.*u_slip,2)));
if tangential_defect > 1e-12
    error('jul16_squirmer_passive_sphere:tangentiality', ...
        'Slip tangentiality defect is %.3e.',tangential_defect);
end

solve_start = tic;
[U,iters,~,~,lambda] = solve_mobility_with_DLP( ...
    q,rvec_in,rvec_out,nout,wout,Fvec,u_slip,opt);
solve_time = toc(solve_start);

V_sq = U(3);
V_ps = U(9);
e_sq = abs(V_sq-exact_velocity(1))/U0;
e_ps = abs(V_ps-exact_velocity(2))/U0;

V1_transverse = U(1:2);
V2_transverse = U(7:8);
Omega1 = U(4:6);
Omega2 = U(10:12);
e_sym = sqrt( ...
    norm(V1_transverse,2)^2 + norm(V2_transverse,2)^2 + ...
    a^2*(norm(Omega1,2)^2+norm(Omega2,2)^2))/U0;

rcheck1 = q(1,:) + a*n_check;
rcheck2 = q(2,:) + a*n_check;
rcheck = [rcheck1; rcheck2];
u_slip_check = zeros(size(rcheck));
u_slip_check(1:size(n_check,1),:) = ...
    evaluate_squirmer_slip(n_check,B1,B2);
u_slip_check_vec = reshape(u_slip_check.',[],1);

u_mfs_check = getStokesletFlow(lambda,rvec_in,rcheck,opt);
u_rigid_check = [ ...
    getKmat(rcheck1,q(1,:))*U(1:6); ...
    getKmat(rcheck2,q(2,:))*U(7:12)];
boundary_residual = ...
    u_mfs_check-u_rigid_check-u_slip_check_vec;
slip_norm = norm(u_slip_check_vec,inf);
e_surface = norm(boundary_residual,inf)/slip_norm;

result = struct( ...
    'U',U, ...
    'V_sq',V_sq, ...
    'V_ps',V_ps, ...
    'e_sq',e_sq, ...
    'e_ps',e_ps, ...
    'e_sym',e_sym, ...
    'e_surface',e_surface, ...
    'iters',iters, ...
    'lambda_total_norm',norm(lambda,inf), ...
    'solve_time',solve_time, ...
    'tangential_defect',tangential_defect);
end

function u_slip = evaluate_squirmer_slip(n,B1,B2)
ez = [0 0 1];
nz = n(:,3);
tangent_ez = repmat(ez,size(n,1),1)-nz.*n;
u_slip = -(B1+B2*nz).*tangent_ez;
end

function [mode_map,info,tolerance_stability] = ...
        validated_exact_map(h,a,series_opts)
[mode_map,info] = papavassiliou_squirmer_passive_map( ...
    h,a,series_opts);
if ~info.converged
    error('jul16_squirmer_passive_sphere:exactConvergence', ...
        'Exact series did not converge at h = %.17g.',h);
end
if info.reciprocity_defect > 1e-10
    error('jul16_squirmer_passive_sphere:reciprocity', ...
        'Resistance reciprocity defect %.3e at h = %.17g.', ...
        info.reciprocity_defect,h);
end

tighter_opts = series_opts;
tighter_opts.tol = max(series_opts.tol/10,1e-15);
[mode_map_tight,info_tight] = ...
    papavassiliou_squirmer_passive_map(h,a,tighter_opts);
tolerance_stability = norm(mode_map-mode_map_tight,inf) / ...
    max(norm(mode_map_tight,inf),realmin);
if ~info_tight.converged || tolerance_stability > 1e-9
    error('jul16_squirmer_passive_sphere:toleranceStability', ...
        ['Exact series changed by %.3e under a tighter tolerance at ', ...
         'h = %.17g.'],tolerance_stability,h);
end
info.tighter_tolerance = info_tight;
info.tolerance_stability = tolerance_stability;
end

function R = axis_angle_rotation(axis,angle)
axis = axis(:)/norm(axis);
K = [0 -axis(3) axis(2); ...
     axis(3) 0 -axis(1); ...
     -axis(2) axis(1) 0];
R = eye(3)+sin(angle)*K+(1-cos(angle))*(K*K);
end

function format_axes(ax)
grid(ax,'on');
box(ax,'on');
set(ax,'TickLabelInterpreter','latex','FontSize',11);
end
