%JUL16_FIG2P7D_TWO_SPHERE_MOBILITY_PROXY_SWEEP
%
% Mobility analogue of Fig. 2.7(c,d) in Broms's PhD thesis:
% two unit spheres in squeezing motion, with a sweep over proxy radius.
% The comparison is between
%   1. solve_mobility
%   2. solve_mobility_with_DLP
%
% The original figure prescribes equal and opposite sphere velocities. For
% the mobility problem, equal and opposite unit forces are prescribed along
% the line of centres instead. Each solver's off-grid surface residual is
% normalized by the maximum magnitude of its computed no-slip velocity.
%
% Full experiment:
%   run('experiments/dev_DLP/jul16_two_sphere_mobility_proxy_sweep.m')
%
%
% Created July 16, 2026.

clear;
close all;
clc;


%% Paths
script_dir = fileparts(mfilename('fullpath'));
repo_root = find_repo_root(script_dir);
run(fullfile(repo_root,'setpath.m'));

%% mobility configuration
delta_vec = [0.2 2];
N_request_vec = [482 686 969];
Rp_vec = 0.10:0.05:0.95;
a_glob = 1.2;
P = 2;

% Squeezing forces: sphere 1 moves toward sphere 2 and vice versa.
Fvec = zeros(6*P,1);
Fvec(1) = 1;
Fvec(7) = -1;

gmres_tol = 1e-13;
maxit = 600;
use_fmm = false;
fmm_tol = 1e-12;

% DLP/SLP solver choices. outer_force=true uses the force/torque map of
% the complete T*W*S representation.
dlp_inner_only = false;
dlp_symmetrize_weighted = false;
dlp_add_rank1 = false;
dlp_outer_force = true;


method_names = {'solve_mobility','solve_mobility_with_DLP'};
method_labels = {'standard MFS','DLP/SLP MFS'};
remove_titles = false; % Set true for title-free figures.

cfg = struct();
cfg.delta_vec = delta_vec;
cfg.N_request_vec = N_request_vec;
cfg.Rp_vec = Rp_vec;
cfg.a_glob = a_glob;
cfg.P = P;
cfg.Fvec = Fvec;
cfg.gmres_tol = gmres_tol;
cfg.maxit = maxit;
cfg.use_fmm = use_fmm;
cfg.fmm_tol = fmm_tol;
cfg.dlp_inner_only = dlp_inner_only;
cfg.dlp_symmetrize_weighted = dlp_symmetrize_weighted;
cfg.dlp_add_rank1 = dlp_add_rank1;
cfg.dlp_outer_force = dlp_outer_force;
cfg.method_names = method_names;

results_file = fullfile(repo_root,'data', ...
    'jul16_two_sphere_mobility_proxy_sweep.mat');
save_results = 0;
resume_results = 1;

%% Storage
n_method = numel(method_names);
n_delta = numel(delta_vec);
n_N = numel(N_request_vec);
n_Rp = numel(Rp_vec);

if resume_results && isfile(results_file)
    saved = load(results_file,'results');
    if ~isfield(saved,'results') || ~isfield(saved.results,'config') || ...
            ~isequaln(saved.results.config,cfg)
        error('jul16_fig2p7d:incompatibleResults', ...
            ['Existing result file has a different configuration:\n  %s\n', ...
             'Move or delete it before starting this sweep.'],results_file);
    end
    results = saved.results;
    fprintf('Resuming compatible results from\n  %s\n\n',results_file);
else
    results = struct();
    results.created = char(datetime('now', ...
        'Format','yyyyMMdd''T''HHmmss'));
    results.description = [ ...
        'Mobility analogue of Fig. 2.7(c,d): two squeezing spheres, ', ...
        'standard MFS versus DLP/SLP MFS.'];
    results.config = cfg;
    results.actual.N = nan(1,n_N);
    results.actual.M = nan(1,n_N);
    results.completed = false(n_method,n_delta,n_N,n_Rp);
    results.U = nan(6*P,n_method,n_delta,n_N,n_Rp);
    results.iters = nan(n_method,n_delta,n_N,n_Rp);
    results.uerr = nan(n_method,n_delta,n_N,n_Rp);
    results.lambda_norm = nan(n_method,n_delta,n_N,n_Rp);
    results.solve_time = nan(n_method,n_delta,n_N,n_Rp);
    results.rel_U_dlp_vs_standard = nan(n_delta,n_N,n_Rp);
end

fprintf('Fig. 2.7(c,d) two-sphere mobility proxy-radius sweep\n');
fprintf('  delta = [%s]\n',num2str(delta_vec));
fprintf('  N request = [%s]\n',num2str(N_request_vec));
fprintf('  Rp = %.2f:%.2f:%.2f\n', ...
    Rp_vec(1),Rp_step(Rp_vec),Rp_vec(end));
fprintf('  M/N target = %.3g\n',a_glob);
fprintf('  squeezing forces: F1 = [+1,0,0], F2 = [-1,0,0]\n');
fprintf('  gmres_tol = %.1e, maxit = %d\n',gmres_tol,maxit);
fprintf('  DLP options: inner_only=%d, symmetrize_weighted=%d, ', ...
    dlp_inner_only,dlp_symmetrize_weighted);
fprintf('add_rank1=%d, outer_force=%d\n\n', ...
    dlp_add_rank1,dlp_outer_force);

fprintf('%7s %5s %5s %6s | %6s %10s %10s %8s | %6s %10s %10s %8s | %10s\n', ...
    'delta','Rp','N','M', ...
    'it std','uerr std','|lam| std','t std', ...
    'it DLP','uerr DLP','|lam| DLP','t DLP','rel U');

%% Main sweep
for id = 1:n_delta
    delta = delta_vec(id);
    centre_offset = 1 + delta/2;
    q = [-centre_offset 0 0; centre_offset 0 0];

    for in = 1:n_N
        N_request = N_request_vec(in);

        for ir = 1:n_Rp
            Rp = Rp_vec(ir);

            [rvec_in,rvec_out,opt] = init_spheres( ...
                q,Rp,N_request,a_glob);
            opt = configure_options( ...
                opt,gmres_tol,maxit,use_fmm,fmm_tol);

            N_actual = size(rvec_in,1)/P;
            M_actual = size(rvec_out,1)/P;
            [nout,wout] = sphere_normals_and_weights(q,rvec_out);

            results.actual.N(in) = N_actual;
            results.actual.M(in) = M_actual;

            if ~results.completed(1,id,in,ir)
                t_start = tic;
                [U_std,it_std,ln_std,uerr_std] = solve_mobility( ...
                    q,rvec_in,rvec_out,Fvec,opt);
                t_std = toc(t_start);

                results.U(:,1,id,in,ir) = U_std;
                results.iters(1,id,in,ir) = it_std;
                results.uerr(1,id,in,ir) = uerr_std;
                results.lambda_norm(1,id,in,ir) = ln_std;
                results.solve_time(1,id,in,ir) = t_std;
                results.completed(1,id,in,ir) = true;
            end

            if ~results.completed(2,id,in,ir)
                opt_dlp = configure_dlp_options(opt, ...
                    dlp_inner_only,dlp_symmetrize_weighted, ...
                    dlp_add_rank1,dlp_outer_force);

                t_start = tic;
                [U_dlp,it_dlp,ln_dlp,uerr_dlp] = ...
                    solve_mobility_with_DLP( ...
                    q,rvec_in,rvec_out,nout,wout,Fvec,[],opt_dlp);
                t_dlp = toc(t_start);

                results.U(:,2,id,in,ir) = U_dlp;
                results.iters(2,id,in,ir) = it_dlp;
                results.uerr(2,id,in,ir) = uerr_dlp;
                results.lambda_norm(2,id,in,ir) = ln_dlp;
                results.solve_time(2,id,in,ir) = t_dlp;
                results.completed(2,id,in,ir) = true;
            end

            U_std = results.U(:,1,id,in,ir);
            U_dlp = results.U(:,2,id,in,ir);
            results.rel_U_dlp_vs_standard(id,in,ir) = ...
                relerr_inf(U_dlp,U_std);

            fprintf('%7.3g %5.2f %5d %6d | %6d %10.3e %10.3e %8.2f | %6d %10.3e %10.3e %8.2f | %10.3e\n', ...
                delta,Rp,N_actual,M_actual, ...
                results.iters(1,id,in,ir), ...
                results.uerr(1,id,in,ir), ...
                results.lambda_norm(1,id,in,ir), ...
                results.solve_time(1,id,in,ir), ...
                results.iters(2,id,in,ir), ...
                results.uerr(2,id,in,ir), ...
                results.lambda_norm(2,id,in,ir), ...
                results.solve_time(2,id,in,ir), ...
                results.rel_U_dlp_vs_standard(id,in,ir));

            if save_results
                save(results_file,'results');
            end
        end
    end
end

%% Plots
plot_solver_figures(results,'uerr', ...
    '$e_{\mathrm{res}}^{\max}$', ...
    'Fig. 2.7(d) mobility analogue: proxy-radius residual', ...
    method_labels,remove_titles);

plot_solver_figures(results,'lambda_norm', ...
    '$\|\lambda\|_\infty$', ...
    'Fig. 2.7(c) mobility analogue: coefficient magnitude', ...
    method_labels,remove_titles);

plot_velocity_difference(results,remove_titles);

%% Local functions

function repo_root = find_repo_root(start_dir)
repo_root = start_dir;
while true
    if isfile(fullfile(repo_root,'setpath.m'))
        return;
    end

    parent_dir = fileparts(repo_root);
    if isempty(parent_dir) || strcmp(parent_dir,repo_root)
        error('jul16_fig2p7d:repoRootNotFound', ...
            'Could not find setpath.m starting from %s.',start_dir);
    end
    repo_root = parent_dir;
end
end

function opt = configure_options(opt,gmres_tol,maxit,use_fmm,fmm_tol)
opt.fmm = use_fmm;
opt.fmm_tol = fmm_tol;
opt.lr = 0;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = 0;
opt.maxit = maxit;
opt.plot = 0;
opt.debug = 0;
opt.profile = 0;
opt.compute_residual = true;
end

function opt = configure_dlp_options( ...
        opt,inner_only,symmetrize_weighted,add_rank1,outer_force)
opt.inner_only = inner_only;
opt.symmetrize_weighted = symmetrize_weighted;
opt.add_rank1 = add_rank1;
opt.outer_force = outer_force;
end

function [nout,wout] = sphere_normals_and_weights(q,rvec_out)
P = size(q,1);
M = size(rvec_out,1)/P;
nout = rvec_out - kron(q,ones(M,1));
nout = nout ./ vecnorm(nout,2,2);
wout = (4*pi/M) * ones(P*M,1);
end

function e = relerr_inf(a,b)
den = norm(b,inf);
if den == 0
    e = norm(a-b,inf);
else
    e = norm(a-b,inf)/den;
end
end

function step = Rp_step(Rp_vec)
if numel(Rp_vec) > 1
    step = Rp_vec(2) - Rp_vec(1);
else
    step = 0;
end
end

function plot_solver_figures( ...
        results,metric_field,y_label,figure_name,method_labels,remove_titles)
cfg = results.config;

for im = 1:numel(cfg.method_names)
    window_name = sprintf('%s - %s',figure_name,method_labels{im});
    figure('Color','w','Name',window_name);
    ax = axes();
    plot_fig2p7_panel(ax,results,metric_field,im);
    xlabel(ax,'$R_p$','Interpreter','latex');
    ylabel(ax,y_label,'Interpreter','latex');
    if ~remove_titles
        title(ax,sprintf('%s: %s',figure_name,method_labels{im}), ...
            'Interpreter','latex');
    end
end
end

function plot_fig2p7_panel(ax,results,metric_field,method_index)
cfg = results.config;
delta_colors = [0.00 0.45 0.74; 1.00 0.10 0.10];
line_styles = {'--','-',':'};
metric = results.(metric_field);

hold(ax,'on');
for id = 1:numel(cfg.delta_vec)
    color = delta_colors(min(id,size(delta_colors,1)),:);
    for in = 1:numel(cfg.N_request_vec)
        line_style = line_styles{min(in,numel(line_styles))};
        y = squeeze(metric(method_index,id,in,:));
        plot(ax,cfg.Rp_vec,y, ...
            'Color',color, ...
            'LineStyle',line_style, ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('$N=%d,\\ \\delta=%.3g$', ...
            results.actual.N(in),cfg.delta_vec(id)));
    end
end

set(ax,'YScale','log');
set(ax,'TickLabelInterpreter','latex');
grid(ax,'on');
box(ax,'on');
if numel(cfg.Rp_vec) > 1
    xlim(ax,[min(cfg.Rp_vec),max(cfg.Rp_vec)]);
    xticks(ax,0.2:0.2:0.8);
end
legend(ax,'Location','best','Interpreter','latex');
end

function plot_velocity_difference(results,remove_titles)
cfg = results.config;
delta_colors = [0.00 0.45 0.74; 1.00 0.10 0.10];
line_styles = {'--','-',':'};

figure('Color','w','Name','DLP versus standard rigid-body velocity');
ax = axes();
hold(ax,'on');

for id = 1:numel(cfg.delta_vec)
    color = delta_colors(min(id,size(delta_colors,1)),:);
    for in = 1:numel(cfg.N_request_vec)
        line_style = line_styles{min(in,numel(line_styles))};
        y = squeeze(results.rel_U_dlp_vs_standard(id,in,:));
        plot(ax,cfg.Rp_vec,y, ...
            'Color',color, ...
            'LineStyle',line_style, ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('$N=%d,\\ \\delta=%.3g$', ...
            results.actual.N(in),cfg.delta_vec(id)));
    end
end

set(ax,'YScale','log');
set(ax,'TickLabelInterpreter','latex');
grid(ax,'on');
box(ax,'on');
xlabel(ax,'$R_p$','Interpreter','latex');
ylabel(ax, ...
    ['$\|\mathbf{U}_{\mathrm{DLP}}-\mathbf{U}_{\mathrm{standard}}\|_\infty/' ...
     '\|\mathbf{U}_{\mathrm{standard}}\|_\infty$'], ...
    'Interpreter','latex');
if ~remove_titles
    title(ax,'DLP/SLP MFS versus standard MFS', ...
        'Interpreter','latex');
end
if numel(cfg.Rp_vec) > 1
    xlim(ax,[min(cfg.Rp_vec),max(cfg.Rp_vec)]);
    xticks(ax,0.2:0.2:0.8);
end
legend(ax,'Location','best','Interpreter','latex');
end
