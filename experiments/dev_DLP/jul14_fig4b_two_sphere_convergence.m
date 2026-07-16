%JUL14_FIG4B_TWO_SPHERE_CONVERGENCE  
%
% This script sets up the two-sphere convergence test from Fig. 4b of
% Broms, Barnett, and Tornberg (J. Comput. Phys., 2025), but for the
% mobility problem only. The standard MFS mobility solver is compared with
% solve_mobility_with_DLP.
%
% Differences from the paper's Fig. 4b:
%   1. Forces/torques are prescribed directly; no resistance solve is used.
%   2. The high-resolution N = 3529 mobility solution is used as Uref.
%   3. The dashed R_acc convergence guide lines are intentionally omitted.
%
% Created July 14, 2026.

clear;
close all;
clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
repo_root = find_repo_root(script_dir);
run(fullfile(repo_root,'setpath.m'));

%% Fig. 4b mobility configuration
delta_vec = [0.1 0.2];
N_request_vec = 100:200:1900;
N_ref_request = 3529;
%N_ref_request = 100;

Rp = 0.6;
a_glob = 1.2;
P = 2;

% Mobility analogue of the translation boundary condition in Fig. 4b:
% both particles receive the same force and no torque.
Fvec = zeros(12,1);
Fvec(1) = 1; 
Fvec(7) = 1; 

gmres_tol = 1e-10;
maxit = 600;
use_fmm = false;

dlp_inner_only = false;
dlp_symmetrize_weighted = false;
dlp_add_rank1 = false;
dlp_outer_force = true; %change here!

method_names = {'solve_mobility','solve_mobility_with_DLP'};

%% Storage
n_method = numel(method_names);
n_delta = numel(delta_vec);
n_N = numel(N_request_vec);

results = struct();
results.created = '2026-07-14';
results.description = ['Two-sphere Fig. 4b mobility analogue: ', ...
    'solve_mobility vs solve_mobility_with_DLP.'];
results.delta_vec = delta_vec;
results.N_request_vec = N_request_vec;
results.N_ref_request = N_ref_request;
results.Rp = Rp;
results.a_glob = a_glob;
results.P = P;
results.Fvec = Fvec;
results.method_names = method_names;
results.reference_method = 'solve_mobility';
results.gmres_tol = gmres_tol;
results.maxit = maxit;
results.use_fmm = use_fmm;
results.fmm_tol = fmm_tol;
results.dlp_inner_only = dlp_inner_only;
results.dlp_symmetrize_weighted = dlp_symmetrize_weighted;
results.dlp_add_rank1 = dlp_add_rank1;
results.dlp_outer_force = dlp_outer_force;
results.compute_Racc_reference_lines = false;

results.reference.N_actual = nan(1,n_delta);
results.reference.M_actual = nan(1,n_delta);
results.reference.U = nan(6*P,n_delta);
results.reference.iters = nan(1,n_delta);
results.reference.uerr = nan(1,n_delta);
results.reference.lambda_norm = nan(1,n_delta);
results.reference.solve_time = nan(1,n_delta);

results.actual.N = nan(n_delta,n_N);
results.actual.M = nan(n_delta,n_N);

results.mobility.iters = nan(n_method,n_delta,n_N);
results.mobility.uerr = nan(n_method,n_delta,n_N);
results.mobility.Uerr = nan(n_method,n_delta,n_N);
results.mobility.lambda_norm = nan(n_method,n_delta,n_N);
results.mobility.solve_time = nan(n_method,n_delta,n_N);
results.mobility.coeff_count = nan(n_method,n_delta,n_N);
results.mobility.gmres_unknown_count = nan(n_method,n_delta,n_N);

fprintf('Jul 14 Fig. 4b two-sphere mobility comparison\n');
fprintf('  delta = [%s]\n',num2str(delta_vec));
fprintf('  N request = [%s]\n',num2str(N_request_vec));
fprintf('  N reference request = %d\n',N_ref_request);
fprintf('  Rp = %.3g, a_glob = %.3g\n',Rp,a_glob);
fprintf('  gmres_tol = %.1e, maxit = %d\n\n',gmres_tol,maxit);

fprintf('%8s %6s %6s | %8s %10s %10s %10s | %8s %10s %10s %10s\n', ...
    'delta','N','M', ...
    'it std','uerr std','Uerr std','|lam| std', ...
    'it DLP','uerr DLP','Uerr DLP','|lam| DLP');

%% Main sweep
for id = 1:n_delta
    delta = delta_vec(id);
    q = [0 0 0; 2 + delta 0 0];

    [rref_in,rref_out,opt_ref] = init_spheres(q,Rp,N_ref_request,a_glob);
    opt_ref = configure_options(opt_ref,gmres_tol,maxit,use_fmm,fmm_tol);

    t_start = tic;
    [Uref,it_ref,ln_ref,uerr_ref] = solve_mobility(q,rref_in,rref_out,Fvec,opt_ref);
    t_ref = toc(t_start);

    results.reference.N_actual(id) = size(rref_in,1)/P;
    results.reference.M_actual(id) = size(rref_out,1)/P;
    results.reference.U(:,id) = Uref;
    results.reference.iters(id) = it_ref;
    results.reference.uerr(id) = uerr_ref;
    results.reference.lambda_norm(id) = ln_ref;
    results.reference.solve_time(id) = t_ref;

    for in = 1:n_N
        N_request = N_request_vec(in);

        [rvec_in,rvec_out,opt] = init_spheres(q,Rp,N_request,a_glob);
        opt = configure_options(opt,gmres_tol,maxit,use_fmm,fmm_tol);

        N_actual = size(rvec_in,1)/P;
        M_actual = size(rvec_out,1)/P;
        [nout,wout] = sphere_normals_and_weights(q,rvec_out);

        results.actual.N(id,in) = N_actual;
        results.actual.M(id,in) = M_actual;

        t_start = tic;
        [U_std,it_std,ln_std,uerr_std] = solve_mobility(q,rvec_in,rvec_out,Fvec,opt);
        t_std = toc(t_start);

        results.mobility.iters(1,id,in) = it_std;
        results.mobility.uerr(1,id,in) = uerr_std;
        results.mobility.Uerr(1,id,in) = relerr_inf(U_std,Uref);
        results.mobility.lambda_norm(1,id,in) = ln_std;
        results.mobility.solve_time(1,id,in) = t_std;
        results.mobility.coeff_count(1,id,in) = 3*N_actual*P;
        results.mobility.gmres_unknown_count(1,id,in) = 3*M_actual*P;

        opt_dlp = configure_dlp_options(opt, ...
            dlp_inner_only,dlp_symmetrize_weighted,dlp_add_rank1,dlp_outer_force);

        t_start = tic;
        [U_dlp,it_dlp,ln_dlp,uerr_dlp] = solve_mobility_with_DLP( ...
            q,rvec_in,rvec_out,nout,wout,Fvec,opt_dlp);
        t_dlp = toc(t_start);

        results.mobility.iters(2,id,in) = it_dlp;
        results.mobility.uerr(2,id,in) = uerr_dlp;
        results.mobility.Uerr(2,id,in) = relerr_inf(U_dlp,Uref);
        results.mobility.lambda_norm(2,id,in) = ln_dlp;
        results.mobility.solve_time(2,id,in) = t_dlp;
        results.mobility.coeff_count(2,id,in) = 3*N_actual*P;
        results.mobility.gmres_unknown_count(2,id,in) = 3*N_actual*P;

        fprintf('%8.3g %6d %6d | %8d %10.3e %10.3e %10.3e | %8d %10.3e %10.3e %10.3e\n', ...
            delta,N_actual,M_actual, ...
            it_std,uerr_std,results.mobility.Uerr(1,id,in),ln_std, ...
            it_dlp,uerr_dlp,results.mobility.Uerr(2,id,in),ln_dlp);

        results_file = fullfile(repo_root,'experiments','jul14_fig4b_two_sphere_convergence.mat');
        save(results_file,'results');
    end
end


%% Plot: no R_acc reference lines
%results_file = fullfile(repo_root,'experiments','jul14_fig4b_two_sphere_mobility_dlp_compare_results.mat');
%figure('Name','Fig. 4b mobility analogue: residual and velocity error');
%results_file = fullfile(repo_root,'data','jul14_fig4b_two_sphere_mobility_dlp_compare_results.mat');
results_file = fullfile(repo_root,'data','jul14_fig4b_two_sphere_convergence_07.mat');
load(results_file);

plot_fig4b_metric(results,'uerr', ...
    '$\epsilon_{\mathrm{res}}^{\max}$', ...
    'Mobility residual comparison',true);

plot_fig4b_metric(results,'Uerr', ...
    '$\epsilon_U$', ...
    'Mobility velocity error comparison',true);

plot_fig4b_metric(results,'lambda_norm', ...
    '$\|\lambda\|_\infty$', ...
    'Mobility coefficient norm comparison',true);

plot_fig4b_metric(results,'iters', ...
    '$n_{\mathrm{GMRES}}$', ...
    'Mobility GMRES iteration comparison',false);

%% Local functions

function repo_root = find_repo_root(start_dir)
repo_root = start_dir;
while true
    if isfile(fullfile(repo_root,'setpath.m'))
        return;
    end

    parent_dir = fileparts(repo_root);
    if isempty(parent_dir) || strcmp(parent_dir,repo_root)
        error('jul14_fig4b_two_sphere_convergence:repoRootNotFound', ...
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

function opt = configure_dlp_options(opt,inner_only,symmetrize_weighted,add_rank1,outer_force)
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

function plot_fig4b_metric(results,metric_field,y_label,figure_name,use_log)
if nargin < 5
    use_log = true;
end
figure('Name',figure_name);
plot_fig4b_metric_panel(results,metric_field,y_label,use_log);
title(figure_name,'Interpreter','latex');
end

function plot_fig4b_metric_panel(results,metric_field,y_label,use_log)
delta_vec = results.delta_vec;
N_actual = results.actual.N;
method_names = results.method_names;
metric_data = results.mobility.(metric_field);
colors = [0.85 0.10 0.10; 0.00 0.45 0.74];
line_styles = {'-','--'};
markers = {'o','s'};

hold on;
for id = 1:numel(delta_vec)
    x = sqrt(N_actual(id,:));
    color = colors(min(id,size(colors,1)),:);

    for im = 1:numel(method_names)
        marker = markers{min(im,numel(markers))};
        line_style = line_styles{min(im,numel(line_styles))};
        if im == 1
            marker_face_color = color;
        else
            marker_face_color = 'none';
        end

        plot(x,squeeze(metric_data(im,id,:)), ...
            line_style, ...
            'Color',color, ...
            'Marker',marker, ...
            'MarkerFaceColor',marker_face_color, ...
            'MarkerEdgeColor',color, ...
            'DisplayName',sprintf('$%s,\\ \\delta = %.3g$', ...
                method_latex_name(method_names{im}),delta_vec(id)));
    end
end

grid('on');
box('on');
ax = gca;
set(ax,'TickLabelInterpreter','latex');
if use_log
    set(ax,'YScale','log');
end
xlabel('$N$','Interpreter','latex');
ylabel(y_label,'Interpreter','latex');
set_sqrtN_ticks(ax,results.N_request_vec);
legend(ax,'Location','southwest','Interpreter','latex');
end

function set_sqrtN_ticks(ax,N_request_vec)
tick_N = unique([N_request_vec(1), 200, 500, 1000, 1500, N_request_vec(end)]);
tick_N = tick_N(tick_N >= min(N_request_vec) & tick_N <= max(N_request_vec));
set(ax,'XTick',sqrt(tick_N));
set(ax,'XTickLabel',compose('$%d$',tick_N));
set(ax,'TickLabelInterpreter','latex');
xlim(ax,sqrt([min(N_request_vec), max(N_request_vec)]));
end

function name = method_latex_name(method_name)
switch method_name
    case 'solve_mobility'
        name = '\mathrm{rectangular}';
    case 'solve_mobility_with_DLP'
        name = '\mathrm{symmetrized}';
    otherwise
        name = strrep(method_name,'_','\_');
end
end
