%JUL08_MOBILITY_SYM_LEFT_SWEEP Compare left/right/standard mobility solves.
%
% This script sweeps sphere-cluster separations and particle counts. It is
% intended to be run manually.

clear;
close all;
clc;

%% Mode and sweep configuration
symmetrize_weighted_vec = [false true];
left_inner_only = true;
right_dlp_inner_only = true;
delta_vec = [0.1, 0.5, 1, 5];
P_vec = 5:5:30;
base_seed = 708;

%% Discretization and solver options
Rp = 0.30;
Ndes = 100;
%Ndes = 200;
a = 2;
gmres_tol = 1e-10;
maxit = 500;
use_fmm = 0;

%% Storage
n_delta = numel(delta_vec);
n_P = numel(P_vec);
n_mode = numel(symmetrize_weighted_vec);

results = struct();
results.symmetrize_weighted_vec = symmetrize_weighted_vec;
results.left_inner_only = left_inner_only;
results.right_dlp_inner_only = right_dlp_inner_only;
results.delta_vec = delta_vec;
results.P_vec = P_vec;
results.Rp = Rp;
results.Ndes = Ndes;
results.a = a;
results.gmres_tol = gmres_tol;
results.maxit = maxit;

results.left.iters = nan(n_mode,n_delta,n_P);
results.left.uerr = nan(n_mode,n_delta,n_P);
results.left.lambda_norm = nan(n_mode,n_delta,n_P);
results.left.solve_time = nan(n_mode,n_delta,n_P);

results.right_dlp.iters = nan(n_mode,n_delta,n_P);
results.right_dlp.uerr = nan(n_mode,n_delta,n_P);
results.right_dlp.lambda_norm = nan(n_mode,n_delta,n_P);
results.right_dlp.solve_time = nan(n_mode,n_delta,n_P);

results.standard.iters = nan(n_delta,n_P);
results.standard.uerr = nan(n_delta,n_P);
results.standard.lambda_norm = nan(n_delta,n_P);
results.standard.solve_time = nan(n_delta,n_P);

results.rel_U_left_right_dlp = nan(n_mode,n_delta,n_P);
results.rel_U_left_standard = nan(n_mode,n_delta,n_P);
results.rel_U_right_dlp_standard = nan(n_mode,n_delta,n_P);

fprintf('Jul 08 mobility sym-left sweep\n');
fprintf('  symmetrize_weighted_vec = [%s]\n',num2str(symmetrize_weighted_vec));
fprintf('  left_inner_only = %d\n',left_inner_only);
fprintf('  right_dlp_inner_only = %d\n',right_dlp_inner_only);
fprintf('  delta = [%s]\n',num2str(delta_vec));
fprintf('  P     = [%s]\n\n',num2str(P_vec));
fprintf('%4s %8s %5s | %8s %10s %8s | %8s %10s %8s | %8s %10s %8s | %11s %11s %11s\n', ...
    'sym','delta','P','it left','uerr left','t left', ...
    'it DLP','uerr DLP','t DLP','it std','uerr std','t std', ...
    'rel L/DLP','rel L/std','rel DLP/std');

%% Main sweep
for id = 1:n_delta
    delta = delta_vec(id);
    for ip = 1:n_P
        P = P_vec(ip);
        rng(base_seed + 1000*id + P);

        [q,~] = grow_cluster(P,delta);
        Fref = rand(6*P,1);

        [rvec_in,rvec_out,opt] = init_spheres(q,Rp,Ndes,a);
        opt.fmm = use_fmm;
        opt.lr = 0;
        opt.gmres_tol = gmres_tol;
        opt.gmres_verbose = 0;
        opt.maxit = maxit;
        opt.compute_residual = true;
        opt.debug = 0;

        M = size(rvec_out,1)/P;
        nout = rvec_out - kron(q,ones(M,1));
        nout = nout ./ vecnorm(nout,2,2);
        wout = (4*pi/M) * ones(P*M,1);

        t_start = tic;
        [U_std,it_std,ln_std,uerr_std] = ...
            solve_mobility(q,rvec_in,rvec_out,Fref,opt);
        t_std = toc(t_start);

        results.standard.iters(id,ip) = it_std;
        results.standard.uerr(id,ip) = uerr_std;
        results.standard.lambda_norm(id,ip) = ln_std;
        results.standard.solve_time(id,ip) = t_std;

        for im = 1:n_mode
            symmetrize_weighted = symmetrize_weighted_vec(im);
            opt_mode = opt;
            opt_mode.inner_only = left_inner_only;
            opt_mode.symmetrize_weighted = symmetrize_weighted;

            t_start = tic;
            [U_left,it_left,ln_left,uerr_left] = ...
                solve_mobility_sym_left(q,rvec_in,rvec_out,nout,wout,Fref,opt_mode);
            t_left = toc(t_start);

            opt_dlp = opt_mode;
            opt_dlp.inner_only = right_dlp_inner_only;
            t_start = tic;
            [U_dlp,it_dlp,ln_dlp,uerr_dlp] = ...
                solve_mobility_with_DLP( ...
                    q,rvec_in,rvec_out,nout,wout,Fref,[],opt_dlp);
            t_dlp = toc(t_start);

            results.left.iters(im,id,ip) = it_left;
            results.left.uerr(im,id,ip) = uerr_left;
            results.left.lambda_norm(im,id,ip) = ln_left;
            results.left.solve_time(im,id,ip) = t_left;

            results.right_dlp.iters(im,id,ip) = it_dlp;
            results.right_dlp.uerr(im,id,ip) = uerr_dlp;
            results.right_dlp.lambda_norm(im,id,ip) = ln_dlp;
            results.right_dlp.solve_time(im,id,ip) = t_dlp;

            results.rel_U_left_right_dlp(im,id,ip) = relerr_inf(U_left,U_dlp);
            results.rel_U_left_standard(im,id,ip) = relerr_inf(U_left,U_std);
            results.rel_U_right_dlp_standard(im,id,ip) = relerr_inf(U_dlp,U_std);

            fprintf('%4d %8.3g %5d | %8d %10.3e %8.3g | %8d %10.3e %8.3g | %8d %10.3e %8.3g | %11.3e %11.3e %11.3e\n', ...
                symmetrize_weighted,delta,P,it_left,uerr_left,t_left, ...
                it_dlp,uerr_dlp,t_dlp,it_std,uerr_std,t_std, ...
                results.rel_U_left_right_dlp(im,id,ip), ...
                results.rel_U_left_standard(im,id,ip), ...
                results.rel_U_right_dlp_standard(im,id,ip));
        end
    end
end

save('experiments/jul08_mobility_sym_left_sweep_results.mat','results');

%% Plots
for im = 1:n_mode
    symmetrize_weighted = symmetrize_weighted_vec(im);
    suffix = sprintf('symmetrize\\_weighted = %d',symmetrize_weighted);

    plot_three_solver_metric(P_vec,delta_vec, ...
        squeeze(results.left.iters(im,:,:)), ...
        squeeze(results.right_dlp.iters(im,:,:)), ...
        results.standard.iters, ...
        'GMRES iterations',sprintf('mobility iterations (%s)',suffix),false);

    plot_three_solver_metric(P_vec,delta_vec, ...
        squeeze(results.left.uerr(im,:,:)), ...
        squeeze(results.right_dlp.uerr(im,:,:)), ...
        results.standard.uerr, ...
        'surface residual',sprintf('mobility surface residuals (%s)',suffix),true);

    plot_three_solver_metric(P_vec,delta_vec, ...
        squeeze(results.left.solve_time(im,:,:)), ...
        squeeze(results.right_dlp.solve_time(im,:,:)), ...
        results.standard.solve_time, ...
        'solve time (s)',sprintf('mobility solve times (%s)',suffix),false);

    figure('Color','w','Name',sprintf('mobility velocity differences (%s)',suffix));
    tiledlayout(1,n_delta,'TileSpacing','compact','Padding','compact');
    for id = 1:n_delta
        ax = nexttile;
        hold(ax,'on');
        semilogy(ax,P_vec,squeeze(results.rel_U_left_right_dlp(im,id,:)),'o-','DisplayName','left / right DLP');
        semilogy(ax,P_vec,squeeze(results.rel_U_left_standard(im,id,:)),'s-','DisplayName','left / standard');
        semilogy(ax,P_vec,squeeze(results.rel_U_right_dlp_standard(im,id,:)),'^-','DisplayName','right DLP / standard');
        grid(ax,'on');
        box(ax,'on');
        title(ax,sprintf('\\delta = %.3g',delta_vec(id)));
        xlabel(ax,'P');
        ylabel(ax,'relative U difference');
        if id == 1
            legend(ax,'Location','best');
        end
    end
    sgtitle(sprintf('Jul 08 mobility velocity differences (%s)',suffix),'Interpreter','none');
end

function e = relerr_inf(a,b)
den = norm(b,inf);
if den == 0
    e = norm(a-b,inf);
else
    e = norm(a-b,inf)/den;
end
end

function plot_three_solver_metric(P_vec,delta_vec,left_data,dlp_data,std_data,y_label,fig_title,use_log)
figure('Color','w','Name',fig_title);
tiledlayout(1,numel(delta_vec),'TileSpacing','compact','Padding','compact');
for id = 1:numel(delta_vec)
    ax = nexttile;
    hold(ax,'on');
    plot(ax,P_vec,left_data(id,:),'o-','DisplayName','left weighted');
    plot(ax,P_vec,dlp_data(id,:),'s-','DisplayName','right DLP');
    plot(ax,P_vec,std_data(id,:),'^-','DisplayName','standard');
    grid(ax,'on');
    box(ax,'on');
    title(ax,sprintf('\\delta = %.3g',delta_vec(id)));
    xlabel(ax,'P');
    ylabel(ax,y_label);
    if use_log
        set(ax,'YScale','log');
    end
    if id == 1
        legend(ax,'Location','best');
    end
end
sgtitle(fig_title,'Interpreter','none');
end
