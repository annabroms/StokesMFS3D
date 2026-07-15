%JUL13_EXAMPLE2_MOBILITY_VS_P 
%
% This script follows the sphere-cluster setup from Example 2 / Fig. 7 of
% Broms, Barnett, and Tornberg (Adv. Comput. Math., 2025), but compares only
% the two mobility solvers in this repository:
%   1. solve_mobility
%   2. solve_mobility_with_DLP
%
% A resistance solve is still used to generate F from a known Uref, since
% the paper's 2-way error is defined by resistance followed by mobility.
% The plotted quantities are only mobility quantities.
%
% Created July 13, 2026. 

clear;
close all;
clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
start_dir = pwd;
cleanup_obj = onCleanup(@() cd(start_dir));

cd(repo_root);
run(fullfile(repo_root,'setpath.m'));

% get_sphdesign.m reads spherical_designs/filelist.txt relative to pwd.
cd(fullfile(repo_root,'geometry','spheres'));

%% Example 2 configuration
P_vec = 10:10:200;
delta_vec = [0.2 1.0];
base_seed = 713;

gmres_tol = 1e-7;
gmres_tol = 1e-10;
maxit = 500;
use_fmm = true;
fmm_tol = 1e-8;

dlp_inner_only = false;
dlp_symmetrize_weighted = false;
dlp_add_rank1 = false;
dlp_outer_force = true;

% In the paper, 2-way error uses distinct resistance and mobility proxy
% offsets with sep_mob = 1.05*sep_res. Fig. 7 reports Rp_res = 0.7.
mobility_offset_factor = 1.05;

resolutions = struct([]);
resolutions(1).name = 'low';
resolutions(1).N_request = 100;
resolutions(1).a = 1.2;
resolutions(1).Rp_mob = 0.30;
resolutions(1).Rp_res = 1 - (1 - resolutions(1).Rp_mob)/mobility_offset_factor;

resolutions(2).name = 'paper high';
resolutions(2).N_request = 700;
resolutions(2).a = 1.2;
resolutions(2).Rp_res = 0.63;
resolutions(2).Rp_mob = 1 - mobility_offset_factor*(1 - resolutions(2).Rp_res);


method_names = {'solve_mobility','solve_mobility_with_DLP'};

%% Storage
n_method = numel(method_names);
n_res = numel(resolutions);
n_delta = numel(delta_vec);
n_P = numel(P_vec);

try
    results_file = fullfile(repo_root,'experiments','jul13_example2_mobility_dlp_resolution_compare_results.mat');
    load(results_file,'results');
catch
    results = struct();
    results.created = '2026-07-13';
    results.description = 'Example 2 sphere-cluster mobility comparison: solve_mobility vs solve_mobility_with_DLP.';
    results.P_vec = P_vec;
    results.delta_vec = delta_vec;
    results.resolutions = resolutions;
    results.method_names = method_names;
    results.gmres_tol = gmres_tol;
    results.maxit = maxit;
    results.use_fmm = use_fmm;
    results.fmm_threshold = fmm_threshold;
    results.fmm_tol = fmm_tol;
    results.dlp_inner_only = dlp_inner_only;
    results.dlp_symmetrize_weighted = dlp_symmetrize_weighted;
    results.dlp_add_rank1 = dlp_add_rank1;
    results.dlp_outer_force = dlp_outer_force;
    
    results.actual.N_res = nan(n_res,n_delta,n_P);
    results.actual.M_res = nan(n_res,n_delta,n_P);
    results.actual.N_mob = nan(n_res,n_delta,n_P);
    results.actual.M_mob = nan(n_res,n_delta,n_P);
    
    results.resistance.iters = nan(n_res,n_delta,n_P);
    results.resistance.uerr = nan(n_res,n_delta,n_P);
    results.resistance.lambda_norm = nan(n_res,n_delta,n_P);
    results.resistance.solve_time = nan(n_res,n_delta,n_P);
    
    results.mobility.iters = nan(n_method,n_res,n_delta,n_P);
    results.mobility.uerr = nan(n_method,n_res,n_delta,n_P);
    results.mobility.two_way = nan(n_method,n_res,n_delta,n_P);
    results.mobility.lambda_norm = nan(n_method,n_res,n_delta,n_P);
    results.mobility.solve_time = nan(n_method,n_res,n_delta,n_P);
    results.mobility.coeff_count = nan(n_method,n_res,n_delta,n_P);
    results.mobility.gmres_unknown_count = nan(n_method,n_res,n_delta,n_P);
end

fprintf('Jul 13 Example 2 mobility comparison\n');
fprintf('  P     = [%s]\n',num2str(P_vec));
fprintf('  delta = [%s]\n',num2str(delta_vec));
fprintf('  gmres_tol = %.1e, maxit = %d\n',gmres_tol,maxit);

fprintf('%12s %8s %5s %6s %6s | %8s %10s %8s | %8s %10s %10s %10s | %8s %10s %10s %10s\n', ...
    'resolution','delta','P','Nmob','Mmob', ...
    'it res','uerr res','t res', ...
    'it std','uerr std','2-way std','|lam| std', ...
    'it DLP','uerr DLP','2-way DLP','|lam| DLP');

%% Main sweep
for id = 1:n_delta
    delta = delta_vec(id);

    for ip = 1:n_P
        P = P_vec(ip);
        rng(base_seed + 10000*id + P);

        [q,~] = grow_cluster(P,delta);
        Uref = rand(6*P,1);

        for ir = 1:n_res
            if ir == 2
                fmm_threshold = 70;
            else
                fmm_threshold = 700;
            end
            if ~isnan(results.actual.N_res(ir,id,ip))
                continue;
            end
            cfg = resolutions(ir);

            [rres_in,rres_out,opt_res] = init_spheres(q,cfg.Rp_res,cfg.N_request,cfg.a);
            opt_res = configure_options(opt_res,P,gmres_tol,maxit,use_fmm,fmm_threshold,fmm_tol);
            opt_res.lr = 0;

            N_res = size(rres_in,1)/P;
            M_res = size(rres_out,1)/P;
            results.actual.N_res(ir,id,ip) = N_res;
            results.actual.M_res(ir,id,ip) = M_res;

            t_start = tic;
            [Fvec,it_res,ln_res,uerr_res] = solve_resistance(q,rres_in,rres_out,Uref,opt_res);
            t_res = toc(t_start);

            results.resistance.iters(ir,id,ip) = it_res;
            results.resistance.uerr(ir,id,ip) = uerr_res;
            results.resistance.lambda_norm(ir,id,ip) = ln_res;
            results.resistance.solve_time(ir,id,ip) = t_res;

            [rmob_in,rmob_out,opt_mob] = init_spheres(q,cfg.Rp_mob,cfg.N_request,cfg.a);
            opt_mob = configure_options(opt_mob,P,gmres_tol,maxit,use_fmm,fmm_threshold,fmm_tol);
            opt_mob.lr = 0;
            opt_mob.compute_residual = true;

            N_mob = size(rmob_in,1)/P;
            M_mob = size(rmob_out,1)/P;
            [nout,wout] = sphere_normals_and_weights(q,rmob_out);

            results.actual.N_mob(ir,id,ip) = N_mob;
            results.actual.M_mob(ir,id,ip) = M_mob;

            t_start = tic;
            [U_std,it_std,ln_std,uerr_std] = solve_mobility(q,rmob_in,rmob_out,Fvec,opt_mob);
            t_std = toc(t_start);

            results.mobility.iters(1,ir,id,ip) = it_std;
            results.mobility.uerr(1,ir,id,ip) = uerr_std;
            results.mobility.two_way(1,ir,id,ip) = relerr_inf(U_std,Uref);
            results.mobility.lambda_norm(1,ir,id,ip) = ln_std;
            results.mobility.solve_time(1,ir,id,ip) = t_std;
            results.mobility.coeff_count(1,ir,id,ip) = 3*N_mob*P;
            results.mobility.gmres_unknown_count(1,ir,id,ip) = 3*M_mob*P;

            opt_dlp = opt_mob;
            opt_dlp.inner_only = dlp_inner_only;
            opt_dlp.symmetrize_weighted = dlp_symmetrize_weighted;
            opt_dlp.add_rank1 = dlp_add_rank1;
            opt_dlp.outer_force = dlp_outer_force;

            t_start = tic;
            [U_dlp,it_dlp,ln_dlp,uerr_dlp] = solve_mobility_with_DLP( ...
                q,rmob_in,rmob_out,nout,wout,Fvec,opt_dlp);
            t_dlp = toc(t_start);

            results.mobility.iters(2,ir,id,ip) = it_dlp;
            results.mobility.uerr(2,ir,id,ip) = uerr_dlp;
            results.mobility.two_way(2,ir,id,ip) = relerr_inf(U_dlp,Uref);
            results.mobility.lambda_norm(2,ir,id,ip) = ln_dlp;
            results.mobility.solve_time(2,ir,id,ip) = t_dlp;
            results.mobility.coeff_count(2,ir,id,ip) = 3*N_mob*P;
            results.mobility.gmres_unknown_count(2,ir,id,ip) = 3*N_mob*P;

            fprintf('%12s %8.3g %5d %6d %6d | %8d %10.3e %8.3g | %8d %10.3e %10.3e %10.3e | %8d %10.3e %10.3e %10.3e\n', ...
                cfg.name,delta,P,N_mob,M_mob, ...
                it_res,uerr_res,t_res, ...
                it_std,uerr_std,results.mobility.two_way(1,ir,id,ip),ln_std, ...
                it_dlp,uerr_dlp,results.mobility.two_way(2,ir,id,ip),ln_dlp);

            results_file = fullfile(repo_root,'experiments','jul13_example2_mobility_dlp_resolution_compare_results.mat');
            save(results_file,'results');
        end
    end
end


%%
results_file = fullfile(repo_root,'data','jul13_example2_mobility_dlp_resolution_compare_results.mat');
%results_file = fullfile(repo_root,'experiments','jul13_example2_mobility_dlp_resolution_compare_results.mat');
load(results_file,'results');

P_vec = results.P_vec;
delta_vec = results.delta_vec;
resolutions = results.resolutions;
method_names = results.method_names;


figure();
plot_metric_panel(P_vec,delta_vec,resolutions,method_names, ...
    results.mobility.uerr,'relative boundary residual',true,true);
figure()
plot_metric_panel(P_vec,delta_vec,resolutions,method_names, ...
    results.mobility.two_way,'2-way error',true,true);
figure()
plot_metric_panel(P_vec,delta_vec,resolutions,method_names, ...
    results.mobility.lambda_norm,'coefficient vector size, ||\lambda||_\infty',true,true);
figure()
plot_metric_panel(P_vec,delta_vec,resolutions,method_names, ...
    results.mobility.iters,'GMRES iterations',false,true);

figure()
plot_metric_panel(P_vec,delta_vec,resolutions,method_names, ...
    results.mobility.solve_time,'Mobility solve time [s]',false,true);


sgtitle('Example 2 mobility: standard MFS vs DLP MFS','Interpreter','none');



function opt = configure_options(opt,P,gmres_tol,maxit,use_fmm,fmm_threshold,fmm_tol)
opt.fmm = use_fmm && P > fmm_threshold;
opt.eps = fmm_tol;
opt.lr = 0;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = 0;
opt.maxit = maxit;
opt.plot = 0;
opt.debug = 0;
opt.profile = 0;
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


function plot_metric_panel(P_vec,delta_vec,resolutions,method_names,data,y_label,use_log,show_legend)
colors = lines(numel(resolutions));
markers = {'o','s'};
line_styles = {'-','--',':'};

for ir = 1:numel(resolutions)
    for im = 1:numel(method_names)
        for id = 1:numel(delta_vec)
            marker = markers{min(im,numel(markers))};
            line_style = line_styles{min(id,numel(line_styles))};
            display_name = sprintf('%s, %s, delta=%.3g', ...
                short_method_name(method_names{im}),resolutions(ir).name,delta_vec(id));

            y = squeeze(data(im,ir,id,:));
            plot(P_vec,y, ...
                'Color',colors(ir,:), ...
                'LineStyle',line_style, ...
                'Marker',marker, ...
                'DisplayName',display_name);
            hold on
        end
    end
end

grid('on');
box('on');
xlabel('P');
ylabel(y_label);
ax = gca; 
if use_log
    set(ax,'YScale','log');
end
if show_legend
    legend(ax,'Location','bestoutside');
end
end

function name = short_method_name(method_name)
switch method_name
    case 'solve_mobility'
        name = 'standard';
    case 'solve_mobility_with_DLP'
        name = 'DLP';
    otherwise
        name = method_name;
end
end
