%AUG05_DEFLATION_PCG_GAP_SWEEP Compare standard DLP GMRES with deflated PCG.
%
% Two unit spheres receive equal and opposite squeezing forces.  At fixed
% spherical-design resolution, the script sweeps the surface gap and plots
% the maximum off-grid surface residual, Krylov iterations, and solve time.
% Both residuals use the same independently rotated spherical design.

clear;
close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
repo_root = find_repo_root(script_dir);
run(fullfile(repo_root,'setpath.m'));

rng(1);
gap_vec = [0.1 0.2 0.3 0.5 1 2 5];
N_request = 200;
Rp = 0.30;
a_glob = 1.2;
gmres_tol = 1.5e-10;
maxit = 600;
method_names = {'standard_dlp_gmres','deflated_split_pcg'};
method_labels = {'standard DLP GMRES','deflated split PCG'};

P = 2;
Fvec = zeros(6*P,1);
Fvec(3) = 1;
Fvec(9) = -1;

results = struct();
results.created = '2026-08-05';
results.gap = gap_vec;
results.N_request = N_request;
results.Rp = Rp;
results.a_glob = a_glob;
results.gmres_tol = gmres_tol;
results.maxit = maxit;
results.method_names = method_names;
results.N_actual = nan(size(gap_vec));
results.M_actual = nan(size(gap_vec));
results.iters = nan(2,numel(gap_vec));
results.max_relative_residual = nan(2,numel(gap_vec));
results.solve_time = nan(2,numel(gap_vec));
results.pcg_call_time = nan(1,numel(gap_vec));
results.pcg_converged = false(1,numel(gap_vec));
results.pcg_failure_reason = repmat({''},1,numel(gap_vec));
results.projected_residual = nan(1,numel(gap_vec));
results.smallest_eigenvalue = nan(1,numel(gap_vec));

fprintf('Aug 05 deflated PCG gap sweep\n');
fprintf('  N request = %d, Rp = %.2f, M/N target = %.2f\n', ...
    N_request,Rp,a_glob);
fprintf('  solver tolerance = %.1e, maxit = %d\n\n',gmres_tol,maxit);
fprintf(['%8s %5s %5s | %8s %12s %10s | ', ...
         '%8s %12s %10s %12s\n'], ...
    'gap','N','M','GMRES it','GMRES res','GMRES sec', ...
    'PCG it','PCG res','PCG sec','projected');

for igap = 1:numel(gap_vec)
    gap = gap_vec(igap);
    centre = 1 + gap/2;
    q = [0 0 -centre;0 0 centre];
    [rin,rout,opt] = init_spheres(q,Rp,N_request,a_glob);
    M = size(rout,1)/P;
    nout = rout - kron(q,ones(M,1));
    nout = nout ./ vecnorm(nout,2,2);
    wout = (4*pi/M) * ones(P*M,1);

    opt.fmm = false;
    opt.fmm_tol = 1e-12;
    opt.gmres_tol = gmres_tol;
    opt.solver_tol = gmres_tol;
    opt.gmres_verbose = 0;
    opt.maxit = maxit;
    opt.compute_residual = false;
    opt.symmetrize_weighted = true;
    opt.inner_only = false;
    opt.add_rank1 = false;
    opt.outer_force = true;

    standard_start = tic;
    [U_gmres,it_gmres,~,~,lambda_gmres] = solve_mobility_with_DLP( ...
        q,rin,rout,nout,wout,Fvec,[],opt);
    standard_time = toc(standard_start);
    residual_gmres = sphere_residual(q,rin,lambda_gmres,U_gmres,opt);

    pcg_opt = opt;
    pcg_opt.deflation_solver = 'pcg_split';
    N_actual = size(rin,1)/P;
    pcg_opt.deflation_check_nmax = max(150,ceil(1.6*N_actual));
    pcg_opt.deflation_verify_nmax = max(180,ceil(2.2*N_actual));
    pcg_start = tic;
    [U_pcg,it_pcg,~,~,lambda_pcg,pcg_out] = solve_with_deflation( ...
        q,rin,rout,nout,wout,Fvec,pcg_opt);
    pcg_time = toc(pcg_start);
    pcg_converged = ~pcg_out.used_fallback && ...
        strcmp(pcg_out.solver,'pcg_split');
    if pcg_converged
        residual_pcg = sphere_residual(q,rin,lambda_pcg,U_pcg,pcg_opt);
        pcg_iterations = it_pcg;
    else
        residual_pcg = NaN;
        pcg_iterations = pcg_out.krylov_iterations;
    end

    results.N_actual(igap) = N_actual;
    results.M_actual(igap) = M;
    results.iters(:,igap) = [it_gmres;pcg_iterations];
    results.max_relative_residual(:,igap) = ...
        [residual_gmres;residual_pcg];
    results.solve_time(1,igap) = standard_time;
    if pcg_converged
        results.solve_time(2,igap) = pcg_time;
    end
    results.pcg_call_time(igap) = pcg_time;
    results.pcg_converged(igap) = pcg_converged;
    results.pcg_failure_reason{igap} = pcg_out.fallback_reason;
    results.projected_residual(igap) = pcg_out.projected_residual;
    results.smallest_eigenvalue(igap) = pcg_out.smallest_eigenvalue;

    fprintf('%8.3g %5d %5d | %8d %12.3e %10.3f | ', ...
        gap,results.N_actual(igap),M,it_gmres,residual_gmres,standard_time);
    if pcg_converged
        fprintf('%8d %12.3e %10.3f %12.3e\n', ...
            pcg_iterations,residual_pcg,pcg_time,pcg_out.projected_residual);
        assert(all(isfinite([residual_pcg,pcg_time, ...
            pcg_out.projected_residual])), ...
            'A nonfinite PCG diagnostic was returned at gap %.3g.',gap);
    else
        fprintf('%8d %12s %10.3f %12s  %s ', ...
            pcg_iterations,'failed',pcg_time,'n/a',pcg_out.fallback_reason);
        fprintf('(lambda_min %.3e, tol %.3e)\n', ...
            pcg_out.smallest_eigenvalue,pcg_out.curvature_tol);
    end
    assert(all(isfinite([residual_gmres,standard_time])), ...
        'A nonfinite GMRES diagnostic was returned at gap %.3g.',gap);
end

if any(~results.pcg_converged)
    fprintf(['\n  PCG failed at %d of %d gaps. Those entries are NaN in ', ...
             'the PCG residual/time curves; fallback GMRES values are not ', ...
             'plotted as PCG.\n'], ...
        nnz(~results.pcg_converged),numel(gap_vec));
end

figure('Color','w','Name','Aug 05 DLP GMRES versus deflated PCG');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

ax = nexttile;
loglog(ax,gap_vec,results.max_relative_residual(1,:),'o-', ...
    'LineWidth',1.2,'DisplayName',method_labels{1});
hold(ax,'on');
loglog(ax,gap_vec,results.max_relative_residual(2,:),'s-', ...
    'LineWidth',1.2,'DisplayName',method_labels{2});
grid(ax,'on');
box(ax,'on');
ylabel(ax,'max relative residual');
legend(ax,'Location','best');

ax = nexttile;
semilogx(ax,gap_vec,results.iters(1,:),'o-', ...
    'LineWidth',1.2,'DisplayName',method_labels{1});
hold(ax,'on');
semilogx(ax,gap_vec,results.iters(2,:),'s-', ...
    'LineWidth',1.2,'DisplayName',method_labels{2});
grid(ax,'on');
box(ax,'on');
ylabel(ax,'Krylov iterations');

ax = nexttile;
loglog(ax,gap_vec,results.solve_time(1,:),'o-', ...
    'LineWidth',1.2,'DisplayName',method_labels{1});
hold(ax,'on');
loglog(ax,gap_vec,results.solve_time(2,:),'s-', ...
    'LineWidth',1.2,'DisplayName',method_labels{2});
grid(ax,'on');
box(ax,'on');
xlabel(ax,'surface gap');
ylabel(ax,'solve time (s)');

sgtitle('Two-sphere mobility: standard DLP GMRES versus deflated split PCG');

function repo_root = find_repo_root(start_dir)
% Find setpath.m without assuming the script was launched from the repo root.
repo_root = start_dir;
while ~isfile(fullfile(repo_root,'setpath.m'))
    parent = fileparts(repo_root);
    if isempty(parent) || strcmp(parent,repo_root)
        error('aug05_deflation_pcg_gap_sweep:repoRootNotFound', ...
            'Could not find setpath.m above %s.',start_dir);
    end
    repo_root = parent;
end
end

function uerr = sphere_residual(q,rin,lambda,U,opt)
% Evaluate a common off-grid residual on a rotated spherical design.
[unit_points,~] = get_sphdesign(180);
alpha = 0.31;
beta = -0.43;
Rz = [cos(alpha) -sin(alpha) 0; ...
      sin(alpha)  cos(alpha) 0; ...
      0           0          1];
Ry = [cos(beta) 0 sin(beta); ...
      0         1 0; ...
     -sin(beta) 0 cos(beta)];
unit_points = ((Rz*Ry) * unit_points')';

P = size(q,1);
Mcheck = size(unit_points,1);
rcheck = zeros(P*Mcheck,3);
ucheck = zeros(3*P*Mcheck,1);
for k = 1:P
    point_rows = (k-1)*Mcheck + (1:Mcheck);
    vector_rows = (k-1)*3*Mcheck + (1:3*Mcheck);
    rcheck(point_rows,:) = q(k,:) + unit_points;
    Kcheck = getKmat(rcheck(point_rows,:),q(k,:));
    ucheck(vector_rows) = Kcheck * U((k-1)*6+(1:6));
end

uboundary = getStokesletFlow(lambda,rin,rcheck,opt);
node_error = vecnorm(reshape(ucheck-uboundary,3,[]),2,1);
velocity_scale = max(vecnorm(reshape(ucheck,3,[]),2,1));
uerr = max(node_error) / max(velocity_scale,eps);
end
