function results = jul16_fmm_breakeven_benchmark(config)
%JUL16_FMM_BREAKEVEN_BENCHMARK Find direct/FMM matvec crossovers.
%
%   RESULTS = JUL16_FMM_BREAKEVEN_BENCHMARK() compares direct evaluation
%   (vars.fmm = false) with FMM3D (vars.fmm = true) for
%
%       1. G:     Stokeslet flow from proxy nodes to surface nodes,
%       2. A:     A = T*W*G, evaluated by apply_Amat,
%       3. A_sym: (A + A')/2, evaluated by apply_sym_Amat.
%
%   The default discretization requests approx N = 700 proxy nodes per sphere and
%   uses a_glob = 1.2 for the outer surface grid. 
%
%   RESULTS = JUL16_FMM_BREAKEVEN_BENCHMARK(CONFIG) overrides any default
%   configuration field. For example, a quick smoke test is
%
%       cfg.P_vec = [1 2];
%       cfg.n_repeats = 1;
%       cfg.make_plots = false;
%       cfg.save_results = false;
%       results = jul16_fmm_breakeven_benchmark(cfg);
%
%   Timings target repeated matvec use: one untimed call per operator and
%   backend is made before the sweep. Direct/FMM timing order alternates
%   between repetitions to reduce ordering bias.
%
%   Important: apply_sym_Amat always evaluates its traction term directly.
%   Thus the A_sym comparison changes its G, T and G' evaluations to FMM but
%   retains one common O(P^2) direct traction evaluation.
%
%   Created July 16, 2026.

if nargin < 1
    config = struct();
end

script_dir = fileparts(mfilename('fullpath'));
repo_root = find_repo_root(script_dir);
run(fullfile(repo_root,'setpath.m'));

defaults = default_config(repo_root);
cfg = merge_config(defaults,config);
validate_config(cfg);
check_dependencies();

rng(cfg.seed);
Pmax = max(cfg.P_vec);
q = compact_lattice_centers(Pmax,cfg.center_spacing);

[rin_all,rout_all,opt] = init_spheres( ...
    q,cfg.proxy_radius,cfg.N_request,cfg.a_glob);

N = size(rin_all,1)/Pmax;
M = size(rout_all,1)/Pmax;
nout_all = rout_all - kron(q,ones(M,1));
nout_all = nout_all ./ vecnorm(nout_all,2,2);
wout_all = (4*pi/M) * ones(M*Pmax,1);
x_all = randn(3*N*Pmax,1);

operators = get_operators();
n_operator = numel(operators);
n_P = numel(cfg.P_vec);

raw_direct = nan(n_operator,n_P,cfg.n_repeats);
raw_fmm = nan(n_operator,n_P,cfg.n_repeats);
relative_error = nan(n_operator,n_P);

fprintf('Jul 16 direct/FMM break-even benchmark\n');
fprintf('  requested N = %d, actual N = %d\n',cfg.N_request,N);
fprintf('  a_glob = %.3g, actual M = %d, actual M/N = %.6f\n', ...
    cfg.a_glob,M,M/N);
fprintf('  proxy radius = %.3g, center spacing = %.3g\n', ...
    cfg.proxy_radius,cfg.center_spacing);
fprintf('  FMM tolerance = %.1e, repetitions = %d\n', ...
    cfg.fmm_tol,cfg.n_repeats);
fprintf('  P = [%s]\n\n',num2str(cfg.P_vec));

if cfg.warm_up
    warm_up_operators(operators,rin_all(1:N,:),rout_all(1:M,:), ...
        nout_all(1:M,:),wout_all(1:M),x_all(1:3*N),opt,cfg.fmm_tol);
end

fprintf('%6s %-8s | %11s %11s %9s %12s\n', ...
    'P','operator','direct (s)','FMM (s)','speedup','relative err');

for ip = 1:n_P
    P = cfg.P_vec(ip);
    rin = rin_all(1:N*P,:);
    rout = rout_all(1:M*P,:);
    nout = nout_all(1:M*P,:);
    wout = wout_all(1:M*P);
    x = x_all(1:3*N*P);

    vars_direct = opt;
    vars_direct.fmm = false;
    vars_direct.fmm_tol = cfg.fmm_tol;

    vars_fmm = vars_direct;
    vars_fmm.fmm = true;

    for io = 1:n_operator
        y_direct = [];
        y_fmm = [];

        for ir = 1:cfg.n_repeats
            if mod(ir,2) == 1
                [raw_direct(io,ip,ir),y_direct] = timed_apply( ...
                    operators(io).apply,x,rin,rout,nout,wout,vars_direct);
                [raw_fmm(io,ip,ir),y_fmm] = timed_apply( ...
                    operators(io).apply,x,rin,rout,nout,wout,vars_fmm);
            else
                [raw_fmm(io,ip,ir),y_fmm] = timed_apply( ...
                    operators(io).apply,x,rin,rout,nout,wout,vars_fmm);
                [raw_direct(io,ip,ir),y_direct] = timed_apply( ...
                    operators(io).apply,x,rin,rout,nout,wout,vars_direct);
            end
        end

        relative_error(io,ip) = relerr(y_fmm,y_direct);
        direct_time = median(raw_direct(io,ip,:),3);
        fmm_time = median(raw_fmm(io,ip,:),3);

        fprintf('%6d %-8s | %11.4g %11.4g %9.3f %12.3e\n', ...
            P,operators(io).short_name,direct_time,fmm_time, ...
            direct_time/fmm_time,relative_error(io,ip));
    end
end

direct_time = median(raw_direct,3);
fmm_time = median(raw_fmm,3);
speedup = direct_time ./ fmm_time;

crossover = repmat(struct( ...
    'operator','', ...
    'bracket',[nan nan], ...
    'estimate',nan, ...
    'status',''),n_operator,1);

fprintf('\nEstimated FMM break-even points\n');
for io = 1:n_operator
    crossover(io) = estimate_crossover( ...
        operators(io).name,cfg.P_vec,direct_time(io,:),fmm_time(io,:));
    print_crossover(crossover(io));
end

results = struct();
results.created = '2026-07-16';
results.description = [ ...
    'Direct versus FMM3D break-even benchmark for G, A = T*W*G, ', ...
    'and A_sym on sphere proxy/surface grids.'];
results.config = cfg;
results.actual.N = N;
results.actual.M = M;
results.actual.M_over_N = M/N;
results.centers = q;
results.operator_names = {operators.name};
results.operator_short_names = {operators.short_name};
results.time.raw_direct = raw_direct;
results.time.raw_fmm = raw_fmm;
results.time.direct = direct_time;
results.time.fmm = fmm_time;
results.speedup = speedup;
results.relative_error = relative_error;
results.crossover = crossover;
results.notes = { ...
    'Timings exclude one warm-up call per operator and backend when warm_up is true.', ...
    'A_sym retains a direct traction evaluation for both backend settings.', ...
    'The compact cubic lattice is deterministic; geometry can shift FMM crossover timing.'};

if cfg.save_results
    save(cfg.results_file,'results');
    fprintf('\nSaved results to %s\n',cfg.results_file);
end

if cfg.make_plots
    plot_results(results);
end

end


function cfg = default_config(repo_root)
cfg = struct();
cfg.N_request = 700;
cfg.a_glob = 1.2;
cfg.proxy_radius = 0.68;
cfg.P_vec = [40 48 56 64 72 80 88 96];
cfg.center_spacing = 2.5;
cfg.fmm_tol = 1e-10;
cfg.n_repeats = 3;
cfg.seed = 716;
cfg.warm_up = true;
cfg.make_plots = true;
cfg.save_results = true;
cfg.results_file = fullfile( ...
    repo_root,'data','jul16_fmm_breakeven_benchmark_results.mat');
end


function cfg = merge_config(defaults,config)
cfg = defaults;
names = fieldnames(config);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(defaults,name)
        error('jul16_fmm_breakeven_benchmark:unknownConfigField', ...
            'Unknown configuration field "%s".',name);
    end
    cfg.(name) = config.(name);
end
end


function validate_config(cfg)
validateattributes(cfg.N_request,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.a_glob,{'numeric'},{'scalar','>',1});
validateattributes(cfg.proxy_radius,{'numeric'},{'scalar','positive','<',1});
validateattributes(cfg.P_vec,{'numeric'},{'vector','integer','positive'});
validateattributes(cfg.center_spacing,{'numeric'},{'scalar','>',2});
validateattributes(cfg.fmm_tol,{'numeric'},{'scalar','positive'});
validateattributes(cfg.n_repeats,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.seed,{'numeric'},{'scalar','integer','nonnegative'});

if any(diff(cfg.P_vec) <= 0)
    error('jul16_fmm_breakeven_benchmark:invalidPVec', ...
        'config.P_vec must be strictly increasing.');
end
end


function check_dependencies()
dependencies = { ...
    'stfmm3d', ...
    'SE0P_Stokeslet_direct_full_ext_mex', ...
    'SE0P_Stresslet_direct_full_ext_mex', ...
    'SE0P_Stokestraction_direct_full_ext_mex'};

for k = 1:numel(dependencies)
    if exist(dependencies{k},'file') == 0
        error('jul16_fmm_breakeven_benchmark:missingDependency', ...
            'Required dependency "%s" is not on the MATLAB path.', ...
            dependencies{k});
    end
end
end


function operators = get_operators()
operators = struct([]);

operators(1).name = 'G: proxy to surface Stokeslet';
operators(1).short_name = 'G';
operators(1).apply = @apply_G;

operators(2).name = 'A = T*W*G';
operators(2).short_name = 'A';
operators(2).apply = @apply_A;

operators(3).name = 'A_sym = (A + A'')/2';
operators(3).short_name = 'A_sym';
operators(3).apply = @apply_A_sym;
end


function y = apply_G(x,rin,rout,~,~,vars)
y = getStokesletFlow(x,rin,rout,vars);
end


function y = apply_A(x,rin,rout,nout,wout,vars)
y = apply_Amat(x,rin,rout,nout,wout,vars);
end


function y = apply_A_sym(x,rin,rout,nout,wout,vars)
y = apply_sym_Amat(x,rin,rout,nout,wout,vars);
end


function warm_up_operators(operators,rin,rout,nout,wout,x,opt,fmm_tol)
for use_fmm = [false true]
    vars = opt;
    vars.fmm = use_fmm;
    vars.fmm_tol = fmm_tol;
    for io = 1:numel(operators)
        operators(io).apply(x,rin,rout,nout,wout,vars);
    end
end
end


function [elapsed,y] = timed_apply(apply,x,rin,rout,nout,wout,vars)
t_start = tic;
y = apply(x,rin,rout,nout,wout,vars);
elapsed = toc(t_start);
end


function e = relerr(a,b)
denominator = norm(b);
if denominator == 0
    e = norm(a-b);
else
    e = norm(a-b)/denominator;
end
end


function crossover = estimate_crossover(operator,P_vec,direct_time,fmm_time)
difference = direct_time - fmm_time;
first_fmm_win = find(difference > 0,1,'first');

crossover = struct();
crossover.operator = operator;
crossover.bracket = [nan nan];
crossover.estimate = nan;

if isempty(first_fmm_win)
    crossover.status = 'FMM was not faster in the tested range.';
    return;
end

if first_fmm_win == 1
    crossover.bracket = [nan P_vec(1)];
    crossover.estimate = P_vec(1);
    crossover.status = 'FMM was already faster at the first tested P.';
    return;
end

k0 = first_fmm_win - 1;
k1 = first_fmm_win;
crossover.bracket = P_vec([k0 k1]);

d0 = difference(k0);
d1 = difference(k1);
if d1 == d0
    crossover.estimate = mean(crossover.bracket);
else
    crossover.estimate = P_vec(k0) - ...
        d0 * (P_vec(k1)-P_vec(k0))/(d1-d0);
end
crossover.status = 'Crossover linearly interpolated within the bracket.';
end


function print_crossover(crossover)
if all(isfinite(crossover.bracket))
    fprintf('  %-30s: P in [%g, %g], estimate P ~= %.1f\n', ...
        crossover.operator,crossover.bracket(1),crossover.bracket(2), ...
        crossover.estimate);
else
    fprintf('  %-30s: %s\n',crossover.operator,crossover.status);
end
end


function q = compact_lattice_centers(P,spacing)
nside = ceil(P^(1/3));
[ix,iy,iz] = ndgrid(0:nside-1,0:nside-1,0:nside-1);
lattice = [ix(:),iy(:),iz(:)];
lattice = lattice - (nside-1)/2;

[~,order] = sort(sum(lattice.^2,2),'ascend');
q = spacing * lattice(order(1:P),:);
q = q - mean(q,1);
end


function plot_results(results)
P_vec = results.config.P_vec;
colors = lines(numel(results.operator_names));

figure('Color','w','Name','Jul 16 direct/FMM break-even benchmark');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile;
hold(ax,'on');
for io = 1:numel(results.operator_names)
    semilogy(ax,P_vec,results.time.direct(io,:),'o-', ...
        'Color',colors(io,:), ...
        'DisplayName',sprintf('%s direct',results.operator_short_names{io}));
    semilogy(ax,P_vec,results.time.fmm(io,:),'s--', ...
        'Color',colors(io,:), ...
        'DisplayName',sprintf('%s FMM',results.operator_short_names{io}));
end
grid(ax,'on');
box(ax,'on');
xlabel(ax,'number of particles, P');
ylabel(ax,'median matvec time (s)');
title(ax,'Direct and FMM timings');
legend(ax,'Location','northwest');

ax = nexttile;
hold(ax,'on');
for io = 1:numel(results.operator_names)
    plot(ax,P_vec,results.speedup(io,:),'o-', ...
        'Color',colors(io,:), ...
        'DisplayName',results.operator_short_names{io});
end
yline(ax,1,'k--','FMM break-even','HandleVisibility','off');
grid(ax,'on');
box(ax,'on');
xlabel(ax,'number of particles, P');
ylabel(ax,'direct time / FMM time');
title(ax,'FMM speedup');
legend(ax,'Location','northwest');
end


function repo_root = find_repo_root(start_dir)
repo_root = start_dir;
while true
    if isfile(fullfile(repo_root,'setpath.m'))
        return;
    end

    parent_dir = fileparts(repo_root);
    if isempty(parent_dir) || strcmp(parent_dir,repo_root)
        error('jul16_fmm_breakeven_benchmark:repoRootNotFound', ...
            'Could not find setpath.m starting from %s.',start_dir);
    end
    repo_root = parent_dir;
end
end
