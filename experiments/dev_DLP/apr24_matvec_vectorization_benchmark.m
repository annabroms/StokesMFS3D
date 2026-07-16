%Benchmark looped and vectorized matvec kernels.
%
%   APR24_MATVEC_VECTORIZATION_BENCHMARK() benchmarks the
%   per-particle preconditioner mapping and self-interaction correction
%   kernels used in matvecStokesMFS and matvecStokesMFS_DLP. The script
%   compares the current looped implementations with candidate vectorized
%   versions for a range of P and reports timing, relative error, and a
%   switch/keep recommendation.

rng(24);

P_list = [8 32 64 100];
speedup_threshold = 1.10;
equiv_tol = 1e-10;

fprintf('Apr 24 matvec vectorization benchmark\n');
fprintf('  P list            : [%s]\n', num2str(P_list));
fprintf('  speedup threshold : %.2fx\n', speedup_threshold);
fprintf('  equivalence tol   : %.1e\n\n', equiv_tol);

sphere = build_sphere_case(max(P_list));
ellipsoid = build_ellipsoid_case(max(P_list));

sep = repmat('-', 1, 70);

results = struct();
fprintf('%s\n', sep);
results.standard_precond = benchmark_precond_family( ...
    'standard_precond', sphere, ellipsoid, P_list, speedup_threshold, equiv_tol);
fprintf('%s\n', sep);
results.dlp_precond = benchmark_precond_family( ...
    'dlp_precond', sphere, ellipsoid, P_list, speedup_threshold, equiv_tol);
fprintf('%s\n', sep);
results.standard_self = benchmark_self_family( ...
    'standard_self', sphere, ellipsoid, P_list, speedup_threshold, equiv_tol);
fprintf('%s\n', sep);
results.dlp_self = benchmark_self_family( ...
    'dlp_self', sphere, ellipsoid, P_list, speedup_threshold, equiv_tol);

fprintf('\n%s\nSummary by branch\n%s\n', sep, sep);
print_branch_summary(results.standard_precond, 'matvecStokesMFS preconditioner');
print_branch_summary(results.dlp_precond, 'matvecStokesMFS_DLP preconditioner');
print_branch_summary(results.standard_self, 'matvecStokesMFS self correction');
print_branch_summary(results.dlp_self, 'matvecStokesMFS_DLP self correction');



function family_results = benchmark_precond_family(name, sphere, ellipsoid, P_list, speedup_threshold, equiv_tol)
family_results = struct();
family_results.name = name;
family_results.sphere.resistance = benchmark_precond_case( ...
    name, 'sphere', 'resistance', sphere, P_list, speedup_threshold, equiv_tol);
family_results.sphere.mobility = benchmark_precond_case( ...
    name, 'sphere', 'mobility', sphere, P_list, speedup_threshold, equiv_tol);
family_results.ellipsoid.resistance = benchmark_precond_case( ...
    name, 'ellipsoid', 'resistance', ellipsoid, P_list, speedup_threshold, equiv_tol);
family_results.ellipsoid.mobility = benchmark_precond_case( ...
    name, 'ellipsoid', 'mobility', ellipsoid, P_list, speedup_threshold, equiv_tol);
end

function family_results = benchmark_self_family(name, sphere, ellipsoid, P_list, speedup_threshold, equiv_tol)
family_results = struct();
family_results.name = name;
family_results.sphere = benchmark_self_case( ...
    name, 'sphere', sphere, P_list, speedup_threshold, equiv_tol);
family_results.ellipsoid = benchmark_self_case( ...
    name, 'ellipsoid', ellipsoid, P_list, speedup_threshold, equiv_tol);
end

function rows = benchmark_precond_case(name, geometry_name, mode_name, data, P_list, speedup_threshold, equiv_tol)
fprintf('%s | %s | %s\n', upper(name), upper(geometry_name), upper(mode_name));
fprintf('  %-5s | %-14s | %-14s | %-10s | %-12s | %s\n', ...
    'P', 'loop (s)', 'vec (s)', 'speedup', 'rel err', 'recommendation');

rows = repmat(struct('P', [], 'loop_time', [], 'vec_time', [], ...
    'speedup', [], 'rel_err', [], 'recommendation', ''), numel(P_list), 1);

for k = 1:numel(P_list)
    P = P_list(k);
    case_data = slice_case_data(data, P);
    if strcmp(name, 'standard_precond')
        loop_fn = @() standard_precond_loop(case_data.mu_standard, ...
            case_data.std.(mode_name).U, case_data.std.(mode_name).Y, ...
            case_data.std.mobility.L, case_data.vars_standard, case_data.R, ...
            strcmp(mode_name, 'resistance'));
        vec_fn = @() standard_precond_vectorized(case_data.mu_standard, ...
            case_data.std.(mode_name).U, case_data.std.(mode_name).Y, ...
            case_data.std.mobility.L, case_data.vars_standard, case_data.R, ...
            strcmp(mode_name, 'resistance'));
    else
        loop_fn = @() dlp_precond_loop(case_data.mu_dlp, ...
            case_data.dlp.(mode_name).U, case_data.dlp.(mode_name).Y, ...
            case_data.dlp.mobility.L, case_data.vars_dlp, case_data.R, ...
            strcmp(mode_name, 'resistance'));
        vec_fn = @() dlp_precond_vectorized(case_data.mu_dlp, ...
            case_data.dlp.(mode_name).U, case_data.dlp.(mode_name).Y, ...
            case_data.dlp.mobility.L, case_data.vars_dlp, case_data.R, ...
            strcmp(mode_name, 'resistance'));
    end

    loop_out = loop_fn();
    vec_out = vec_fn();
    rel_err = relative_error(loop_out, vec_out);

    loop_fn();
    vec_fn();
    loop_time = timeit(loop_fn);
    vec_time = timeit(vec_fn);
    speedup = loop_time / vec_time;
    recommendation = recommend_switch(rel_err, speedup, speedup_threshold, equiv_tol);

    rows(k) = make_row(P, loop_time, vec_time, speedup, rel_err, recommendation);
    print_row(rows(k));
end

fprintf('\n');
end

function rows = benchmark_self_case(name, geometry_name, data, P_list, speedup_threshold, equiv_tol)
fprintf('%s | %s\n', upper(name), upper(geometry_name));
fprintf('  %-5s | %-14s | %-14s | %-10s | %-12s | %s\n', ...
    'P', 'loop (s)', 'vec (s)', 'speedup', 'rel err', 'recommendation');

rows = repmat(struct('P', [], 'loop_time', [], 'vec_time', [], ...
    'speedup', [], 'rel_err', [], 'recommendation', ''), numel(P_list), 1);

for k = 1:numel(P_list)
    P = P_list(k);
    case_data = slice_case_data(data, P);
    if strcmp(name, 'standard_self')
        loop_fn = @() standard_self_loop(case_data.lambda_stokes, ...
            case_data.rin, case_data.rout, case_data.vars_standard, case_data.Sblock_loop);
        vec_fn = @() standard_self_vectorized(case_data.lambda_stokes, ...
            case_data.vars_standard, case_data.R, case_data.Sblock_vec);
    else
        loop_fn = @() dlp_self_loop(case_data.lambda_stokes, ...
            case_data.rin, case_data.rout, case_data.vars_dlp, case_data.Tblock_loop, ...
            case_data.R, case_data.TSblock_vec);
        vec_fn = @() dlp_self_vectorized(case_data.lambda_stokes, ...
            case_data.vars_dlp, case_data.R, case_data.TSblock_vec);
    end

    loop_out = loop_fn();
    vec_out = vec_fn();
    rel_err = relative_error(loop_out, vec_out);

    loop_fn();
    vec_fn();
    loop_time = timeit(loop_fn);
    vec_time = timeit(vec_fn);
    speedup = loop_time / vec_time;
    recommendation = recommend_switch(rel_err, speedup, speedup_threshold, equiv_tol);

    rows(k) = make_row(P, loop_time, vec_time, speedup, rel_err, recommendation);
    print_row(rows(k));
end

fprintf('\n');
end

function data = build_sphere_case(Pmax)
q = [(0:Pmax-1)' * 4, zeros(Pmax, 2)];
[rin, rout, opt] = init_spheres(q, 0.68, 50, 1.2);

N = size(rin, 1) / Pmax;
M = size(rout, 1) / Pmax;

rin1 = rin(1:N,:);
rout1 = rout(1:M,:);
nout1 = rout1;
nout1 = nout1 ./ vecnorm(nout1,2,2);
wout1 = (4*pi/M) * ones(M,1);

[Yres_std, Ures_std] = oneBodyPrecondRes(rin1, rout1);
[Ymob_std, Umob_std, Lstd] = oneBodyPrecondMob(rin1, rout1, [0 0 0]);
[Yres_dlp, Ures_dlp] = oneBodyPrecondResDLP(rin1, rout1, rout1);
[Ymob_dlp, Umob_dlp, Ldlp] = oneBodyPrecondMobDLP( ...
    rin1, rout1, struct(), [0 0 0], wout1, nout1);

Sblock = stokes_SLP_mat(rin1, rout1);
Tblock = stokes_DLP_mat(rout1, rin1, rout1);

data.Pmax = Pmax;
data.q = q;
data.R = repmat({eye(3)}, 1, Pmax);
data.rin = rin;
data.rout = rout;
data.N = N;
data.M = M;
data.std.resistance.U = Ures_std;
data.std.resistance.Y = Yres_std;
data.std.mobility.U = Umob_std;
data.std.mobility.Y = Ymob_std;
data.std.mobility.L = Lstd;
data.dlp.resistance.U = Ures_dlp;
data.dlp.resistance.Y = Yres_dlp;
data.dlp.mobility.U = Umob_dlp;
data.dlp.mobility.Y = Ymob_dlp;
data.dlp.mobility.L = Ldlp;
data.Sblock_vec = Sblock;
data.Sblock_loop = Sblock;
data.Tblock_loop = Tblock;
data.TSblock_vec = Tblock * Sblock;
data.mu_standard = randn(3*M*Pmax, 1);
data.mu_dlp = randn(3*N*Pmax, 1);
data.lambda_stokes = randn(3*N*Pmax, 1);
data.vars_standard = benchmark_vars(opt, M, N, false, Sblock, []);
data.vars_dlp = benchmark_vars(opt, M, N, false, Sblock, Tblock * Sblock);
end

function data = build_ellipsoid_case(Pmax)
E0 = [0.4 0.6 1.0];
Nv = 15;
N1 = 0.75 * Nv;
sep = 0.125;
q = [(0:Pmax-1)' * 4, zeros(Pmax, 2)];
R = random_rotations(Pmax);
[rin, rout, ~, ~, ~, wout, nout] = getEllipsoidGrids(E0, Pmax, 1, N1, Nv, sep, R, q);

N = size(rin, 1) / Pmax;
M = size(rout, 1) / Pmax;

rin_body = (R{1}' * (rin(1:N,:) - q(1,:))')';
rout_body = (R{1}' * (rout(1:M,:) - q(1,:))')';
nout_body = (R{1}' * nout(1:M,:)')';

[Yres_std, Ures_std] = oneBodyPrecondRes(rin_body, rout_body);
[Ymob_std, Umob_std, Lstd] = oneBodyPrecondMob(rin_body, rout_body, [0 0 0]);
[Ymob_dlp, Umob_dlp, Ldlp] = oneBodyPrecondMobDLP( ...
    rin_body, rout_body, struct(), [0 0 0], wout(1:M), nout_body);

Sblock = stokes_SLP_mat(rin_body, rout_body);
Tbody = stokes_DLP_mat(rout_body, rin_body, rout_body);
[Yres_dlp, Ures_dlp] = getPseudoFactors(Tbody * Sblock, 1e-10, 0);
Ures_dlp = Ures_dlp';

Tloop = stokes_DLP_mat(rout(1:M,:), rin(1:N,:), rout(1:M,:) - q(1,:));

opt = init_MFS();
opt.ellipsoid = 1;

data.Pmax = Pmax;
data.q = q;
data.R = R;
data.rin = rin;
data.rout = rout;
data.N = N;
data.M = M;
data.std.resistance.U = Ures_std;
data.std.resistance.Y = Yres_std;
data.std.mobility.U = Umob_std;
data.std.mobility.Y = Ymob_std;
data.std.mobility.L = Lstd;
data.dlp.resistance.U = Ures_dlp;
data.dlp.resistance.Y = Yres_dlp;
data.dlp.mobility.U = Umob_dlp;
data.dlp.mobility.Y = Ymob_dlp;
data.dlp.mobility.L = Ldlp;
data.Sblock_vec = Sblock;
data.Sblock_loop = Sblock;
data.Tblock_loop = Tloop;
data.TSblock_vec = Tbody * Sblock;
data.mu_standard = randn(3*M*Pmax, 1);
data.mu_dlp = randn(3*N*Pmax, 1);
data.lambda_stokes = randn(3*N*Pmax, 1);
data.vars_standard = benchmark_vars(opt, M, N, true, Sblock, []);
data.vars_dlp = benchmark_vars(opt, M, N, true, Sblock, Tbody * Sblock);
end

function vars = benchmark_vars(opt, M, N, is_ellipsoid, Sblock, TSblock)
vars = opt;
vars.fmm = 0;
vars.profile = 0;
vars.fmm_tol = 1e-8;
vars.M = M;
vars.N = N;
vars.ellipsoid = is_ellipsoid;
vars.Sblock = Sblock;
if ~isempty(TSblock)
    vars.TSblock = TSblock;
end
end

function case_data = slice_case_data(data, P)
case_data = data;
case_data.P = P;
case_data.R = data.R(1:P);
case_data.rin = data.rin(1:data.N*P,:);
case_data.rout = data.rout(1:data.M*P,:);
case_data.mu_standard = data.mu_standard(1:3*data.M*P);
case_data.mu_dlp = data.mu_dlp(1:3*data.N*P);
case_data.lambda_stokes = data.lambda_stokes(1:3*data.N*P);
case_data.vars_standard = data.vars_standard;
case_data.vars_standard.M = data.M;
case_data.vars_standard.N = data.N;
case_data.vars_dlp = data.vars_dlp;
case_data.vars_dlp.M = data.M;
case_data.vars_dlp.N = data.N;
end

function lambda_stokes = standard_precond_loop(mu, U, Y, L, vars, R, resistance_flag)
P = numel(R);
M = vars.M;
N = vars.N;
lambda_stokes = zeros(3*P*N, 1);

for i = 1:P
    rows_m = (i-1)*3*M + (1:3*M);
    rows_n = (i-1)*3*N + (1:3*N);

    if resistance_flag
        if vars.ellipsoid
            Ri = R{i};
            step0 = rotate_vector(mu(rows_m), Ri');
            step1 = U * step0;
            lambda_i = rotate_vector(Y * step1, Ri);
        else
            step1 = U * mu(rows_m);
            lambda_i = Y * step1;
        end
    else
        if vars.ellipsoid
            Ri = R{i};
            step0 = rotate_vector(mu(rows_m), Ri');
            tau_body = Y * (U * step0);
            lambda_i = rotate_vector(tau_body - L * tau_body, Ri);
        else
            tau_body = Y * (U * mu(rows_m));
            lambda_i = tau_body - L * tau_body;
        end
    end

    lambda_stokes(rows_n) = lambda_i;
end
end

function lambda_stokes = standard_precond_vectorized(mu, U, Y, L, vars, R, resistance_flag)
M = vars.M;
N = vars.N;
P = numel(R);

if vars.ellipsoid
    Rpages = cat(3, R{:});
    mu_body = rotate_vector_pages(mu, permute(Rpages, [2 1 3]), M);
    step1 = U * reshape(mu_body, 3*M, P);
    tau_body = Y * step1;

    if resistance_flag
        lambda_stokes = rotate_vector_pages(tau_body, Rpages, N);
    else
        lambda_stokes = rotate_vector_pages(tau_body - L * tau_body, Rpages, N);
    end

    lambda_stokes = lambda_stokes(:);
    return;
end

step1 = U * reshape(mu, 3*M, P);
tau_body = Y * step1;
if resistance_flag
    lambda_stokes = tau_body(:);
else
    lambda_stokes = reshape(tau_body - L * tau_body, [], 1);
end
end

function lambda_stokes = dlp_precond_loop(mu, U, Y, L, vars, R, resistance_flag)
P = numel(R);
N = vars.N;
lambda_stokes = zeros(3*P*N, 1);

for i = 1:P
    rows_n = (i-1)*3*N + (1:3*N);

    if resistance_flag
        if vars.ellipsoid
            Ri = R{i};
            step0 = rotate_vector(mu(rows_n), Ri');
            step1 = U * step0;
            lambda_i = rotate_vector(Y * step1, Ri);
        else
            step1 = U * mu(rows_n);
            lambda_i = Y * step1;
        end
    else
        if vars.ellipsoid
            Ri = R{i};
            step0 = rotate_vector(mu(rows_n), Ri');
            tau_body = Y * (U * step0);
            lambda_i = rotate_vector(tau_body - L * tau_body, Ri);
        else
            tau_body = Y * (U * mu(rows_n));
            lambda_i = tau_body - L * tau_body;
        end
    end

    lambda_stokes(rows_n) = lambda_i;
end
end

function lambda_stokes = dlp_precond_vectorized(mu, U, Y, L, vars, R, resistance_flag)
N = vars.N;
P = numel(R);

if vars.ellipsoid
    Rpages = cat(3, R{:});
    mu_body = rotate_vector_pages(mu, permute(Rpages, [2 1 3]), N);
    step1 = U * reshape(mu_body, 3*N, P);
    tau_body = Y * step1;

    if resistance_flag
        lambda_stokes = rotate_vector_pages(tau_body, Rpages, N);
    else
        lambda_stokes = rotate_vector_pages(tau_body - L * tau_body, Rpages, N);
    end

    lambda_stokes = lambda_stokes(:);
    return;
end

step1 = U * reshape(mu, 3*N, P);
tau_body = Y * step1;
if resistance_flag
    lambda_stokes = tau_body(:);
else
    lambda_stokes = reshape(tau_body - L * tau_body, [], 1);
end
end

function u_self = standard_self_loop(lambda_stokes, rin, rout, vars, ~)
P = size(rin, 1) / vars.N;
vars_local = vars;
vars_local.fmm = 0;
u_self = zeros(3*vars.M*P, 1);

for i = 1:P
    rows_n = (i-1)*vars.N + (1:vars.N);
    rows_m = (i-1)*vars.M + (1:vars.M);
    vec_n = (i-1)*3*vars.N + (1:3*vars.N);
    vec_m = (i-1)*3*vars.M + (1:3*vars.M);
    u_self(vec_m) = getStokesletFlow(lambda_stokes(vec_n), rin(rows_n,:), rout(rows_m,:), vars_local);
end
end

function u_self = standard_self_vectorized(lambda_stokes, vars, R, Sblock)
N = vars.N;
M = vars.M;
P = numel(R);

if vars.ellipsoid
    Rpages = cat(3, R{:});
    lambda_body = rotate_vector_pages(lambda_stokes, permute(Rpages, [2 1 3]), N);
    u_body = Sblock * reshape(lambda_body, 3*N, P);
    u_self = rotate_vector_pages(u_body, Rpages, M);
    u_self = u_self(:);
    return;
end

u_self = reshape(Sblock * reshape(lambda_stokes, 3*N, P), [], 1);
end

function u_self = dlp_self_loop(lambda_stokes, rin, rout, vars, Tblock, R, TSblock)
% For ellipsoids, operates in the body frame using the precomputed body-frame
% TSblock (= Tbody*Sblock) and per-particle rotations R.
P = size(rin, 1) / vars.N;
vars_local = vars;
vars_local.fmm = 0;
u_self = zeros(3*vars.N*P, 1);

for i = 1:P
    rows_n = (i-1)*vars.N + (1:vars.N);
    vec_n  = (i-1)*3*vars.N + (1:3*vars.N);

    if vars.ellipsoid
        Ri = R{i};
        lambda_body_i = rotate_vector(lambda_stokes(vec_n), Ri');
        u_body_i      = TSblock * lambda_body_i;
        u_self(vec_n) = rotate_vector(u_body_i, Ri);
    else
        rows_m = (i-1)*vars.M + (1:vars.M);
        u_local       = getStokesletFlow(lambda_stokes(vec_n), rin(rows_n,:), rout(rows_m,:), vars_local);
        u_self(vec_n) = Tblock * u_local;
    end
end
end

function u_self = dlp_self_vectorized(lambda_stokes, vars, R, TSblock)
N = vars.N;
P = numel(R);

if vars.ellipsoid
    Rpages = cat(3, R{:});
    lambda_body = rotate_vector_pages(lambda_stokes, permute(Rpages, [2 1 3]), N);
    u_body = TSblock * reshape(lambda_body, 3*N, P);
    u_self = rotate_vector_pages(u_body, Rpages, N);
    u_self = u_self(:);
    return;
end

u_self = reshape(TSblock * reshape(lambda_stokes, 3*N, P), [], 1);
end

function R = random_rotations(P)
R = cell(1, P);
Rz = @(t) [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
Ry = @(t) [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)];

for k = 1:P
    u = randn(3,1);
    u = u / norm(u);
    alpha = atan2(u(2), u(1));
    beta = acos(u(3));
    gamma = 2*pi*rand;
    R{k} = Rz(alpha) * Ry(beta) * Rz(gamma);
end
end

function rel_err = relative_error(a, b)
rel_err = norm(a - b, inf) / max(norm(a, inf), eps);
end

function recommendation = recommend_switch(rel_err, speedup, speedup_threshold, equiv_tol)
if rel_err <= equiv_tol && speedup >= speedup_threshold
    recommendation = 'switch';
else
    recommendation = 'keep loop';
end
end

function row = make_row(P, loop_time, vec_time, speedup, rel_err, recommendation)
row.P = P;
row.loop_time = loop_time;
row.vec_time = vec_time;
row.speedup = speedup;
row.rel_err = rel_err;
row.recommendation = recommendation;
end

function print_row(row)
fprintf('  P=%-3d | %-14.4e | %-14.4e | %6.2fx     | %.3e     | %s\n', ...
    row.P, row.loop_time, row.vec_time, row.speedup, row.rel_err, row.recommendation);
end

function print_branch_summary(branch_results, label)
fprintf('%s\n', label);

if isfield(branch_results, 'sphere') && isfield(branch_results.sphere, 'resistance')
    summarise_case(branch_results.sphere.resistance, '  sphere resistance');
    summarise_case(branch_results.sphere.mobility, '  sphere mobility');
    summarise_case(branch_results.ellipsoid.resistance, '  ellipsoid resistance');
    summarise_case(branch_results.ellipsoid.mobility, '  ellipsoid mobility');
else
    summarise_case(branch_results.sphere, '  sphere');
    summarise_case(branch_results.ellipsoid, '  ellipsoid');
end
end

function summarise_case(rows, label)
keep_count = sum(strcmp({rows.recommendation}, 'keep loop'));
if keep_count == 0
    result = 'switch';
else
    result = 'keep loop';
end
fprintf('%s: %s\n', label, result);
end
