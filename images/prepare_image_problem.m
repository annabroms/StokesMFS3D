function prepared = prepare_image_problem(q,opt,modes)
%PREPARE_IMAGE_PROBLEM Build one image geometry and requested one-body SVDs.
%   prepared = prepare_image_problem(q,opt,["resistance","mobility"])
% Pass the result as opt.prepared to either image solver at exactly this
% geometry and with exactly these options. No global dense matrix is built.

if nargin < 3
    modes = ["resistance","mobility"];
end
modes = string(modes);
assert(all(ismember(modes,["resistance","mobility"])), ...
    'Unknown image solve mode.');
opt = image_solver_options(opt);
if isfield(opt,'prepared')
    opt = rmfield(opt,'prepared');
end
assert(opt.solverMode == "oneBody", ...
    'Prepared factors apply only to oneBody image solves.');
timer = tic;
grid = build_image_discretization(q,opt);
P = size(q,1);
prepared = struct('q',q,'options',opt,'grid',grid, ...
    'resistanceFactors',[],'mobilityFactors',[], ...
    'gridAssemblies',1,'resistanceFactorizations',0, ...
    'mobilityFactorizations',0,'setupSeconds',0);
if any(modes == "resistance")
    factors = repmat(struct('Y',[],'U',[],'rank',0, ...
        'singularValues',[],'N_total',0,'collocationCount',0),P,1);
    for k = 1:P
        factors(k) = one_body_image_preconditioner( ...
            grid.body(k),opt.svdCutoff);
    end
    prepared.resistanceFactors = factors;
    prepared.resistanceFactorizations = P;
end
if any(modes == "mobility")
    factors = repmat(struct('H',[],'C',[],'Y',[],'U',[], ...
        'rank',0,'singularValues',[],'sourceCount',0, ...
        'collocationCount',0),P,1);
    for k = 1:P
        factors(k) = one_body_image_mobility_preconditioner( ...
            grid.body(k),q(k,:),opt.svdCutoff);
    end
    prepared.mobilityFactors = factors;
    prepared.mobilityFactorizations = P;
end
prepared.setupSeconds = toc(timer);
end
