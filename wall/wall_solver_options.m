function opt = wall_solver_options(input)
%WALL_SOLVER_OPTIONS Complete and validate options once at the wall entry point.
if nargin == 0, self_test(); opt = []; return; end
opt = struct('radius',1,'mu',1,'wall_height',0,'kernel',"blake", ...
    'backend',"fmm",'preconditioner',"wall",'solverMode',"oneBody", ...
    'N',700,'proxyRadius',0.63,'collocationRatio',1.2, ...
    'svdCutoff',1e-12,'gmresTolerance',1e-9,'fmmTolerance',1e-12, ...
    'maxIterations',400,'validationResolution',[46,56], ...
    'tractionResolution',[64,32],'computeResidual',true, ...
    'preparedProblem',[]);
assert(isstruct(input) && isscalar(input),'Options must be a scalar struct.');
names = fieldnames(input);
for j = 1:numel(names)
    assert(isfield(opt,names{j}),'wall:UnknownOption','Unknown option: %s',names{j});
    opt.(names{j}) = input.(names{j});
end
validateattributes(opt.radius,{'double'},{'real','finite','vector','positive'});
assert(all(opt.radius == opt.radius(1)), ...
    'wall:UnequalRadii','Only equal spheres are supported.');
validateattributes(opt.mu,{'double'},{'real','finite','scalar','positive'});
validateattributes(opt.wall_height,{'double'},{'real','finite','scalar'});
validateattributes(opt.N,{'double'},{'scalar','integer','>=',6});
validateattributes(opt.proxyRadius,{'double'},{'scalar','>',0,'<',1});
validateattributes(opt.collocationRatio,{'double'},{'scalar','>',1,'finite'});
validateattributes(opt.svdCutoff,{'double'},{'scalar','>',0,'<',1});
validateattributes(opt.gmresTolerance,{'double'},{'scalar','>',0,'<',1});
validateattributes(opt.fmmTolerance,{'double'},{'scalar','>',0,'<',1});
validateattributes(opt.maxIterations,{'double'},{'scalar','integer','positive'});
validateattributes(opt.validationResolution,{'double'},{'numel',2,'integer','>=',4});
validateattributes(opt.tractionResolution,{'double'},{'numel',2,'integer','>=',4});
validateattributes(opt.computeResidual,{'logical'},{'scalar'});
opt.kernel = string(validatestring(opt.kernel,{'blake'}));
opt.backend = string(validatestring(opt.backend,{'blake','ggv-direct','fmm'}));
opt.preconditioner = string(validatestring(opt.preconditioner,{'wall','freeSpace'}));
opt.solverMode = string(validatestring(opt.solverMode,{'oneBody','denseSVD'}));
end

function self_test()
% Check documented defaults and reject misspelled options.
o = wall_solver_options(struct());
assert(o.svdCutoff == 1e-12 && o.proxyRadius == 0.63);
failed = false;
try
    wall_solver_options(struct('typo',1));
catch
    failed = true;
end
assert(failed);
fprintf('wall_solver_options: PASS; Rp/a=%.2f, SVD cutoff=%.1e; unknown option rejected.\n', ...
    o.proxyRadius,o.svdCutoff);
end
