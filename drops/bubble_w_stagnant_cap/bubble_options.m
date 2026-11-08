function opt = bubble_options(user)
%BUBBLE_OPTIONS Complete controls for the exterior-only stagnant-cap bubble.
if nargin == 0
    bubble_self_test();
    opt = [];
    return
end
if isempty(user), user = struct(); end
assert(isstruct(user) && isscalar(user));
if isfield(user,'geometryOptions')
    user = rmfield(user,'geometryOptions');
end
allowed = {'radius','centre','muOut','phi','innerProxyRatio', ...
    'sourceRequest','collocationRatio','ringLevels','ringAzimuthCount', ...
    'useExtraCollocation','edgeBandAngle','svdTolerance','batchSize', ...
    'validationBaseRequest','validationMultiplier','gravity', ...
    'densityDifference','completionScale','n_level', ...
    'sourceClusterRatio','sourceClustering','sourceClusterSigma', ...
    'innerSourceStartDistance','n_az_source', ...
    'n_az_colloc','collocAngularExtent','collocClusterRatio', ...
    'collocClustering','collocClusterSigma', ...
    'refinementFamilies','derivative_direction','use_orth_dir_images', ...
    'useColumnScaling','useAreaWeights','useSourceRefinement', ...
    'sourceRefinementLevels','sourceRefinementRadialFactors', ...
    'sourceRefinementAzimuthFactors','conditionWeights', ...
    'useConditionCountScaling','maxSvdMemoryGB'};
assert(all(ismember(fieldnames(user),allowed)), ...
    'bubble_options:UnknownOption','Unknown bubble option.');
if isfield(user,'refinementFamilies')
    families = image_family_list(user.refinementFamilies);
else
    families = "S";
end
if isfield(user,'derivative_direction')
    derivativeDirection = string(user.derivative_direction);
else
    derivativeDirection = "none";
end
assert(isscalar(derivativeDirection) && ...
    ismember(derivativeDirection,["none","radial","cone_normal"]), ...
    'bubble_options:BadDerivativeDirection', ...
    'Use none, radial, or cone_normal for derivative_direction.');
if isfield(user,'use_orth_dir_images')
    useOrthogonal = user.use_orth_dir_images;
else
    useOrthogonal = false;
end
assert(islogical(useOrthogonal) && isscalar(useOrthogonal));
% The old options are used only to reuse geometry and validation grids.
geometryInput = user;
solverOptions = struct('useColumnScaling',false,'useAreaWeights',false, ...
    'useConditionCountScaling',false);
solverNames = fieldnames(solverOptions);
for k = 1:numel(solverNames)
    name = solverNames{k};
    if isfield(user,name)
        solverOptions.(name) = user.(name);
        geometryInput = rmfield(geometryInput,name);
    end
    assert(islogical(solverOptions.(name)) && isscalar(solverOptions.(name)), ...
        'bubble_options:BadSolverFlag','%s must be a logical scalar.',name);
end
conditionWeights = [1 1 1];
if isfield(user,'conditionWeights')
    conditionWeights = user.conditionWeights;
    geometryInput = rmfield(geometryInput,'conditionWeights');
end
validateattributes(conditionWeights,{'numeric'}, ...
    {'vector','numel',3,'real','finite','positive'});
maxSvdMemoryGB = inf;
if isfield(user,'maxSvdMemoryGB')
    maxSvdMemoryGB = user.maxSvdMemoryGB;
    geometryInput = rmfield(geometryInput,'maxSvdMemoryGB');
end
validateattributes(maxSvdMemoryGB,{'numeric'},{'scalar','real','positive'});
% These bubble-only controls add collocation patches without moving sources.
refinement = struct('useSourceRefinement',false,'sourceRefinementLevels',1, ...
    'sourceRefinementRadialFactors',[0.25 0.5 1 2], ...
    'sourceRefinementAzimuthFactors',[-2 -1 0 1 2]);
names = fieldnames(refinement);
for k = 1:numel(names)
    name = names{k};
    if isfield(user,name)
        refinement.(name) = user.(name);
        geometryInput = rmfield(geometryInput,name);
    end
end
assert(islogical(refinement.useSourceRefinement) && ...
    isscalar(refinement.useSourceRefinement), ...
    'bubble_options:BadRefinementFlag','useSourceRefinement must be logical.');
validateattributes(refinement.sourceRefinementLevels,{'numeric'}, ...
    {'scalar','real','finite','integer','positive'});
validateattributes(refinement.sourceRefinementRadialFactors,{'numeric'}, ...
    {'vector','nonempty','real','finite','positive'});
validateattributes(refinement.sourceRefinementAzimuthFactors,{'numeric'}, ...
    {'vector','nonempty','real','finite'});
refinement.sourceRefinementRadialFactors = ...
    unique(refinement.sourceRefinementRadialFactors(:).');
refinement.sourceRefinementAzimuthFactors = ...
    unique(refinement.sourceRefinementAzimuthFactors(:).');
geometryInput.muIn = 1;
geometryInput.outerProxyRatio = 2;
geometryInput.refinementFamilies = "S";
geometryInput.derivative_direction = "none";
if ~isfield(geometryInput,'n_level')
    geometryInput.n_level = 0;
end
if geometryInput.n_level==0
    assert(~refinement.useSourceRefinement, ...
        'bubble_options:NoSourcesToRefine','Source patches require n_level>0.');
    assert(isequal(families,"S") && derivativeDirection=="none", ...
        'bubble_options:NoConeForFamilies', ...
        'Extra source families require positive n_level.');
end
if geometryInput.n_level>0
    if ~isfield(geometryInput,'ringLevels')
        geometryInput.ringLevels = 2*geometryInput.n_level;
    end
    assert(geometryInput.ringLevels>=2*geometryInput.n_level, ...
        'bubble_options:TooFewRings', ...
        'Cone geometry needs at least two collocation rings per source level.');
    % The shared refined builder uses a factor; ringLevels is the bubble API.
    geometryInput.colloc_factor = ...
        geometryInput.ringLevels/geometryInput.n_level;
end
if ~isfield(geometryInput,'collocClusterRatio') && geometryInput.n_level==0
    % Preserve the uniform-grid ring spacing: 0.5 per level by default.
    geometryInput.collocClusterRatio = 0.5;
end
geometryOptions = stagnant_cap_options(geometryInput);
removed = {'muIn','outerProxyRatio','outerSourceStartDistance', ...
    'refinementFamilies','derivative_direction','use_orth_dir_images', ...
    'plotThetaCount','plotAzimuthCount','colloc_factor'};
opt = rmfield(geometryOptions,removed);
opt.refinementFamilies = families;
opt.derivative_direction = derivativeDirection;
opt.use_orth_dir_images = useOrthogonal;
opt.useColumnScaling = solverOptions.useColumnScaling;
opt.useAreaWeights = solverOptions.useAreaWeights;
opt.useConditionCountScaling = solverOptions.useConditionCountScaling;
opt.conditionWeights = conditionWeights(:).';
opt.maxSvdMemoryGB = maxSvdMemoryGB;
for k = 1:numel(names)
    opt.(names{k}) = refinement.(names{k});
end
opt.geometryOptions = geometryOptions;
end

function bubble_self_test()
% Check bubble-only options and the geometry adapter's fixed source family.
o = bubble_options(struct('phi',pi/3,'n_level',5));
assert(o.n_level==5 && o.ringLevels==10 && o.muOut==1);
assert(~isfield(o,'colloc_factor'));
assert(~isfield(o,'muIn') && o.geometryOptions.refinementFamilies=="S");
assert(o.geometryOptions.derivative_direction=="none");
assert(bubble_options(o).n_level==5);
assert(~o.useColumnScaling && ~o.useAreaWeights);
scaled = bubble_options(struct('useColumnScaling',true,'useAreaWeights',true));
scaled = bubble_options(scaled);
assert(scaled.useColumnScaling && scaled.useAreaWeights && ...
    ~isfield(scaled.geometryOptions,'useAreaWeights') && ...
    ~isfield(scaled.geometryOptions,'useColumnScaling'));
for name = ["useColumnScaling","useAreaWeights"]
    try
        bubble_options(struct(name,1));
        error('bubble_options:SelfTest','Numeric solver flag was accepted.');
    catch ME
        assert(strcmp(ME.identifier,'bubble_options:BadSolverFlag'));
    end
end
mixed = bubble_options(struct('n_level',1, ...
    'refinementFamilies',["S","T_t","T_n","T_n2"], ...
    'derivative_direction',"radial",'use_orth_dir_images',true));
assert(numel(mixed.refinementFamilies)==4 && ...
    mixed.derivative_direction=="radial" && ...
    mixed.use_orth_dir_images && ...
    isequal(mixed.geometryOptions.refinementFamilies,"S"));
customCone = bubble_options(struct('n_level',5,'ringLevels',12));
assert(customCone.ringLevels==12 && ...
    abs(customCone.geometryOptions.colloc_factor-12/5)<eps);
uniform = bubble_options(struct('n_level',0));
assert(uniform.collocClusterRatio==0.5);
custom = bubble_options(struct('n_level',0,'collocClusterRatio',0.3));
assert(custom.collocClusterRatio==0.3);
try
    bubble_options(struct('muIn',0));
    error('bubble_options:SelfTest','Interior viscosity was accepted.');
catch ME
    assert(~strcmp(ME.identifier,'bubble_options:SelfTest'));
end
try
    bubble_options(struct('n_level',0,'refinementFamilies',"R"));
    error('bubble_options:SelfTest','A cone family without cone nodes was accepted.');
catch ME
    assert(~strcmp(ME.identifier,'bubble_options:SelfTest'));
end
fprintf('bubble_options self-test passed.\n');
end
