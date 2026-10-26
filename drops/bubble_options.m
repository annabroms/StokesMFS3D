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
    'sourceClusterRatio','innerSourceStartDistance','n_az_source', ...
    'n_az_colloc','collocAngularExtent','collocClusterRatio', ...
    'refinementFamilies','derivative_direction','use_orth_dir_images'};
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
geometryInput.muIn = 1;
geometryInput.outerProxyRatio = 2;
geometryInput.refinementFamilies = "S";
geometryInput.derivative_direction = "none";
if ~isfield(geometryInput,'n_level')
    geometryInput.n_level = 0;
end
if geometryInput.n_level==0
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
