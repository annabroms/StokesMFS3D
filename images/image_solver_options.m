function opt = image_solver_options(opt)
%IMAGE_SOLVER_OPTIONS Complete and check shared image-solver controls.
%   image_solver_options() runs a small options regression check.

if nargin == 0
    self_test();
    opt = struct();
    return
end
if isempty(opt)
    opt = struct();
end
defaults = struct('proxyRadius',0.63,'N',686, ...
    'viscosity',1, ...
    'collocationRatio',1.2,'nearThreshold',0.15,'N_image',[], ...
    'imageOffsetFraction',0.05,'imageFamilies',"SRD",'useImages',true, ...
    'use_orth_dir_images',false, ...
    'capMode',"fibonacci",'capAngles',[pi/5,pi/60], ...
    'capUpsampling',[6,6],'capPointCounts',[], ...
    'capLayers',[18,18],'capExponent',[3,3], ...
    'capMinAzimuth',8,'rowWeightMode',"",'useRowScaling',true, ...
    'svdCutoff',5e-12,'gmresTolerance',1e-6,'maxIterations',800, ...
    'solverMode',"oneBody",'batchSize',128, ...
    'validationCount',1000,'validationGapCount',401);
names = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(opt,name) || isempty(opt.(name))
        opt.(name) = defaults.(name);
    end
end
opt.capMode = string(opt.capMode);
opt.rowWeightMode = string(opt.rowWeightMode);
opt.imageFamilies = string(opt.imageFamilies);
opt.solverMode = string(opt.solverMode);
if opt.rowWeightMode == ""
    if opt.capMode == "fibonacci"
        opt.rowWeightMode = "paper";
    else
        opt.rowWeightMode = "local";
    end
end
assert(ismember(opt.capMode,["fibonacci","exponential"]));
assert(ismember(opt.rowWeightMode,["paper","local"]));
assert(islogical(opt.useRowScaling) && isscalar(opt.useRowScaling), ...
    'useRowScaling must be a logical scalar.');
assert(isscalar(opt.solverMode) && ...
    ismember(opt.solverMode,["oneBody","denseSVD"]), ...
    'solverMode must be "oneBody" or "denseSVD".');
assert(opt.capMode ~= "exponential" || opt.rowWeightMode == "local", ...
    'Exponentially clustered caps require local merged-grid row areas.');
image_family_list(opt.imageFamilies);
assert(numel(opt.capAngles) == 2 && all(opt.capAngles > 0 & opt.capAngles < pi));
assert(numel(opt.capUpsampling) == 2 && all(opt.capUpsampling > 0));
nCaps = 2;
if opt.capMode == "exponential"
    nCaps = 1;
end
if ~isempty(opt.capPointCounts)
    validateattributes(opt.capPointCounts,{'numeric'}, ...
        {'vector','integer','>=',4});
    assert(numel(opt.capPointCounts) == nCaps || ...
        (nCaps == 1 && numel(opt.capPointCounts) == 2));
end
assert(numel(opt.capLayers) >= nCaps && all(opt.capLayers >= 1));
assert(numel(opt.capExponent) >= nCaps && all(opt.capExponent > 0));
validateattributes(opt.capMinAzimuth,{'numeric'}, ...
    {'scalar','integer','>=',3});
assert(islogical(opt.useImages) && isscalar(opt.useImages));
validateattributes(opt.use_orth_dir_images,{'logical'},{'scalar'});
validateattributes(opt.svdCutoff,{'numeric'},{'scalar','positive','<',1});
validateattributes(opt.viscosity,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(opt.nearThreshold,{'numeric'},{'scalar','positive'});
if ~isempty(opt.N_image)
    validateattributes(opt.N_image,{'numeric'}, ...
        {'scalar','integer','nonnegative'});
end
validateattributes(opt.batchSize,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opt.validationCount,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opt.validationGapCount,{'numeric'}, ...
    {'scalar','integer','>=',4});
end

function self_test()
%SELF_TEST checks defaults and cap-dependent row weights.
opt = image_solver_options(struct());
assert(opt.imageFamilies == "SRD" && opt.rowWeightMode == "paper");
opt = image_solver_options(struct('capMode',"exponential"));
assert(opt.rowWeightMode == "local");
assert(opt.useRowScaling);
opt = image_solver_options(struct('useRowScaling',false));
assert(~opt.useRowScaling && opt.rowWeightMode == "paper");
assert(~opt.use_orth_dir_images);
opt = image_solver_options(struct('use_orth_dir_images',true));
assert(opt.use_orth_dir_images);
opt = image_solver_options(struct('viscosity',1.7));
if opt.viscosity ~= 1.7
    error('image_solver_options:SelfTestFailed', ...
        'Viscosity option was not preserved.');
end
fprintf('image_solver_options self-test passed.\n');
end
