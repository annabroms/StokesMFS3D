function opt = lap3d_complete_grid_options(userOptions)
%LAP3D_COMPLETE_GRID_OPTIONS Fill and validate two-body grid controls.

defaults.coarseSourceResolution = [28,18];
defaults.basicCollocationResolution = []; % quadrature override, otherwise repository default
defaults.fineSourceResolution = [32,20];
defaults.basicGridDistribution = "ellipsoidquadrature";
defaults.basicSourceCount = 400;
defaults.basicCollocationCount = 1350;
defaults.enhancingSourceCount = 640;
defaults.enhancingLaplaceBasis = "SLP";
defaults.laplaceDipoleLengthScale = 1;
defaults.columnScaling = false;
defaults.surfaceElementPreconditioning = false;
defaults.sourceOffset = 0.125;
defaults.sourceChopOffset = 0; % signed depth beyond the basic proxy cut
defaults.sourceBeta = 0.3;
defaults.sourceWidthScale = 1;
defaults.enhancingSourceShape = "spheroid";
defaults.sourceSphereRadius = [];
defaults.enhancingSourceDistribution = "tensor";
defaults.enhancingSourceExtraAzimuthNodes = 0;
defaults.enhancingSourceMinAzimuthNodes = 8;
defaults.removeDuplicates = true; % source nodes only; collocation is unchanged
defaults.collocationPatchAngle = 0.5;
defaults.collocationPatchLayers = 18;
defaults.collocationExponent = 3;
defaults.collocationMinAzimuth = 8;
defaults.collocationAzimuthOversampling = 1;
defaults.basicSvdCutoff = 1e-13;
defaults.enhancedSvdCutoff = 1e-13;
defaults.testResolution = [64,40];
defaults.testPatchLayers = 0; % independent evaluation rings; 0 preserves old grid
defaults.testPatchAzimuthOversampling = 2;
defaults.visualise = true;
defaults.confirmBeforeSolve = true;

if nargin < 1 || isempty(userOptions)
    userOptions = struct();
end
if ~isstruct(userOptions) || ~isscalar(userOptions)
    error('lap3d_complete_grid_options:BadOptions', ...
        'opt must be a scalar struct.');
end

knownFields = fieldnames(defaults);
userFields = fieldnames(userOptions);
unknownFields = setdiff(userFields,knownFields);
if ~isempty(unknownFields)
    error('lap3d_complete_grid_options:UnknownOption', ...
        'Unknown option field "%s".',unknownFields{1});
end

opt = defaults;
for k = 1:numel(userFields)
    opt.(userFields{k}) = userOptions.(userFields{k});
end

booleanNames = {'removeDuplicates','confirmBeforeSolve', ...
    'columnScaling','surfaceElementPreconditioning'};
for k = 1:numel(booleanNames)
    name = booleanNames{k};
    validateattributes(opt.(name),{'logical','numeric'}, ...
        {'real','scalar','binary'},mfilename,name);
    opt.(name) = logical(opt.(name));
end
if ~isempty(opt.basicCollocationResolution)
    validateattributes(opt.basicCollocationResolution,{'numeric'}, ...
        {'real','finite','vector','numel',2,'integer','>=',4}, ...
        mfilename,'basicCollocationResolution');
end

opt.basicGridDistribution = lower(string(opt.basicGridDistribution));
opt.enhancingLaplaceBasis = upper(string(opt.enhancingLaplaceBasis));
if ~isscalar(opt.enhancingLaplaceBasis) || ...
        ~ismember(opt.enhancingLaplaceBasis,["SLP","SLP+DLP"])
    error('lap3d_complete_grid_options:BadLaplaceBasis', ...
        'enhancingLaplaceBasis must be "SLP" or "SLP+DLP".');
end
if ~isscalar(opt.basicGridDistribution) || ...
        ~ismember(opt.basicGridDistribution, ...
        ["ellipsoidquadrature","sphericaldesign"])
    error('lap3d_complete_grid_options:BadBasicDistribution', ...
        ['basicGridDistribution must be "ellipsoidquadrature" or ', ...
         '"sphericaldesign".']);
end
countNames = {'basicSourceCount','basicCollocationCount', ...
    'enhancingSourceCount'};
for k = 1:numel(countNames)
    name = countNames{k};
    value = opt.(name);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value < 3 || value ~= round(value)
        error('lap3d_complete_grid_options:BadDesignCount', ...
            '%s must be an integer of at least three.',name);
    end
end

resolutionNames = {'coarseSourceResolution','fineSourceResolution', ...
    'testResolution'};
for k = 1:numel(resolutionNames)
    name = resolutionNames{k};
    value = opt.(name);
    if ~isnumeric(value) || numel(value) ~= 2 || any(~isfinite(value)) || ...
            any(value < 4) || any(value ~= round(value))
        error('lap3d_complete_grid_options:BadResolution', ...
            '%s must contain two integers of at least four.',name);
    end
end

positiveScalarNames = {'sourceOffset','collocationPatchAngle', ...
    'collocationExponent','collocationAzimuthOversampling','sourceWidthScale', ...
    'testPatchAzimuthOversampling','laplaceDipoleLengthScale'};
for k = 1:numel(positiveScalarNames)
    name = positiveScalarNames{k};
    value = opt.(name);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value <= 0
        error('lap3d_complete_grid_options:BadPositiveScalar', ...
            '%s must be a positive finite scalar.',name);
    end
end

validateattributes(opt.sourceChopOffset,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'sourceChopOffset');
if ~isfinite(opt.sourceOffset+opt.sourceChopOffset) || ...
        opt.sourceOffset+opt.sourceChopOffset <= 0
    error('lap3d_complete_grid_options:BadChopDepth', ...
        'sourceOffset + sourceChopOffset must be positive and finite.');
end
validateattributes(opt.testPatchLayers,{'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative'},mfilename,'testPatchLayers');

if ~isscalar(opt.sourceBeta) || ~isfinite(opt.sourceBeta) || ...
        opt.sourceBeta <= 0 || opt.sourceBeta >= 1
    error('lap3d_complete_grid_options:BadSourceBeta', ...
        'sourceBeta must lie strictly between zero and one.');
end
opt.enhancingSourceShape = lower(string(opt.enhancingSourceShape));
if ~isscalar(opt.enhancingSourceShape) || ...
        ~ismember(opt.enhancingSourceShape,["spheroid","sphere"])
    error('lap3d_complete_grid_options:BadSourceShape', ...
        'enhancingSourceShape must be "spheroid" or "sphere".');
end
if opt.enhancingSourceShape == "sphere"
    if ~isnumeric(opt.sourceSphereRadius) || ...
            ~isscalar(opt.sourceSphereRadius) || ...
            ~isreal(opt.sourceSphereRadius) || ...
            ~isfinite(opt.sourceSphereRadius) || opt.sourceSphereRadius <= 0
        error('lap3d_complete_grid_options:BadSourceSphereRadius', ...
            'sourceSphereRadius must be positive for a spherical source surface.');
    end
elseif ~isempty(opt.sourceSphereRadius)
    error('lap3d_complete_grid_options:UnusedSourceSphereRadius', ...
        'sourceSphereRadius is used only when enhancingSourceShape is "sphere".');
end
opt.enhancingSourceDistribution = lower(string( ...
    opt.enhancingSourceDistribution));
if ~isscalar(opt.enhancingSourceDistribution) || ...
        ~ismember(opt.enhancingSourceDistribution, ...
        ["tensor","quasiuniform","sphericaldesign"])
    error('lap3d_complete_grid_options:BadSourceDistribution', ...
        ['enhancingSourceDistribution must be "tensor", ', ...
         '"quasiuniform", or "sphericaldesign".']);
end
if ~isscalar(opt.enhancingSourceExtraAzimuthNodes) || ...
        ~isfinite(opt.enhancingSourceExtraAzimuthNodes) || ...
        opt.enhancingSourceExtraAzimuthNodes < 0 || ...
        opt.enhancingSourceExtraAzimuthNodes ~= ...
        round(opt.enhancingSourceExtraAzimuthNodes)
    error('lap3d_complete_grid_options:BadSourceExtraAzimuthNodes', ...
        'enhancingSourceExtraAzimuthNodes must be a nonnegative integer.');
end
if ~isscalar(opt.enhancingSourceMinAzimuthNodes) || ...
        ~isfinite(opt.enhancingSourceMinAzimuthNodes) || ...
        opt.enhancingSourceMinAzimuthNodes < 2 || ...
        opt.enhancingSourceMinAzimuthNodes ~= ...
        round(opt.enhancingSourceMinAzimuthNodes)
    error('lap3d_complete_grid_options:BadSourceMinAzimuthNodes', ...
        'enhancingSourceMinAzimuthNodes must be an integer of at least two.');
end
if ~isscalar(opt.collocationPatchLayers) || ...
        opt.collocationPatchLayers < 1 || ...
        opt.collocationPatchLayers ~= round(opt.collocationPatchLayers)
    error('lap3d_complete_grid_options:BadPatchLayers', ...
        'collocationPatchLayers must be a positive integer.');
end
if ~isscalar(opt.collocationMinAzimuth) || ...
        opt.collocationMinAzimuth < 4 || ...
        opt.collocationMinAzimuth ~= round(opt.collocationMinAzimuth)
    error('lap3d_complete_grid_options:BadPatchAzimuth', ...
        'collocationMinAzimuth must be an integer of at least four.');
end
if opt.collocationPatchAngle >= pi
    error('lap3d_complete_grid_options:BadPatchAngle', ...
        'collocationPatchAngle must be smaller than pi.');
end
if ~isscalar(opt.basicSvdCutoff) || opt.basicSvdCutoff <= 0 || ...
        opt.basicSvdCutoff >= 1 || ~isfinite(opt.basicSvdCutoff)
    error('lap3d_complete_grid_options:BadBasicCutoff', ...
        'basicSvdCutoff must lie strictly between zero and one.');
end
if ~isscalar(opt.enhancedSvdCutoff) || opt.enhancedSvdCutoff <= 0 || ...
        opt.enhancedSvdCutoff >= 1 || ~isfinite(opt.enhancedSvdCutoff)
    error('lap3d_complete_grid_options:BadEnhancedCutoff', ...
        'enhancedSvdCutoff must lie strictly between zero and one.');
end
if ~isscalar(opt.visualise)
    error('lap3d_complete_grid_options:BadVisualise', ...
        'visualise must be a scalar logical value.');
end
opt.visualise = logical(opt.visualise);
end
