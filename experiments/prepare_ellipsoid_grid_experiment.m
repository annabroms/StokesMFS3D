function [pair,opt,pairBuildTime] = prepare_ellipsoid_grid_experiment( ...
        caseName,delta,semiAxes,gridControls,opt)
%PREPARE_ELLIPSOID_GRID_EXPERIMENT Shared Laplace/Stokes options and geometry.
% Caller installs repository paths.  opt overrides the editable controls;
% the reflection experiment supplies geometry only, without its figures.
% "random" reuses its exact-gap support construction for configuration sweeps.

assert(isstruct(opt) && isscalar(opt),'opt must be a scalar struct.');
if isfield(opt,'semiAxes')
    semiAxes = opt.semiAxes;
    opt = rmfield(opt,'semiAxes');
end
controlNames = fieldnames(gridControls);
for k = 1:numel(controlNames)
    name = controlNames{k};
    if ~isfield(opt,name)
        opt.(name) = gridControls.(name);
    end
end
if ~isscalar(caseName) || ~ismember(caseName,["tip","side","oblique","random"])
    error('prepare_ellipsoid_grid_experiment:BadCase', ...
        'caseName must be "tip", "side", "oblique", or "random".');
end
validateattributes(delta,{'numeric'},{'real','finite','scalar','positive'});
opt = lap3d_complete_grid_options(opt);
pairStart = tic;
pair = sep09_ellipsoid_reflection_accumulation_experiment( ...
    caseName,delta,false,semiAxes);
pairBuildTime = toc(pairStart);
fprintf('Selected semiaxes: [%g, %g, %g]\n',pair.E1.axes);
end
