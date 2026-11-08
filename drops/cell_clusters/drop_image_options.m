function opt = drop_image_options(opt,P,a)
%DROP_IMAGE_OPTIONS Complete the public drop image controls once.
if nargin == 0
    self_test(); opt = struct(); return
end
if nargin < 1 || isempty(opt), opt = struct(); end
if nargin < 2, P = 1; end
if nargin < 3, a = ones(P,1); end
if ~isfield(opt,'useRowScaling'), opt.useRowScaling = false; end
opt = image_solver_options(opt);
defaults = struct('outerProxyRatio',2,'outerImageOffsetFraction',0.01, ...
    'interiorViscosity',2, ...
    'mirror_ext',false,'velocityScale',1,'tractionScale',1, ...
    'tractionRowWeight',1,'scale_columns',false,'plot_grid',false);
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opt,names{k}) || isempty(opt.(names{k}))
        opt.(names{k}) = defaults.(names{k});
    end
end
if ~isfield(opt,'completionScale') || isempty(opt.completionScale)
    opt.completionScale = 1./(8*pi*a(:));
end
if isscalar(opt.interiorViscosity)
    opt.interiorViscosity = repmat(opt.interiorViscosity,P,1);
end
if isscalar(opt.completionScale)
    opt.completionScale = repmat(opt.completionScale,P,1);
end
opt.interiorViscosity = opt.interiorViscosity(:);
opt.completionScale = opt.completionScale(:);
validateattributes(opt.outerProxyRatio,{'numeric'},{'scalar','>',1});
validateattributes(opt.proxyRadius,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(opt.imageOffsetFraction,{'numeric'}, ...
    {'scalar','>=',0,'<',1});
validateattributes(opt.outerImageOffsetFraction,{'numeric'}, ...
    {'scalar','>',0,'<',1});
validateattributes(opt.mirror_ext,{'logical'},{'scalar'});
validateattributes(opt.interiorViscosity,{'numeric'}, ...
    {'vector','numel',P,'positive','finite'});
validateattributes(opt.completionScale,{'numeric'}, ...
    {'vector','numel',P,'positive','finite'});
validateattributes(opt.velocityScale,{'numeric'},{'scalar','positive'});
validateattributes(opt.tractionScale,{'numeric'},{'scalar','positive'});
validateattributes(opt.tractionRowWeight,{'numeric'},{'scalar','positive'});
validateattributes(opt.scale_columns,{'logical'},{'scalar'});
validateattributes(opt.plot_grid,{'numeric','logical'}, ...
    {'scalar','binary'});
opt.plot_grid = logical(opt.plot_grid);
end

function self_test()
% Verify per-drop completion and disabled row scaling defaults.
opt = drop_image_options(struct('N',30),2,[1;2]);
assert(all(opt.completionScale==1./(8*pi*[1;2])));
assert(~opt.useRowScaling && ~opt.mirror_ext);
assert(~opt.plot_grid);
assert(opt.outerImageOffsetFraction==0.01);
opt = drop_image_options(struct('plot_grid',1),1,1);
assert(opt.plot_grid);
opt = drop_image_options(struct('outerImageOffsetFraction',0.005),1,1);
assert(opt.outerImageOffsetFraction==0.005);
fprintf('drop_image_options self-test passed.\n');
end
