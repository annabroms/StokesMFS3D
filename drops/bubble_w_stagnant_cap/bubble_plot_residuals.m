function figureHandle = bubble_plot_residuals( ...
        report,validation,solution,opt,visibility,pointSetLabel,plotOptions)
%BUBBLE_PLOT_RESIDUALS Region-specific normalized exterior residual maps.
% Optional plotOptions selects axes, paired limits and conditions (5=full cap).
if nargin == 0
    bubble_self_test();
    figureHandle = [];
    return
end
if nargin<5, visibility = 'on'; end
if nargin<6, pointSetLabel = 'Independent validation'; end
if nargin<7, plotOptions = struct(); end
assert(all(ismember(fieldnames(plotOptions), ...
    {'axes','conditionIndices','colorLimits','logFloor'})));
if ~isfield(plotOptions,'axes'), plotOptions.axes = []; end
if ~isfield(plotOptions,'conditionIndices'), plotOptions.conditionIndices = 1:4; end
if ~isfield(plotOptions,'colorLimits'), plotOptions.colorLimits = []; end
if ~isfield(plotOptions,'logFloor'), plotOptions.logFloor = 1e-16; end
indices = plotOptions.conditionIndices;
assert(all(ismember(indices,1:5)) && ~isempty(indices));
assert(isscalar(plotOptions.logFloor) && plotOptions.logFloor>0);
externalAxes = ~isempty(plotOptions.axes);
if externalAxes
    assert(numel(plotOptions.axes)==numel(indices) && all(isgraphics(plotOptions.axes,'axes')));
    figureHandle = ancestor(plotOptions.axes(1),'figure');
else
    figureHandle = figure('Name',['Bubble boundary residuals: ',pointSetLabel], ...
        'Color','w','Visible',visibility,'Position',[50 50 1400 850]);
    layout = tiledlayout(figureHandle,ceil(numel(indices)/2),2, ...
        'TileSpacing','compact','Padding','compact');
end
if ~isempty(plotOptions.colorLimits)
    assert(isequal(size(plotOptions.colorLimits),[numel(indices),2]) && ...
        all(plotOptions.colorLimits(:,2)>plotOptions.colorLimits(:,1)));
end
masks = {validation.cap,validation.mobile, ...
    validation.cap,validation.mobile};
values = {report.point.normal,report.point.normal, ...
    report.point.capTangent,report.point.mobileTraction};
scales = [solution.Uref,solution.Uref,solution.Uref, ...
    solution.tractionScale];
titles = {'Cap: normal velocity','Mobile: normal velocity', ...
    'Cap: tangential velocity','Mobile: tangential traction'};
maxima = report.maxima;
if any(indices==5)
    masks{5} = validation.cap;
    values{5} = report.point.capVelocity;
    scales(5) = solution.Uref;
    titles{5} = 'Cap: full velocity';
    maxima(5) = report.maxCapVelocity;
end
psi = linspace(0,2*pi,240).';
edge = opt.centre+opt.radius*[sin(opt.phi)*cos(psi), ...
    sin(opt.phi)*sin(psi),cos(opt.phi)*ones(numel(psi),1)];
for panel = 1:numel(indices)
    k = indices(panel);
    if externalAxes
        ax = plotOptions.axes(panel);
    else
        ax = nexttile(layout,panel);
    end
    mask = masks{k};
    x = validation.points(mask,:);
    scaled = values{k}(mask)/scales(k);
    colors = log10(max(scaled,plotOptions.logFloor));
    scatter3(ax,x(:,1),x(:,2),x(:,3),9,colors,'filled');
    hold(ax,'on');
    if opt.phi>0 && opt.phi<pi
        plot3(ax,edge(:,1),edge(:,2),edge(:,3),'k-','LineWidth',1);
    end
    if ~isempty(plotOptions.colorLimits)
        clim(ax,plotOptions.colorLimits(panel,:));
    elseif ~isempty(colors)
        limits = [min(colors),max(colors)];
        if limits(1)==limits(2), limits = limits+[-0.5 0.5]; end
        clim(ax,limits);
    end
    colormap(ax,turbo);
    axis(ax,'equal'); grid(ax,'on'); view(ax,35,25);
    title(ax,sprintf('%s: max %.3g (%d points)', ...
        titles{k},maxima(k),nnz(mask)));
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    if k~=4
        ylabel(colorbar(ax),'log_{10}(residual/U_{ref})');
    else
        ylabel(colorbar(ax),'log_{10}(residual/(\mu_{out}U_{ref}/a))');
    end
end
if ~externalAxes
    sgtitle(layout,sprintf('%s: %d points; closest edge %.3g a', ...
        pointSetLabel,validation.count,validation.minimumEdgeDistance/opt.radius));
end
end

function bubble_self_test()
% Verify masks and all four plotted residual fields on a small sample.
o = bubble_options(struct('phi',pi/2));
v.points = [0 0 1;1 0 0;0 0 -1];
v.cap = [true;false;false];
v.mobile = ~v.cap;
v.count = 3;
v.minimumEdgeDistance = 0.1;
r.point.normal = [1;2;3];
r.point.capTangent = [4;0;0];
r.point.mobileTraction = [0;5;6];
r.maxima = [0.5,1.5,2,3];
s.Uref = 2;
s.tractionScale = 2;
f = bubble_plot_residuals(r,v,s,o,'off');
assert(numel(findall(f,'Type','Scatter'))==4);
close(f);
% A paired view must respect caller axes, common limits and the full cap mask.
r.point.capVelocity = [5;0;0];
r.maxCapVelocity = 2.5;
f = figure('Visible','off');
layout = tiledlayout(f,1,3);
ax = gobjects(1,3);
for k = 1:3, ax(k) = nexttile(layout,k); end
settings = struct('axes',ax,'conditionIndices',[5 2 4], ...
    'colorLimits',[-16 1;-16 1;-16 1]);
bubble_plot_residuals(r,v,s,o,'off','paired test',settings);
for k = 1:3, assert(isequal(clim(ax(k)),[-16 1])); end
points = findall(ax(1),'Type','Scatter');
assert(numel(points.XData)==nnz(v.cap));
close(f);
fprintf('bubble_plot_residuals self-test passed.\n');
end
