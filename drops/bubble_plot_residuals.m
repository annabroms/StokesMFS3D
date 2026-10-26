function figureHandle = bubble_plot_residuals(report,validation,solution,opt,visibility)
%BUBBLE_PLOT_RESIDUALS Four normalized exterior boundary residual maps.
if nargin == 0
    bubble_self_test();
    figureHandle = [];
    return
end
if nargin<5, visibility = 'on'; end
figureHandle = figure('Name','Bubble boundary residuals','Color','w', ...
    'Visible',visibility,'Position',[50 50 1400 850]);
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact', ...
    'Padding','compact');
masks = {validation.cap,validation.mobile, ...
    validation.cap,validation.mobile};
values = {report.point.normal,report.point.normal, ...
    report.point.capTangent,report.point.mobileTraction};
scales = [solution.Uref,solution.Uref,solution.Uref, ...
    solution.tractionScale];
titles = {'Cap: normal velocity','Mobile: normal velocity', ...
    'Cap: tangential velocity','Mobile: tangential traction'};
psi = linspace(0,2*pi,240).';
edge = opt.centre+opt.radius*[sin(opt.phi)*cos(psi), ...
    sin(opt.phi)*sin(psi),cos(opt.phi)*ones(numel(psi),1)];
for k = 1:4
    ax = nexttile(layout,k);
    mask = masks{k};
    x = validation.points(mask,:);
    scaled = values{k}(mask)/scales(k);
    colors = log10(max(scaled,1e-16));
    scatter3(ax,x(:,1),x(:,2),x(:,3),9,colors,'filled');
    hold(ax,'on');
    if opt.phi>0 && opt.phi<pi
        plot3(ax,edge(:,1),edge(:,2),edge(:,3),'k-','LineWidth',1);
    end
    if ~isempty(colors)
        limits = [min(colors),max(colors)];
        if limits(1)==limits(2), limits = limits+[-0.5 0.5]; end
        clim(ax,limits);
    end
    colormap(ax,turbo);
    axis(ax,'equal'); grid(ax,'on'); view(ax,35,25);
    title(ax,sprintf('%s: max %.3g (%d points)', ...
        titles{k},report.maxima(k),nnz(mask)));
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    if k<4
        ylabel(colorbar(ax),'log_{10}(residual/U_{ref})');
    else
        ylabel(colorbar(ax),'log_{10}(residual/(\mu_{out}U_{ref}/a))');
    end
end
sgtitle(layout,sprintf('Independent validation: %d points; closest edge %.3g a', ...
    validation.count,validation.minimumEdgeDistance/opt.radius));
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
fprintf('bubble_plot_residuals self-test passed.\n');
end
