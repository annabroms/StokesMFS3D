function [figureHandle,profiles] = bubble_plot_residuals_distance( ...
        report,validation,solution,opt,visibility)
%BUBBLE_PLOT_RESIDUALS_DISTANCE Locate boundary errors by surface arc distance.
% Each panel retains all azimuths as a point cloud and plots the maximum and
% median in logarithmic distance bins. The exact sampled maximum is marked.
if nargin==0
    self_test();
    figureHandle = [];
    profiles = [];
    return
end
if nargin<5, visibility = 'on'; end
distance = validation.edgeDistance/opt.radius;
masks = {validation.cap,validation.mobile,validation.cap,validation.mobile};
values = {report.point.normal,report.point.normal, ...
    report.point.capTangent,report.point.mobileTraction};
scales = [solution.Uref,solution.Uref,solution.Uref,solution.tractionScale];
titles = {'Cap: normal velocity','Mobile: normal velocity', ...
    'Cap: tangential velocity','Mobile: tangential traction'};
figureHandle = figure('Name','Bubble residuals versus edge distance', ...
    'Color','w','Visible',visibility,'Position',[60 60 1400 850]);
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact','Padding','compact');
profiles = repmat(struct('distanceOverA',[],'residual',[], ...
    'binDistanceOverA',[],'binMaximum',[],'binMedian',[], ...
    'peakDistanceOverA',nan,'peakResidual',nan,'peakAzimuth',nan),4,1);
for k = 1:4
    ax = nexttile(layout,k);
    hold(ax,'on'); grid(ax,'on');
    set(ax,'XScale','log','YScale','log');
    mask = masks{k} & distance>0;
    s = distance(mask);
    residual = values{k}(mask)/scales(k);
    xlabel(ax,'Surface distance to edge s/a = |\theta-\phi|');
    if k<4
        ylabel(ax,'Residual / U_{ref}');
    else
        ylabel(ax,'Residual / (\mu_{out} U_{ref}/a)');
    end
    if isempty(s)
        title(ax,[titles{k},': no points in region']);
        continue
    end
    edges = logspace(log10(min(s))-1e-10,log10(max(s))+1e-10,61);
    bins = discretize(s,edges);
    centers = sqrt(edges(1:end-1).*edges(2:end));
    maxima = nan(60,1);
    medians = nan(60,1);
    for j = 1:60
        inBin = bins==j;
        if any(inBin)
            maxima(j) = max(residual(inBin));
            medians(j) = median(residual(inBin));
        end
    end
    occupied = isfinite(maxima);
    centers = centers(occupied).';
    maxima = maxima(occupied);
    medians = medians(occupied);
    scatter(ax,s,max(residual,1e-16),7,[0.35 0.55 0.8],'filled', ...
        'MarkerFaceAlpha',0.18,'DisplayName','all sampled azimuths');
    plot(ax,centers,max(maxima,1e-16),'r-','LineWidth',1.4, ...
        'DisplayName','maximum per distance bin');
    plot(ax,centers,max(medians,1e-16),'k--','LineWidth',1, ...
        'DisplayName','median per distance bin');
    [peak,index] = max(residual);
    plot(ax,s(index),max(peak,1e-16),'rp','MarkerSize',10, ...
        'MarkerFaceColor','r','DisplayName','sampled maximum');
    depth = solution.grid.minimumSourceDistance/opt.radius;
    xline(ax,depth,':','d_{min}/a','HandleVisibility','off');
    xlim(ax,[min(s)/1.2,max(s)*1.2]);
    title(ax,sprintf('%s: max %.3g at s/a=%.3g',titles{k},peak,s(index)));
    legend(ax,'Location','best');
    profiles(k).distanceOverA = s;
    profiles(k).residual = residual;
    profiles(k).binDistanceOverA = centers;
    profiles(k).binMaximum = maxima;
    profiles(k).binMedian = medians;
    profiles(k).peakDistanceOverA = s(index);
    profiles(k).peakResidual = peak;
    points = validation.points(mask,:)-opt.centre;
    profiles(k).peakAzimuth = atan2(points(index,2),points(index,1));
end
sgtitle(layout,sprintf(['Boundary residuals versus surface distance: %d points; ', ...
    'dotted line marks closest source depth'],validation.count));
end

function self_test()
% Check arc distances, normalization, exact maxima, and empty endpoint panels.
o = bubble_options(struct('radius',2,'phi',pi/3));
theta = [o.phi-0.1;o.phi-0.01;o.phi+0.02;o.phi+0.2];
v.points = o.radius*[sin(theta),zeros(4,1),cos(theta)];
v.edgeDistance = o.radius*abs(theta-o.phi);
v.cap = theta<o.phi;
v.mobile = ~v.cap;
v.count = 4;
r.point.normal = [2;8;6;4];
r.point.capTangent = [10;2;0;0];
r.point.mobileTraction = [0;0;3;15];
s.Uref = 2;
s.tractionScale = 3;
s.grid.minimumSourceDistance = 0.05;
[f,p] = bubble_plot_residuals_distance(r,v,s,o,'off');
assert(max(abs([p.peakResidual]-[4 3 5 5]))<1e-14);
assert(max(abs([p.peakDistanceOverA]-[0.01 0.02 0.1 0.2]))<1e-14);
assert(numel(findall(f,'Type','Scatter'))==4);
close(f);
v.cap(:) = false;
v.mobile = ~v.cap;
[f,p] = bubble_plot_residuals_distance(r,v,s,o,'off');
assert(isempty(p(1).residual) && isempty(p(3).residual));
close(f);
fprintf('bubble_plot_residuals_distance self-test passed.\n');
end
