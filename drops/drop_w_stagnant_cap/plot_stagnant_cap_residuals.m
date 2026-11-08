function figureHandle = plot_stagnant_cap_residuals( ...
    report,validation,solution,opt,visibility,colorLimits,label)
%PLOT_STAGNANT_CAP_RESIDUALS Surface maps of independent boundary residuals.
% Colours are log10 of pointwise residuals normalized by the velocity or
% traction scale used in stagnant_cap_validate. Each panel scales itself
% unless explicit colour limits are supplied.
if nargin == 0
    self_test();
    return
end
if nargin < 5 || isempty(visibility), visibility = 'on'; end
if nargin < 6, colorLimits = []; end
if nargin < 7, label = ""; end
opt = stagnant_cap_options(opt);
if ~isempty(colorLimits)
    assert(isequal(size(colorLimits),[1 2]) || ...
        isequal(size(colorLimits),[4 2]));
    assert(all(colorLimits(:,1)<colorLimits(:,2)));
end
assert(validation.count==numel(report.point.continuity));
assert(solution.velocityScale>0 && solution.tractionScale>0);

figureHandle = figure('Visible',visibility,'Color','w', ...
    'Position',[50 50 1450 840]);
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact', ...
    'Padding','compact');
names = {'(30) Velocity continuity; all interface', ...
    '(31) No penetration; all interface', ...
    '(32) Zero drop-relative tangential velocity; cap', ...
    '(33) Tangential traction continuity; mobile'};
formulas = {'||u_{out}-u_{in}||', ...
    '|n \cdot (u_{out}-U)|', ...
    '||\Pi(u_{out}-U)||', ...
    '||\Pi(\sigma_{out}-\sigma_{in})n||'};
fields = {'continuity','normal','capOuter','mobileTraction'};
masks = {true(validation.count,1),true(validation.count,1), ...
    validation.cap,validation.mobile};
scales = [repmat(solution.velocityScale,1,3),solution.tractionScale];
edgeAzimuth = linspace(0,2*pi,240);
edge = opt.centre + opt.radius*[sin(opt.phi)*cos(edgeAzimuth(:)), ...
    sin(opt.phi)*sin(edgeAzimuth(:)), ...
    cos(opt.phi)*ones(numel(edgeAzimuth),1)];
for k = 1:4
    ax = nexttile(layout,k);
    mask = masks{k};
    points = validation.points(mask,:);
    values = report.point.(fields{k})(mask)/scales(k);
    assert(all(isfinite(values) & values>=0));
    if isempty(colorLimits)
        positive = values(values>0);
        if isempty(positive)
            % A log of zero is undefined; this display floor is not data.
            floorLog = log10(eps);
        else
            floorLog = min(log10(positive))-1;
        end
        colors = floorLog*ones(size(values));
        colors(values>0) = log10(values(values>0));
        if isempty(colors)
            panelLimits = floorLog+[-0.5 0.5];
        else
            panelLimits = [min(colors),max(colors)];
            if panelLimits(1)==panelLimits(2)
                panelLimits = panelLimits+[-0.5 0.5];
            end
        end
    else
        if size(colorLimits,1)==1
            panelLimits = colorLimits;
        else
            panelLimits = colorLimits(k,:);
        end
        colors = log10(max(values,10^panelLimits(1)));
    end
    scatter3(ax,points(:,1),points(:,2),points(:,3),8,colors, ...
        'filled','MarkerFaceAlpha',0.9);
    hold(ax,'on');
    if opt.phi>0 && opt.phi<pi
        plot3(ax,edge(:,1),edge(:,2),edge(:,3),'k-', ...
            'LineWidth',0.8);
    end
    axis(ax,'equal'); axis(ax,'tight');
    view(ax,35,23); grid(ax,'on');
    colormap(ax,turbo);
    clim(ax,panelLimits);
    cb = colorbar(ax);
    if k<4
        cb.Label.String = 'log_{10}(residual / U_*)';
    else
        cb.Label.String = 'log_{10}(residual / (\mu_{out}U_*/a))';
    end
    title(ax,{sprintf('%s (%d points)',names{k},nnz(mask)), ...
        formulas{k}},'Interpreter','tex');
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
end
sgtitle(layout,sprintf(['%s Note boundary conditions (30)-(33): %s solve; ', ...
    'N=%d, M=%d; %d independent points; nearest edge %.3g a'], ...
    label,solution.kind,solution.grid.N,solution.grid.collocationCount, ...
    validation.count,validation.minimumEdgeDistance/opt.radius), ...
    'Interpreter','none');
end

function self_test()
%SELF_TEST Verify masks, normalization, and plotted logarithmic magnitudes.
o = stagnant_cap_options(struct('phi',pi/2));
v.points = [0 0 1;1 0 0;0 0 -1];
v.cap = [true;false;false];
v.mobile = ~v.cap;
v.count = 3;
v.minimumEdgeDistance = 0.1;
r.point.continuity = [0;0.1;1];
r.point.normal = [0.02;0.2;2];
r.point.capOuter = [0.04;0;0];
r.point.mobileTraction = [0;3;6];
s.kind = "speed";
s.grid.N = 50;
s.grid.collocationCount = 100;
s.velocityScale = 2;
s.tractionScale = 3;
f = plot_stagnant_cap_residuals(r,v,s,o,'off');
plots = findobj(f,'Type','Scatter');
assert(numel(plots)==4);
counts = sort(arrayfun(@(p) numel(p.CData),plots));
assert(isequal(counts(:),[1;2;3;3]));
allColors = vertcat(plots.CData);
assert(all(isfinite(allColors)) && any(abs(allColors)<1e-12));
for k = 1:numel(plots)
    range = clim(ancestor(plots(k),'axes'));
    assert(all(plots(k).CData>=range(1)-1e-12 & ...
        plots(k).CData<=range(2)+1e-12));
    if numel(unique(plots(k).CData))>1
        assert(abs(range(1)-min(plots(k).CData))<1e-12 && ...
            abs(range(2)-max(plots(k).CData))<1e-12);
    end
end
assert(any(abs(allColors-(log10(0.05)-1))<1e-12));
close(f);
limits = [-16 0;-15 1;-14 2;-13 3];
f = plot_stagnant_cap_residuals(r,v,s,o,'off',limits);
axesHandles = findall(f,'Type','Axes');
actualLimits = zeros(numel(axesHandles),2);
for k = 1:numel(axesHandles)
    actualLimits(k,:) = clim(axesHandles(k));
end
assert(isequal(sortrows(actualLimits),sortrows(limits)));
close(f);
fprintf('plot_stagnant_cap_residuals self-test passed.\n');
end
