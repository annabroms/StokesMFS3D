function fig = plot_drop_pressure_sweep(file)
%PLOT_DROP_PRESSURE_SWEEP Display final sweep results from one MAT file.
% fig = plot_drop_pressure_sweep(file) returns motion/pressure, residual,
% and gap-pressure figures. It reads saved data without solves or file writes.
% plot_drop_pressure_sweep() runs a small deterministic plotting self-test.
if nargin == 0
    self_test(); fig = gobjects(0); return
end

% Load only the final tables and sampled pressure; source coefficients are unused.
d = load(file,'sweep','xGrid','zGrid','planePressure','exterior', ...
    'pressureLimit','pressureColours','q','radius','representativeGap');
sweep = d.sweep;
requestedN = unique(sweep.NRequested,'stable');
actualN = zeros(size(requestedN));
for k = 1:numel(requestedN)
    actualN(k) = sweep.NActual(find(sweep.NRequested==requestedN(k),1));
end
xGrid = d.xGrid; zGrid = d.zGrid; planePressure = d.planePressure;
exterior = d.exterior; pressureLimit = d.pressureLimit;
pressureColours = d.pressureColours; q = d.q; radius = d.radius;
representativeGap = d.representativeGap;
xExtent = max(abs(xGrid)); zExtent = max(abs(zGrid));
fig = gobjects(3,1);

% Compare translations, midpoint pressures, and independent surface residuals.
colours = lines(numel(requestedN));
fig(1) = figure('Name','Axial drop motion and midpoint pressure versus gap', ...
    'Color','w','Position',[100 100 1250 480]);
tiledlayout(1,2);
ax = nexttile; hold(ax,'on');
for resolution = 1:numel(requestedN)
    rows = sweep.NRequested==requestedN(resolution);
    plot(ax,sweep.Gap(rows),sweep.U1Toward(rows),'o-', ...
        'Color',colours(resolution,:),'DisplayName',sprintf('U1x, N=%d',actualN(resolution)));
    plot(ax,sweep.Gap(rows),sweep.U2Toward(rows),'s--', ...
        'Color',colours(resolution,:),'DisplayName',sprintf('-U2x, N=%d',actualN(resolution)));
end
yline(ax,0,'k:','HandleVisibility','off');
set(ax,'XScale','log','XGrid','on','YGrid','on');
xlabel(ax,'Gap delta'); ylabel(ax,'Speed toward the other drop');
title(ax,'Positive speeds mean approach'); legend(ax,'Location','best');
ax = nexttile; hold(ax,'on');
for resolution = 1:numel(requestedN)
    rows = sweep.NRequested==requestedN(resolution);
    loglog(ax,sweep.Gap(rows),sweep.PressureMagnitude(rows),'o-', ...
        'Color',colours(resolution,:),'DisplayName',sprintf('N=%d',actualN(resolution)));
end
set(ax,'XScale','log','YScale','log','XGrid','on','YGrid','on');
xlabel(ax,'Gap delta'); ylabel(ax,'|p_0(0,0,0)|');
title(ax,'Exterior midpoint pressure'); legend(ax,'Location','best');

fig(2) = figure('Name','Independent boundary residuals versus gap', ...
    'Color','w','Position',[100 100 1400 480]);
tiledlayout(1,3);
conditions = {'ContinuityRel','NormalRel','TractionRel'};
titles = {'Velocity continuity','Normal velocity','Tangential traction'};
for k = 1:3
    ax = nexttile; hold(ax,'on');
    for resolution = 1:numel(requestedN)
        rows = sweep.NRequested==requestedN(resolution);
        loglog(ax,sweep.Gap(rows),sweep.(conditions{k})(rows),'o-', ...
            'Color',colours(resolution,:),'DisplayName',sprintf('N=%d',actualN(resolution)));
    end
    set(ax,'XScale','log','YScale','log','XGrid','on','YGrid','on');
    xlabel(ax,'Gap delta'); ylabel(ax,'Maximum scaled residual');
    title(ax,titles{k}); legend(ax,'Location','best');
end

% Preserve the saved zoom, common colour limits, and white zero level.
fig(3) = figure('Name','Gap pressure with images','Color','w', ...
    'Position',[100 100 700 800]);
tiledlayout(1,numel(requestedN),'TileSpacing','compact');
theta = linspace(0,2*pi,1201);
for resolution = 1:numel(requestedN)
    ax = nexttile;
    h = imagesc(ax,xGrid,zGrid,reshape(planePressure(:,resolution),[numel(zGrid),numel(xGrid)]));
    set(h,'AlphaData',reshape(exterior,[numel(zGrid),numel(xGrid)]));
    set(ax,'YDir','normal','Color',[0.82 0.82 0.82]);
    hold(ax,'on');
    for i = 1:2
        plot(ax,q(i,1)+radius*cos(theta),radius*sin(theta),'k-','LineWidth',1.2);
    end
    plot(ax,0,0,'ko','MarkerSize',3);
    axis(ax,'equal'); xlim(ax,[-xExtent xExtent]); ylim(ax,[-zExtent zExtent]);
    clim(ax,[-pressureLimit pressureLimit]); colormap(ax,pressureColours);
    cb = colorbar(ax); cb.Label.String = 'Exterior pressure p_0';
    cb.Ticks = linspace(-pressureLimit,pressureLimit,5);
    cb.TickLabels = compose('%.3g',cb.Ticks);
    xlabel(ax,'x'); ylabel(ax,'z');
    title(ax,sprintf('N=%d, delta=%.3g',actualN(resolution),representativeGap));
end
sgtitle('Gap pressure: white = zero; grey = drop interior');
end

function self_test()
% Check saved-data plotting, pressure layout, masking, and colour limits.
sweep = table([20;20;40;40],[18;18;38;38],[0.1;0.01;0.1;0.01], ...
    [1;2;1;2],[1;2;1;2],[2;4;2;4], ...
    0.01*ones(4,1),0.02*ones(4,1),0.03*ones(4,1), ...
    'VariableNames',{'NRequested','NActual','Gap','U1Toward','U2Toward', ...
    'PressureMagnitude','ContinuityRel','NormalRel','TractionRel'});
xGrid = [-0.06 0 0.06]; zGrid = linspace(-0.3,0.3,5);
planePressure = [(1:15).',-(1:15).'];
exterior = true(15,1); exterior([1 15]) = false;
planePressure(~exterior,:) = NaN;
pressureLimit = 15; pressureColours = [0 0 1;1 1 1;1 0 0];
q = [-1.05 0 0;1.05 0 0]; radius = 1; representativeGap = 0.1;
file = [tempname '.mat'];
save(file,'sweep','xGrid','zGrid','planePressure','exterior', ...
    'pressureLimit','pressureColours','q','radius','representativeGap');
fileCleanup = onCleanup(@() delete(file));
visible = get(groot,'defaultFigureVisible');
visibilityCleanup = onCleanup(@() set(groot,'defaultFigureVisible',visible));
set(groot,'defaultFigureVisible','off');
fig = plot_drop_pressure_sweep(file);
figureCleanup = onCleanup(@() close(fig));
if numel(fig)~=3 || numel(findobj(fig,'Type','axes'))~=7
    error('plot_drop_pressure_sweep:SelfTest','Unexpected plot count.');
end
h = findobj(fig(3),'Type','image');
for k = 1:2
    ax = h(k).Parent;
    col = 1+(h(k).CData(3,2)<0);
    if ~isequaln(h(k).CData,reshape(planePressure(:,col),5,3)) || ...
            ~isequal(h(k).AlphaData,reshape(exterior,5,3)) || ...
            ~isequal(ax.CLim,[-pressureLimit pressureLimit]) || ...
            ~isequal(ax.Colormap,pressureColours)
        error('plot_drop_pressure_sweep:SelfTest','Saved pressure data changed.');
    end
end
fprintf('plot_drop_pressure_sweep self-test passed.\n');
end
