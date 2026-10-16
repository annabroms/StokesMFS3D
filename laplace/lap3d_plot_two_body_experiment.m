function lap3d_plot_two_body_experiment(result)
%LAP3D_PLOT_TWO_BODY_EXPERIMENT Plot grids, residuals, and source strengths.

displayedBody = 1;
accumulationPoint = nearest_accumulation_point(result.geometry,displayedBody);
plot_comparison(result,displayedBody,accumulationPoint);
plot_source_distribution(result,displayedBody,accumulationPoint);
plot_pointwise_residuals(result,displayedBody,accumulationPoint);
plot_source_magnitudes(result,displayedBody,accumulationPoint);
if isfield(result,'sourceShapeComparisons')
    plot_source_shape_comparison(result,displayedBody,accumulationPoint);
end
end


function plot_comparison(result,body,accumulationPoint)
%PLOT_COMPARISON shows one particle's diagnostics and the full configuration.

figure('Name',sprintf('Laplace MFS comparison: %s',result.caseName));
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
plot_grid(ax,result.geometry,result.basic.grid,false,body,accumulationPoint);
title(ax,sprintf('Particle %d: basic grid',body))

ax = nexttile(layout);
plot_grid(ax,result.geometry,result.enhanced.grid,true,body,accumulationPoint);
title(ax,sprintf('Particle %d: enhanced grid',body))

ax = nexttile(layout);
testMask = result.basic.grid.testBody == body;
testPoints = result.basic.grid.testPoints(testMask,:);
contactPoint = result.geometry.p1;
if body == 2
    contactPoint = result.geometry.p2;
end
distanceFromContact = vecnorm(testPoints-contactPoint.',2,2);
semilogy(ax,distanceFromContact, ...
    abs(result.basic.solution.testError(testMask)),'.', ...
    'DisplayName','basic');
hold(ax,'on')
semilogy(ax,distanceFromContact, ...
    abs(result.enhanced.solution.testError(testMask)),'.', ...
    'DisplayName','enhanced');
xlabel(ax,'distance from closest point')
ylabel(ax,'absolute boundary error')
grid(ax,'on')
legend(ax,'Location','best')
title(ax,'Independent boundary check')

ax = nexttile(layout);
% Both bodies belong in this configuration overview.  All source and
% residual panels remain restricted to the displayed particle above.
pair = result.geometry;
hold(ax,'on')
particle1 = draw_ellipsoid(ax,pair.E1);
particle2 = draw_ellipsoid(ax,pair.E2);
particle1.FaceColor = [0.20,0.45,0.85];
particle2.FaceColor = [0.85,0.40,0.20];
particle1.FaceAlpha = 0.35;
particle2.FaceAlpha = 0.35;
plot3(ax,[pair.p1(1),pair.p2(1)],[pair.p1(2),pair.p2(2)], ...
    [pair.p1(3),pair.p2(3)],'k-o','MarkerFaceColor',[1,0.85,0], ...
    'MarkerSize',5,'LineWidth',1.2);
mark_accumulation_point(ax,nearest_accumulation_point(pair,1));
mark_accumulation_point(ax,nearest_accumulation_point(pair,2));
legend(ax,[particle1,particle2],{'particle 1','particle 2'},'Location','best');
axis(ax,'equal')
grid(ax,'on')
view(ax,3)
xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')
title(ax,'Particle configuration')

title(layout,sprintf('%s, gap %.3g',result.geometry.label,result.delta))
end


function plot_source_distribution(result,body,accumulationPoint)
%PLOT_SOURCE_DISTRIBUTION isolates basic and enhancing source locations.

figure('Name',sprintf('Laplace MFS source distribution: %s',result.caseName));
ax = axes;
hold(ax,'on')
basicMask = result.basic.grid.sourceBody == body;
basicSources = result.basic.grid.sources(basicMask,:);
scatter3(ax,basicSources(:,1),basicSources(:,2),basicSources(:,3), ...
    12,[0.10,0.35,0.85],'filled','DisplayName','basic proxy');

capSources = result.enhanced.grid.sourceCapBlocks{body};
scatter3(ax,capSources(:,1),capSources(:,2),capSources(:,3), ...
    22,[0.90,0.15,0.60],'filled','DisplayName','enhancing cap');
mark_accumulation_point(ax,accumulationPoint);

axis(ax,'equal')
grid(ax,'on')
view(ax,3)
xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')
legend(ax,'Location','best')
title(ax,sprintf('%s: particle %d source distribution', ...
    result.geometry.label,body))
end


function plot_pointwise_residuals(result,body,accumulationPoint)
%PLOT_POINTWISE_RESIDUALS colors collocation and test nodes by local error.

basicCollocationMask = result.basic.grid.collocationBody == body;
enhancedCollocationMask = result.enhanced.grid.collocationBody == body;
testMask = result.basic.grid.testBody == body;
residuals = { ...
    result.basic.solution.collocationError(basicCollocationMask), ...
    result.enhanced.solution.collocationError(enhancedCollocationMask), ...
    result.basic.solution.testError(testMask), ...
    result.enhanced.solution.testError(testMask)};
logResiduals = cell(size(residuals));
for k = 1:numel(residuals)
    logResiduals{k} = log_magnitudes(residuals{k});
end

figure('Name',sprintf('Laplace MFS pointwise residuals: %s',result.caseName));
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

plot_colored_nodes(nexttile(layout), ...
    result.basic.grid.collocation(basicCollocationMask,:), ...
    logResiduals{1},'Basic collocation residual', ...
    'log_{10}|u-u_{prescribed}|',accumulationPoint);
plot_colored_nodes(nexttile(layout), ...
    result.enhanced.grid.collocation(enhancedCollocationMask,:), ...
    logResiduals{2},'Enhanced collocation residual', ...
    'log_{10}|u-u_{prescribed}|',accumulationPoint);
plot_colored_nodes(nexttile(layout),result.basic.grid.testPoints(testMask,:), ...
    logResiduals{3},'Basic evaluation residual', ...
    'log_{10}|u-u_{prescribed}|',accumulationPoint);
plot_colored_nodes(nexttile(layout),result.enhanced.grid.testPoints(testMask,:), ...
    logResiduals{4},'Enhanced evaluation residual', ...
    'log_{10}|u-u_{prescribed}|',accumulationPoint);

title(layout,sprintf('%s: particle %d pointwise boundary residuals', ...
    result.geometry.label,body))
end


function plot_source_magnitudes(result,body,accumulationPoint)
%PLOT_SOURCE_MAGNITUDES colors every proxy point by its coefficient size.

grids = {result.basic.grid};
solutions = {result.basic.solution};
panelTitles = "Basic sources";
if isfield(result,'sourceShapeComparisons')
    methods = result.sourceShapeComparisons;
    for k = 1:numel(methods)
        grids{end+1} = methods(k).grid; %#ok<AGROW>
        solutions{end+1} = methods(k).solution; %#ok<AGROW>
        panelTitles(end+1) = methods(k).name; %#ok<AGROW>
    end
else
    grids{end+1} = result.enhanced.grid;
    solutions{end+1} = result.enhanced.solution;
    panelTitles(end+1) = "Enhanced sources";
end

figure('Name',sprintf('Laplace MFS source magnitudes: %s',result.caseName));
if numel(grids) == 2
    layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
else
    layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
end

for k = 1:numel(grids)
    sourceMask = grids{k}.sourceBody == body;
    solution = solutions{k};
    colorbarLabel = 'log_{10}|q_j|';
    if isfield(solution,'sourceCoefficientMax')
        nodeMagnitudes = solution.sourceCoefficientMax(sourceMask);
        if ~isempty(solution.dipoleCoefficients)
            colorbarLabel = 'log_{10} max(|q_j|,|d_j|)';
        end
    else
        nodeMagnitudes = abs(solution.coefficients(sourceMask));
    end
    logCoefficients = log_magnitudes(nodeMagnitudes);
    plot_colored_nodes(nexttile(layout),grids{k}.sources(sourceMask,:), ...
        logCoefficients,panelTitles(k),colorbarLabel, ...
        accumulationPoint);
end

title(layout,sprintf('Particle %d MFS coefficient magnitudes',body))
end


function plot_source_shape_comparison(result,body,accumulationPoint)
%PLOT_SOURCE_SHAPE_COMPARISON contrasts the spheroid and two sphere caps.

methods = result.sourceShapeComparisons;
basicMask = result.basic.grid.sourceBody == body;
basicSources = result.basic.grid.sources(basicMask,:);

figure('Name',sprintf('Laplace MFS source-shape comparison: %s', ...
    result.caseName));
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for k = 1:numel(methods)
    ax = nexttile(layout);
    hold(ax,'on')
    scatter3(ax,basicSources(:,1),basicSources(:,2),basicSources(:,3), ...
        9,[0.10,0.35,0.85],'filled');
    capSources = methods(k).grid.sourceCapBlocks{body};
    scatter3(ax,capSources(:,1),capSources(:,2),capSources(:,3), ...
        18,[0.90,0.15,0.60],'filled');
    mark_accumulation_point(ax,accumulationPoint);
    axis(ax,'equal')
    grid(ax,'on')
    view(ax,3)
    xlabel(ax,'x')
    ylabel(ax,'y')
    zlabel(ax,'z')
    title(ax,methods(k).name)
end

ax = nexttile(layout);
collocationRatio = arrayfun( ...
    @(m) m.comparisonToBasic.collocationMaxResidualRatio,methods);
testRatio = arrayfun( ...
    @(m) m.comparisonToBasic.testMaxResidualRatio,methods);
semilogy(ax,1:numel(methods),collocationRatio,'o-', ...
    'DisplayName','collocation nodes');
hold(ax,'on')
semilogy(ax,1:numel(methods),testRatio,'s-', ...
    'DisplayName','independent test nodes');
yline(ax,1,'k:','DisplayName','basic grid');
xticks(ax,1:numel(methods))
xticklabels(ax,{'spheroid','small sphere','large sphere'})
ylabel(ax,'maximum relative residual ratio')
grid(ax,'on')
legend(ax,'Location','best')
title(ax,'Independent boundary-error comparison')

title(layout,sprintf('%s: particle %d enhancing source surfaces', ...
    result.geometry.label,body))
end


function plot_colored_nodes(ax,points,colorValues,plotTitle, ...
        colorbarLabel,accumulationPoint)
%PLOT_COLORED_NODES draws a logarithmic point cloud with its own color scale.

scatter3(ax,points(:,1),points(:,2),points(:,3),16,colorValues,'filled');
hold(ax,'on')
mark_accumulation_point(ax,accumulationPoint);
axis(ax,'equal')
grid(ax,'on')
view(ax,3)
xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')
clim(ax,'auto');
cb = colorbar(ax);
cb.Label.String = colorbarLabel;
title(ax,plotTitle)
end


function logValues = log_magnitudes(values)
%LOG_MAGNITUDES supplies finite log colors for one independently scaled plot.

magnitudes = abs(values);
maximumMagnitude = max(magnitudes);
if maximumMagnitude == 0
    magnitudeFloor = realmin;
else
    magnitudeFloor = max(realmin,maximumMagnitude*1e-15);
end
logValues = log10(max(magnitudes,magnitudeFloor));
end


function plot_grid(ax,pair,mfsGrid,showEnhancement,body,accumulationPoint)
%PLOT_GRID draws one physical body and its source/target node families.

hold(ax,'on')
draw_ellipsoid(ax,pair.E1);

sourceMask = mfsGrid.sourceBody == body;
collocationMask = mfsGrid.collocationBody == body;
scatter3(ax,mfsGrid.sources(sourceMask,1),mfsGrid.sources(sourceMask,2), ...
    mfsGrid.sources(sourceMask,3), ...
    5,[0.8,0.15,0.15]);
scatter3(ax,mfsGrid.collocation(collocationMask,1), ...
    mfsGrid.collocation(collocationMask,2), ...
    mfsGrid.collocation(collocationMask,3),5,[0.1,0.3,0.85]);

if showEnhancement
    scatter3(ax,mfsGrid.sourceCapBlocks{body}(:,1), ...
        mfsGrid.sourceCapBlocks{body}(:,2), ...
        mfsGrid.sourceCapBlocks{body}(:,3), ...
        14,[0.85,0.1,0.65],'filled');
    scatter3(ax,mfsGrid.collocationPatchBlocks{body}(:,1), ...
        mfsGrid.collocationPatchBlocks{body}(:,2), ...
        mfsGrid.collocationPatchBlocks{body}(:,3), ...
        10,[0.1,0.55,0.1],'filled');
end

plot3(ax,pair.Pminus(1,:),pair.Pminus(2,:),pair.Pminus(3,:), ...
    '-','Color',[0.35,0.35,0.35],'LineWidth',1.2);
mark_accumulation_point(ax,accumulationPoint);
plot3(ax,pair.p1(1),pair.p1(2),pair.p1(3), ...
    'o','Color',[0.9,0.65,0],'MarkerFaceColor','y','MarkerSize',6);
axis(ax,'equal')
grid(ax,'on')
view(ax,3)
xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')
end


function point = nearest_accumulation_point(pair,body)
%NEAREST_ACCUMULATION_POINT selects the focal point nearest one surface.

if body == 1
    candidates = pair.Pminus;
    surfacePoint = pair.p1;
else
    candidates = pair.Pplus;
    surfacePoint = pair.p2;
end
[~,index] = min(vecnorm(candidates-surfacePoint,2,1));
point = candidates(:,index);
end


function mark_accumulation_point(ax,point)
%MARK_ACCUMULATION_POINT overlays the selected focal location in black.

ax.SortMethod = 'childorder';
scatter3(ax,point(1),point(2),point(3),50,[0,0,0],'filled', ...
    'Marker','d','MarkerEdgeColor',[1,1,1],'LineWidth',1.5, ...
    'DisplayName','nearest accumulation point');
end


function surfaceHandle = draw_ellipsoid(ax,E)
%DRAW_ELLIPSOID adds a translucent physical body to a supplied axes.

[x,y,z] = sphere(36);
points = E.Q*[E.axes(1)*x(:),E.axes(2)*y(:),E.axes(3)*z(:)].'+E.c;
surfaceHandle = surf(ax,reshape(points(1,:),size(x)),reshape(points(2,:),size(y)), ...
    reshape(points(3,:),size(z)),'EdgeColor','none','FaceAlpha',0.12, ...
    'FaceColor',[0.4,0.4,0.4]);
end
