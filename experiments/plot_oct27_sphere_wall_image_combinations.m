function [resultFigure,nodeFigure] = plot_oct27_sphere_wall_image_combinations(savedResults,outputStem)
%PLOT_OCT27_SPHERE_WALL_IMAGE_COMBINATIONS Replot saved metrics and nodes.
%   plot_oct27_sphere_wall_image_combinations('visuals/oct27_sphere_wall_image_combinations.mat')
%   [f,g] = plot_oct27_sphere_wall_image_combinations(filename,outputStem)
%   An optional outputStem exports both figures as FIG, PNG and PDF.
%   savedResults may also be a struct with the same saved arrays.
%   No kernels, factorizations or velocity evaluations are called.
%   A no-argument call runs a small synthetic-data plotting regression.
%
% Both result panels use ImageScalars: the scalar coefficients on image
% nodes, excluding the fixed 3*N standard proxy coefficients. Different
% families therefore have different abscissae for the same image-point count.
% Geometry shows the largest saved image count, shared by all combinations.

if nargin == 0
    self_test(); resultFigure = []; nodeFigure = []; return
end
if nargin < 2, outputStem = ''; end
if ischar(savedResults) || isstring(savedResults)
    % Load only plotting data, not the saved coefficient and basis arrays.
    data = load(savedResults,'parameters','results','centre', ...
        'proxyPoints','collocationPoints','linePoints');
else
    data = savedResults;
end
parameters = data.parameters;
results = data.results;
combinationNames = unique(results.Combination,'stable');
motion = parameters.rigidMotion;
radius = parameters.radius;
wallHeight = parameters.wallHeight;
centre = data.centre;
[~,lineIndex] = max(cellfun(@(points) size(points,1),data.linePoints));
imagePoints = data.linePoints{lineIndex};

%% ACCURACY AND SOURCE STRENGTH VERSUS SCALAR IMAGE UNKNOWNS
resultFigure = figure('Name','Oct 27 sphere-wall image combinations', ...
    'Position',[60,80,1500,760]);
layout = tiledlayout(resultFigure,1,2,'TileSpacing','compact','Padding','compact');
residualAxes = nexttile(layout); hold(residualAxes,'on');
strengthAxes = nexttile(layout); hold(strengthAxes,'on');
colors = lines(5);
styles = {'-','--',':'};
markers = {'o','s','^','d','v','>','<'};
for k = 1:numel(combinationNames)
    curve = results(results.Combination==combinationNames(k),:);
    curve = sortrows(curve,'ImageScalars');
    colorIndex = mod(k-1,5)+1;
    styleIndex = floor((k-1)/5)+1;
    lineStyle = styles{mod(styleIndex-1,numel(styles))+1};
    marker = markers{mod(styleIndex-1,numel(markers))+1};
    % realmin is a plotting floor only; measurements in the MAT stay exact.
    semilogy(residualAxes,curve.ImageScalars,max(curve.SurfaceMaxRelative,realmin), ...
        'Color',colors(colorIndex,:),'LineStyle',lineStyle,'Marker',marker, ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName',combinationNames(k));
    semilogy(strengthAxes,curve.ImageScalars,max(curve.MaxSourceMagnitude,realmin), ...
        'Color',colors(colorIndex,:),'LineStyle',lineStyle,'Marker',marker, ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName',combinationNames(k));
end
set([residualAxes,strengthAxes],'YScale','log','FontSize',11);
grid(residualAxes,'on'); grid(strengthAxes,'on');
xlabel(residualAxes,'Scalar unknowns on image nodes (proxy unknowns excluded)');
xlabel(strengthAxes,'Scalar unknowns on image nodes (proxy unknowns excluded)');
ylabel(residualAxes,'Maximum relative surface velocity residual');
ylabel(strengthAxes,'Maximum radius-scaled source magnitude');
title(residualAxes,'Independent surface validation');
title(strengthAxes,'Proxy and image source strengths');
legendHandle = legend(residualAxes,'Location','southoutside', ...
    'NumColumns',min(5,numel(combinationNames)),'Interpreter','none');
legendHandle.Layout.Tile = 'south';
title(layout,sprintf('N=%d, delta_{wall}/a=%g, R_p/a=%.2f, V=(%g,%g,%g), Omega=(%g,%g,%g)', ...
    parameters.N,parameters.delta_wall/radius,parameters.proxyRadius/radius,motion));
subtitle(layout,sprintf(['Proxy S always present (%d scalar unknowns); T = T_n + T_{n2} + T_t; ', ...
    'relative SVD cutoff %.0e'],3*parameters.N,parameters.svdCutoff));

%% GEOMETRY: COMPLETE SPHERE AND WALL-FACING CLOSE-UP
nodeFigure = figure('Name','Oct 27 sphere-wall image and collocation nodes', ...
    'Position',[100,100,1250,650]);
nodeLayout = tiledlayout(nodeFigure,1,2,'TileSpacing','compact','Padding','compact');
fullAxes = nexttile(nodeLayout); hold(fullAxes,'on');
[sx,sy,sz] = sphere(28);
surf(fullAxes,centre(1)+radius*sx,centre(2)+radius*sy,centre(3)+radius*sz, ...
    'FaceColor','none','EdgeColor',[0.7,0.7,0.7],'EdgeAlpha',0.2, ...
    'HandleVisibility','off');
extent = 1.15*radius;
patch(fullAxes,centre(1)+extent*[-1,1,1,-1],centre(2)+extent*[-1,-1,1,1], ...
    wallHeight*ones(1,4),[0.8,0.8,0.8],'FaceAlpha',0.2,'EdgeColor','none', ...
    'DisplayName','Wall');
scatter3(fullAxes,data.proxyPoints(:,1),data.proxyPoints(:,2),data.proxyPoints(:,3), ...
    5,[0.6,0.6,0.6],'filled','DisplayName','Proxy Stokeslets');
scatter3(fullAxes,data.collocationPoints(:,1),data.collocationPoints(:,2), ...
    data.collocationPoints(:,3),5,[0.1,0.4,0.8],'filled','DisplayName','Collocation');
scatter3(fullAxes,imagePoints(:,1),imagePoints(:,2),imagePoints(:,3), ...
    32,[0.9,0.25,0.05],'filled','DisplayName','Image-line nodes','Tag','ImageNodesFull');
axis(fullAxes,'equal'); grid(fullAxes,'on'); view(fullAxes,35,20);
xlabel(fullAxes,'x'); ylabel(fullAxes,'y'); zlabel(fullAxes,'z');
title(fullAxes,sprintf('Full sphere: %d proxy, %d collocation, %d image points', ...
    size(data.proxyPoints,1),size(data.collocationPoints,1),size(imagePoints,1)));
legend(fullAxes,'Location','southoutside');

closeAxes = nexttile(nodeLayout); hold(closeAxes,'on');
% A thin slab avoids superimposing nodes from the entire sphere in projection.
slabHalfWidth = 0.1*radius;
near = abs(data.collocationPoints(:,2)-centre(2)) <= slabHalfWidth;
proxyNear = abs(data.proxyPoints(:,2)-centre(2)) <= slabHalfWidth;
boundaryX = centre(1)+linspace(-0.3,0.3,200)*radius;
boundaryZ = centre(3)-sqrt(radius^2-(boundaryX-centre(1)).^2);
plot(closeAxes,boundaryX,boundaryZ,'k-','DisplayName','Sphere surface');
yline(closeAxes,wallHeight,'k--','DisplayName','Wall');
scatter(closeAxes,data.proxyPoints(proxyNear,1),data.proxyPoints(proxyNear,3), ...
    12,[0.6,0.6,0.6],'filled','DisplayName','Proxy Stokeslets');
scatter(closeAxes,data.collocationPoints(near,1),data.collocationPoints(near,3), ...
    14,[0.1,0.4,0.8],'filled','DisplayName','Collocation');
scatter(closeAxes,imagePoints(:,1),imagePoints(:,3),32,[0.9,0.25,0.05], ...
    'filled','DisplayName','Image-line nodes');
axis(closeAxes,'equal'); grid(closeAxes,'on');
xlim(closeAxes,centre(1)+[-0.25,0.25]*radius);
ylim(closeAxes,wallHeight+[-0.02,0.42]*radius);
xlabel(closeAxes,'x'); ylabel(closeAxes,'z');
title(closeAxes,'Wall-facing close-up');
subtitle(closeAxes,'x-z projection of nodes with |y-c_y| <= 0.1a');
legend(closeAxes,'Location','southoutside');
title(nodeLayout,sprintf('Image-line geometry at n_{images}=%d; fitting grid fixed for all solves', ...
    size(imagePoints,1)));

%% OPTIONAL EXPORT -- NUMERICAL RESULTS ARE NEVER CHANGED
if strlength(string(outputStem)) > 0
    figures = [resultFigure,nodeFigure];
    suffixes = {'_image_combinations','_nodes'};
    for k = 1:2
        stem = [char(outputStem),suffixes{k}];
        savefig(figures(k),[stem '.fig']);
        exportgraphics(figures(k),[stem '.png'],'Resolution',200);
        if k == 1
            exportgraphics(figures(k),[stem '.pdf'],'ContentType','vector');
        else
            % Rasterize the transparent 3-D scene to avoid slow PDF tessellation.
            exportgraphics(figures(k),[stem '.pdf'],'ContentType','image','Resolution',200);
        end
    end
end
end

function self_test()
% Check saved-file replotting, family-dependent abscissae and node selection.
parameters = struct('N',6,'radius',1,'mu',1,'wallHeight',0,'delta_wall',1e-3, ...
    'proxyRadius',0.63,'rigidMotion',[1;0;-1;0;0;0],'svdCutoff',1e-12);
names = ["S";"T_t";"T_nT_n2"];
imageCounts = [0;5;10];
results = table(repmat(names,3,1),repelem(imageCounts,3), ...
    'VariableNames',{'Combination','ImagePoints'});
results.ImageScalars = results.ImagePoints.*repmat([3;3;6],3,1);
results.SurfaceMaxRelative = repmat([1e-2;1e-3;1e-4],3,1);
results.MaxSourceMagnitude = (1:9).';
centre = [0,0,1.001];
proxyPoints = centre+0.63*[eye(3);-eye(3)];
collocationPoints = centre+[eye(3);-eye(3)];
linePoints = {zeros(0,3);[0,0,0.2];[0,0,0.1;0,0,0.2]};
filename = [tempname,'.mat'];
fileCleanup = onCleanup(@() delete(filename));
save(filename,'parameters','results','centre','proxyPoints','collocationPoints','linePoints');
[f,g] = plot_oct27_sphere_wall_image_combinations(filename);
figureCleanup = onCleanup(@() close([f,g]));
axesHandles = findall(f,'Type','axes');
for ax = axesHandles.'
    axial = findobj(ax,'Type','line','DisplayName','T_t');
    transverse = findobj(ax,'Type','line','DisplayName','T_nT_n2');
    if ~isequal(axial.XData,[0,15,30]) || ~isequal(transverse.XData,[0,30,60]) ...
            || ~strcmp(ax.YScale,'log')
        error('oct27:PlotTest','Result curves must use image scalar counts, not point counts.');
    end
end
nodes = findobj(g,'Tag','ImageNodesFull');
if numel(nodes.ZData) ~= 2
    error('oct27:PlotTest','Geometry must use the largest saved image count.');
end
fprintf('plot_oct27_sphere_wall_image_combinations PASS: saved-file reload, scalar counts, node plot.\n');
end
