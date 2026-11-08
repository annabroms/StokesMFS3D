function figureHandle = plot_stagnant_cap_source_strengths( ...
    solution,label,visibility,referenceForce,colorLimits)
%PLOT_STAGNANT_CAP_SOURCE_STRENGTHS Show each source family without overlap.
% Rows separate the two flow representations. Columns separate base proxies
% and added families, whose source locations may coincide exactly.
if nargin==0
    self_test();
    figureHandle = [];
    return
end
if nargin<2, label = ''; end
if nargin<3, visibility = 'on'; end
if nargin<4, referenceForce = []; end
if nargin<5, colorLimits = []; end
if ~isempty(referenceForce)
    assert(isscalar(referenceForce) && referenceForce>0 && ...
        isfinite(referenceForce));
end

bases = {solution.exteriorBasis,solution.interiorBasis};
strengths = {solution.lambda,solution.gamma};
families = strings(1,0);
for side = 1:2
    for b = 2:numel(bases{side})
        family = bases{side}(b).family;
        if ~any(families==family)
            families(end+1) = family; %#ok<AGROW>
        end
    end
end
nColumns = 1+numel(families);
if ~isempty(colorLimits)
    assert(iscell(colorLimits) && isequal(size(colorLimits),[2 nColumns]));
end
figureHandle = figure('Name','Stagnant-cap source strengths', ...
    'Color','w','Visible',visibility);
layout = tiledlayout(figureHandle,2,nColumns, ...
    'TileSpacing','compact','Padding','compact');
sideNames = {'Exterior flow: sources inside drop', ...
    'Interior flow: sources outside drop'};
for side = 1:2
    basis = bases{side};
    coefficients = image_source_coefficients(basis,strengths{side});
    for familyIndex = 0:numel(families)
        ax = nexttile(layout,(side-1)*nColumns+familyIndex+1);
        hold(ax,'on');
        panelMinimum = inf;
        panelMaximum = -inf;
        for b = 1:numel(basis)
            if familyIndex==0
                selected = b==1;
            else
                selected = b>1 && basis(b).family==families(familyIndex);
            end
            if ~selected || isempty(basis(b).points), continue; end
            points = basis(b).points;
            magnitude = vecnorm(coefficients{b},2,2);
            if isempty(referenceForce)
                value = log10(max(magnitude,realmin));
            else
                value = log10(max(magnitude/referenceForce,1e-16));
            end
            panelMinimum = min(panelMinimum,min(value));
            panelMaximum = max(panelMaximum,max(value));
            scatter3(ax,points(:,1),points(:,2),points(:,3), ...
                18,value,'filled');
        end
        if familyIndex==0
            familyName = 'base S';
        else
            familyName = sprintf('added %s',families(familyIndex));
        end
        title(ax,sprintf('%s: %s',sideNames{side},familyName));
        axis(ax,'equal'); view(ax,3); grid(ax,'on');
        xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
        if isfinite(panelMinimum)
            if ~isempty(colorLimits)
                clim(ax,colorLimits{side,familyIndex+1});
            elseif panelMinimum==panelMaximum
                clim(ax,panelMinimum+[-0.5 0.5]);
            else
                clim(ax,[panelMinimum panelMaximum]);
            end
            cb = colorbar(ax);
            if isempty(referenceForce)
                cb.Label.String = 'log_{10} ||scaled coefficient||_2';
            else
                cb.Label.String = 'log_{10} max(||b||_2/||F_g||, 10^{-16})';
            end
        end
    end
end
colormap(figureHandle,turbo);
sgtitle(layout,sprintf('%s source coefficients',label), ...
    'Interpreter','none');
end

function self_test()
%SELF_TEST Check mixed-family block ordering, scaling, and separate panels.
block = struct('family',"S",'points',[0 0 0;0 0 1], ...
    'normalDirections',zeros(0,3),'scale',1);
basis = repmat(block,4,1);
basis(2).points = [1 0 0];
basis(3).family = "R";
basis(3).points = basis(2).points;
basis(3).scale = 2;
basis(4).family = "D";
basis(4).points = basis(2).points;
basis(4).scale = 0.5;
solution.exteriorBasis = basis;
solution.interiorBasis = basis;
% Each source has interleaved x,y,z coefficients, in basis-block order.
solution.lambda = [3;4;0;0;0;0;0;0;7;0;3;0;10;0;0];
solution.gamma = solution.lambda;
figureHandle = plot_stagnant_cap_source_strengths( ...
    solution,'mixed self-test','off');
axesHandles = findall(figureHandle,'Type','Axes');
assert(numel(axesHandles)==8);
expected = containers.Map({'base S','added S','added R','added D'}, ...
    {log10(5),log10(7),log10(6),log10(5)});
for side = 1:2
    for name = string(keys(expected))
        if side==1
            titleText = "Exterior flow: sources inside drop: "+name;
        else
            titleText = "Interior flow: sources outside drop: "+name;
        end
        matching = axesHandles(arrayfun( ...
            @(ax) string(ax.Title.String)==titleText,axesHandles));
        assert(isscalar(matching));
        plotted = findobj(matching,'Type','Scatter');
        assert(isscalar(plotted));
        count = 1;
        if name=="base S", count = 2; end
        assert(numel(plotted.XData)==count);
        if name~="base S"
            assert(plotted.XData==1 && plotted.YData==0 && ...
                plotted.ZData==0);
        end
        assert(abs(max(plotted.CData)-expected(char(name)))<1e-12);
        limits = clim(matching);
        assert(all(plotted.CData>=limits(1) & plotted.CData<=limits(2)));
        if name~="base S"
            assert(max(abs(limits-(expected(char(name))+[-0.5 0.5]))) ...
                <1e-12);
        end
    end
end
close(figureHandle);
limits = repmat({[-16 2]},2,4);
figureHandle = plot_stagnant_cap_source_strengths( ...
    solution,'normalized self-test','off',5,limits);
plotted = findall(figureHandle,'Type','Scatter');
assert(all(cellfun(@(v) isequal(v,[-16 2]), ...
    arrayfun(@(p) clim(ancestor(p,'axes')),plotted, ...
    'UniformOutput',false))));
close(figureHandle);
solution.exteriorBasis = basis(1);
solution.interiorBasis = basis(1);
solution.lambda = solution.lambda(1:6);
solution.gamma = solution.gamma(1:6);
figureHandle = plot_stagnant_cap_source_strengths( ...
    solution,'uniform self-test','off');
assert(numel(findall(figureHandle,'Type','Axes'))==2);
assert(numel(findall(figureHandle,'Type','Scatter'))==2);
close(figureHandle);
fprintf('plot_stagnant_cap_source_strengths self-test passed.\n');
end
