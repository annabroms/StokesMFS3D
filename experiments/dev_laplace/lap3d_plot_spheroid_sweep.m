function sweepFigure = lap3d_plot_spheroid_sweep( ...
        rows,controls,delta,basicSolution,sweepFigure)
%LAP3D_PLOT_SPHEROID_SWEEP Create/update the node-beta-width comparison plot.
%   fig = lap3d_plot_spheroid_sweep(rows,controls,delta,basicSolution)
%   fig = lap3d_plot_spheroid_sweep(rows,controls,delta,basicSolution,fig)
%   lap3d_plot_spheroid_sweep   % synthetic-data graphics self-test
%
% rows is the sweep's struct array.  Each width gets two logarithmic panels:
% coefficient max norm and independent residual max norm, versus retained
% enhancing nodes per particle.  Colours distinguish requested beta values.
% Existing line objects are updated, so the figure can show partial results.

if nargin == 0
    self_test;
    sweepFigure = [];
    return
end

if nargin < 5 || isempty(sweepFigure)
    sweepFigure = figure('Name','Chopped-spheroid node/beta/width sweep');
    layout = tiledlayout(sweepFigure,numel(controls.widthScales),2, ...
        'TileSpacing','compact','Padding','compact');
    coefficientLines = gobjects(numel(controls.widthScales),numel(controls.betas));
    residualLines = gobjects(size(coefficientLines));
    coefficientAxes = gobjects(numel(controls.widthScales),1);
    residualAxes = gobjects(size(coefficientAxes));
    lineColors = lines(numel(controls.betas));
    for w = 1:numel(controls.widthScales)
        axCoefficient = nexttile(layout);
        axResidual = nexttile(layout);
        coefficientAxes(w) = axCoefficient;
        residualAxes(w) = axResidual;
        % Set the scales explicitly before hold, including for empty curves.
        set([axCoefficient,axResidual],'YScale','log');
        hold(axCoefficient,'on')
        hold(axResidual,'on')
        for b = 1:numel(controls.betas)
            label = sprintf('beta %.2g',controls.betas(b));
            coefficientLines(w,b) = plot(axCoefficient,NaN,NaN,'o-', ...
                'Color',lineColors(b,:),'DisplayName',label);
            residualLines(w,b) = plot(axResidual,NaN,NaN,'o-', ...
                'Color',lineColors(b,:),'DisplayName',label);
        end
        yline(axCoefficient,basicSolution.coefficientNorm,'k--','basic', ...
            'HandleVisibility','off');
        yline(axResidual,basicSolution.testMaxResidual,'k--','basic', ...
            'HandleVisibility','off');
        yline(axCoefficient,controls.coefficientTarget,'k:','target', ...
            'HandleVisibility','off');
        yline(axResidual,controls.residualTarget,'k:','target', ...
            'HandleVisibility','off');
        xlabel(axCoefficient,'enhancing nodes per particle, after chopping')
        xlabel(axResidual,'enhancing nodes per particle, after chopping')
        ylabel(axCoefficient,'max |c_j|, all sources')
        ylabel(axResidual,'independent relative residual, inf norm')
        title(axCoefficient,sprintf('Width scale %.3g',controls.widthScales(w)))
        title(axResidual,sprintf('Width scale %.3g',controls.widthScales(w)))
        grid(axCoefficient,'on')
        grid(axResidual,'on')
    end
    betaLegend = legend(coefficientLines(1,:),'Orientation','horizontal');
    betaLegend.Layout.Tile = 'south';
    sweepFigure.UserData = struct('layout',layout, ...
        'coefficientLines',coefficientLines,'residualLines',residualLines, ...
        'coefficientAxes',coefficientAxes,'residualAxes',residualAxes);
elseif ~isgraphics(sweepFigure)
    return
end

graphics = sweepFigure.UserData;
for w = 1:numel(controls.widthScales)
    for b = 1:numel(controls.betas)
        selected = [rows.valid] & [rows.widthScale] == controls.widthScales(w) & ...
            [rows.beta] == controls.betas(b);
        set(graphics.coefficientLines(w,b), ...
            'XData',[rows(selected).capNodes1], ...
            'YData',[rows(selected).coefficientInf]);
        set(graphics.residualLines(w,b), ...
            'XData',[rows(selected).capNodes1], ...
            'YData',[rows(selected).testMaxResidual]);
    end
    % Leave room above the basic/target reference lines so labels stay visible.
    selectedWidth = [rows.valid] & [rows.widthScale] == controls.widthScales(w);
    coefficientValues = [basicSolution.coefficientNorm,controls.coefficientTarget, ...
        rows(selectedWidth).coefficientInf];
    residualValues = [basicSolution.testMaxResidual,controls.residualTarget, ...
        rows(selectedWidth).testMaxResidual];
    coefficientValues = coefficientValues(coefficientValues > 0);
    residualValues = residualValues(residualValues > 0);
    ylim(graphics.coefficientAxes(w), ...
        [min(coefficientValues)/1.5,max(coefficientValues)*1.5]);
    ylim(graphics.residualAxes(w), ...
        [min(residualValues)/1.5,max(residualValues)*1.5]);
end
completed = nnz([rows.valid] | strlength([rows.message]) > 0);
title(graphics.layout,sprintf( ...
    'Chopped spheroids: radius %.3g, gap %.3g; %d/%d cases completed', ...
    controls.radius,delta,completed,numel(rows)))
drawnow limitrate
end


function self_test()
%SELF_TEST checks incremental updates and logarithmic axes without a solve.

controls = struct('radius',0.63,'widthScales',[0.5,1], ...
    'betas',[0.1,0.3],'coefficientTarget',1,'residualTarget',1e-3);
basicSolution = struct('coefficientNorm',10,'testMaxResidual',0.1);
[counts,betas,widths] = ndgrid([20,40,80],controls.betas,controls.widthScales);
row = struct('valid',false,'widthScale',NaN,'beta',NaN,'capNodes1',NaN, ...
    'coefficientInf',NaN,'testMaxResidual',NaN,'message',"");
rows = repmat(row,numel(counts),1);
sweepFigure = lap3d_plot_spheroid_sweep(rows,controls,0.01,basicSolution);
for k = 1:numel(counts)
    rows(k).valid = true;
    rows(k).widthScale = widths(k);
    rows(k).beta = betas(k);
    rows(k).capNodes1 = counts(k);
    rows(k).coefficientInf = 2*widths(k)*(1+betas(k))/sqrt(counts(k));
    rows(k).testMaxResidual = 0.01*(1+betas(k))/counts(k);
end
sweepFigure = lap3d_plot_spheroid_sweep( ...
    rows,controls,0.01,basicSolution,sweepFigure);
axesHandles = findall(sweepFigure,'Type','axes');
assert(numel(axesHandles) == 4);
assert(all(strcmp(get(axesHandles,'YScale'),'log')));
assert(isscalar(findall(sweepFigure,'Type','legend')));
graphics = sweepFigure.UserData;
assert(numel(graphics.coefficientLines(1,1).XData) == 3);
assert(isequal(graphics.coefficientLines(1,1).XData,[20,40,80]));
assert(isequal(graphics.residualLines(1,1).YData,[rows(1:3).testMaxResidual]));
title(graphics.layout,'Graphics self-test: synthetic metrics, not experiment results')
fprintf('lap3d_plot_spheroid_sweep self-test passed: incremental updates, logarithmic axes, and plotted values.\n');
end
