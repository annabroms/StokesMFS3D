function ring_band_report(results,x,labels,xLabel,mapIndices,Uref,mu,a)
%RING_BAND_REPORT Report and plot comparable mixed-band case diagnostics.
% results is a cell vector from ring_family_case; x is its numeric abscissa,
% labels are case labels, xLabel describes x, and mapIndices select audit maps.
% Uref is speed, mu viscosity and a radius. Velocities are divided by Uref,
% traction by mu*Uref/a, force by |Fg|, and torque by a*|Fg|. Coefficient
% magnitudes use each family's physical scale; maps use Cartesian coordinates.
if nargin==0
    self_test();
    return
end
n = numel(results);
if any(mapIndices<1 | mapIndices>n | mapIndices~=round(mapIndices))
    error('ring_band_report:MapCases','Map indices must refer to existing cases.');
end
metrics = zeros(n,4);
moments = zeros(n,3);
families = strings(1,0);
for j = 1:n
    r = results{j};
    metrics(j,:) = [r.Uerror,r.auditResidual];
    moments(j,:) = [r.forceResidual,r.torqueResidual,r.swirl];
    families = union(families,r.familyNames,'stable');
    fprintf('\n%s: U=[%.6g %.6g %.6g], U relative error %.3e\n',labels(j),r.U,r.Uerror);
    fprintf('  check cap/normal/traction: %.3e %.3e %.3e\n',r.residual);
    fprintf('  audit cap/normal/traction: %.3e %.3e %.3e; rank %d/%d\n', ...
        r.auditResidual,r.rank,r.unknowns);
    fprintf('  force %.3e, torque %.3e, swirl %.3e\n',moments(j,:));
    fprintf('  regional audit maxima [near edge, far edge, polar], rows cap/normal/traction:\n');
    disp(r.regional);
    if any(r.auditResidual>1.5*max(r.residual,eps))
        fprintf('  WARNING: audit increases a boundary maximum by more than 50%%.\n');
    end
    for k = 1:numel(r.familyNames)
        fprintf('  %s max physical coefficient / scale = %.3e (scale %.3e)\n', ...
            r.familyNames(k),r.coefficientMax(k),r.coefficientScale(k));
    end
end
fig = figure('Color','w','Position',[40 80 1400 500]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
semilogy(ax,x,max(metrics,eps),'o-'); grid(ax,'on');
legend(ax,{'U error','cap velocity','mobile normal','mobile traction'},'Location','best');
ylabel(ax,'relative error / normalized maximum'); title(ax,'Independent audit');
format_axis(ax,x,labels,xLabel);
ax = nexttile(tl);
semilogy(ax,x,max(moments,eps),'o-'); grid(ax,'on');
legend(ax,{'force / |F_g|','torque / (a|F_g|)','swirl / Uref'},'Location','best');
title(ax,'Moment and swirl diagnostics'); format_axis(ax,x,labels,xLabel);
ax = nexttile(tl); hold(ax,'on');
for k = 1:numel(families)
    values = nan(n,1);
    for j = 1:n
        index = find(results{j}.familyNames==families(k),1);
        if ~isempty(index), values(j) = max(results{j}.coefficientMax(index),eps); end
    end
    semilogy(ax,x,values,'o-','DisplayName',families(k));
end
grid(ax,'on'); ylabel(ax,'max per-node coefficient / family scale');
title(ax,'Physical coefficient magnitudes'); legend(ax,'Location','best');
format_axis(ax,x,labels,xLabel);
for j = unique(mapIndices(:).')
    r = results{j};
    ring_plot_residuals(r.mapPoints,r.mapResiduals, ...
        {char(labels(j)+" cap |u-U|/Uref"),char(labels(j)+" mobile |(u-U).n|/Uref"), ...
        char(labels(j)+" mobile tangential traction/(mu Uref/a)")});
end
fprintf('Uref=%.6g; traction scale=%.6g. Maps use independent audit nodes.\n',Uref,mu*Uref/a);
end

function format_axis(ax,x,labels,xLabel)
%FORMAT_AXIS Label family cases without obscuring numerical depth axes.
xlabel(ax,xLabel);
if xLabel=="Family combination"
    xticks(ax,x); xticklabels(ax,labels); xtickangle(ax,45);
end
end

function self_test()
%SELF_TEST Exercise reporting and figures without solving a physical problem.
r = struct('U',[0;0;-1],'Uerror',0.1,'auditResidual',[0.1 0.2 0.3], ...
    'residual',[0.1 0.2 0.3],'forceResidual',0,'torqueResidual',0, ...
    'swirl',0,'rank',3,'unknowns',3,'regional',ones(3,3), ...
    'familyNames',["GL S","S"],'coefficientMax',[1 2], ...
    'coefficientScale',[0.1 0.1],'mapPoints',{{[1 0 0],[0 1 0],[0 0 1]}}, ...
    'mapResiduals',{{0.1,0.2,0.3}});
old = get(groot,'defaultFigureVisible');
cleanup = onCleanup(@() set(groot,'defaultFigureVisible',old));
set(groot,'defaultFigureVisible','off');
existing = findall(groot,'Type','figure');
ring_band_report({r,r},[0.1 0.2],["a","b"],"Depth / a",1,1,1,1);
close(setdiff(findall(groot,'Type','figure'),existing));
fprintf('ring_band_report self-test passed.\n');
end
