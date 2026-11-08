function figureHandle = plot_stagnant_cap_collocation(geometry,opt,visibility)
%PLOT_STAGNANT_CAP_COLLOCATION Show every collocation node on the full sphere.
if nargin == 0
    self_test();
    return
end
if nargin < 3, visibility = 'on'; end
opt = stagnant_cap_options(opt);
assert(isequal(size(geometry.points),[geometry.collocationCount 3]));
assert(nnz(geometry.baseCollocation)+nnz(geometry.ringCollocation)== ...
    geometry.collocationCount);

a = opt.radius;
c = opt.centre;
[sx,sy,sz] = sphere(48);
capSurface = acos(max(-1,min(1,sz)))<opt.phi;
xc = c(1)+a*sx; yc = c(2)+a*sy; zc = c(3)+a*sz;
xm = xc; ym = yc; zm = zc;
xc(~capSurface)=NaN; yc(~capSurface)=NaN; zc(~capSurface)=NaN;
xm(capSurface)=NaN; ym(capSurface)=NaN; zm(capSurface)=NaN;
psi = linspace(0,2*pi,360);
edge = c+a*[sin(opt.phi)*cos(psi(:)), ...
    sin(opt.phi)*sin(psi(:)),cos(opt.phi)*ones(numel(psi),1)];
baseCap = geometry.points(geometry.baseCollocation & geometry.cap,:);
baseMobile = geometry.points(geometry.baseCollocation & geometry.mobile,:);
ringCap = geometry.points(geometry.ringCollocation & geometry.cap,:);
ringMobile = geometry.points(geometry.ringCollocation & geometry.mobile,:);

figureHandle = figure('Name',sprintf('Collocation geometry N=%d',geometry.N), ...
    'Color','w','Visible',visibility,'Position',[80 80 1400 650]);
layout = tiledlayout(figureHandle,1,2,'TileSpacing','compact', ...
    'Padding','compact');
viewAngles = [35 20;215 20];
for k = 1:2
    ax = nexttile(layout,k);
    hold(ax,'on');
    surf(ax,xc,yc,zc,'FaceColor',[0.66 0.38 0.73], ...
        'FaceAlpha',0.10,'EdgeColor','none');
    surf(ax,xm,ym,zm,'FaceColor',[0.17 0.72 0.67], ...
        'FaceAlpha',0.08,'EdgeColor','none');
    if opt.phi>0 && opt.phi<pi
        plot3(ax,edge(:,1),edge(:,2),edge(:,3),'k-', ...
            'LineWidth',1.6);
    end
    baseCapHandle = scatter3(ax,baseCap(:,1),baseCap(:,2), ...
        baseCap(:,3),14,[0.55 0.04 0.65],'s');
    baseMobileHandle = scatter3(ax,baseMobile(:,1),baseMobile(:,2), ...
        baseMobile(:,3),12,[0.02 0.52 0.42],'o');
    ringCapHandle = scatter3(ax,ringCap(:,1),ringCap(:,2), ...
        ringCap(:,3),24,[0.8 0.07 0.48],'d','filled');
    ringMobileHandle = scatter3(ax,ringMobile(:,1),ringMobile(:,2), ...
        ringMobile(:,3),22,[0.0 0.68 0.32],'^','filled');
    axis(ax,'equal');
    xlim(ax,c(1)+a*[-1.1 1.1]);
    ylim(ax,c(2)+a*[-1.1 1.1]);
    zlim(ax,c(3)+a*[-1.1 1.1]);
    view(ax,viewAngles(k,:));
    set(ax,'XGrid','on','YGrid','on','ZGrid','on');
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    if k==1
        title(ax,'View from one side');
        legend(ax,[baseCapHandle,baseMobileHandle, ...
            ringCapHandle,ringMobileHandle], ...
            {'base cap','base mobile','added cap rings', ...
            'added mobile rings'},'Location','southoutside', ...
            'NumColumns',2);
    else
        title(ax,'Opposite view');
    end
end
sgtitle(layout,sprintf(['Full-sphere collocation: N=%d, M=%d ', ...
    '(base %d + rings %d), \\phi=%.3g rad'], ...
    geometry.N,geometry.collocationCount,geometry.baseCount, ...
    nnz(geometry.ringCollocation),opt.phi));
end

function self_test()
%SELF_TEST Verify both views include all base and added ring nodes.
o = stagnant_cap_options(struct('sourceRequest',50,'phi',pi/3, ...
    'ringLevels',3));
g = stagnant_cap_geometry(o);
f = plot_stagnant_cap_collocation(g,o,'off');
assert(isgraphics(f));
scatters = findall(f,'Type','Scatter');
assert(numel(scatters)==8);
assert(sum(arrayfun(@(h) numel(h.XData),scatters))== ...
    2*g.collocationCount);
assert(numel(findall(f,'Type','Surface'))==4);
close(f);
fprintf('plot_stagnant_cap_collocation self-test passed.\n');
end
