function [figure3D,figureMeridional] = bubble_plot_strengths( ...
    solution,opt,visibility)
%BUBBLE_PLOT_STRENGTHS Separate coincident cone families in two views.
if nargin == 0
    bubble_self_test();
    figure3D = [];
    figureMeridional = [];
    return
end
if nargin<3, visibility = 'on'; end
basis = solution.basis;
assert(basis(1).family=="S" && norm(solution.Fg)>0);
coefficients = image_source_coefficients(basis,solution.lambda);
families = [basis.family];
stokesletBlocks = find(families=="S");
otherBlocks = find(families~="S");
panelBlocks = [{stokesletBlocks},num2cell(otherBlocks)];
panelCount = numel(panelBlocks);
values = cell(numel(basis),1);
for b = 1:numel(basis)
    family = basis(b).family;
    if ismember(family,["R","T","T_t","T_n","T_n2"])
        power = 1;
    elseif ismember(family,["D","Q_s"])
        power = 2;
    else
        power = 0;
    end
    reference = norm(solution.Fg)*opt.radius^power;
    values{b} = log10(max( ...
        vecnorm(coefficients{b},2,2)/reference,1e-16));
end
figure3D = figure('Name','Bubble source strengths: 3D', ...
    'Color','w','Visible',visibility, ...
    'Position',[40 80 max(700,340*panelCount) 560]);
figureMeridional = figure( ...
    'Name','Bubble source strengths: meridional', ...
    'Color','w','Visible',visibility, ...
    'Position',[60 100 max(700,340*panelCount) 560]);
figures = [figure3D,figureMeridional];
for viewIndex = 1:2
    layout = tiledlayout(figures(viewIndex),1,panelCount, ...
        'TileSpacing','compact','Padding','compact');
    for panel = 1:panelCount
        ax = nexttile(layout,panel);
        hold(ax,'on');
        blocks = panelBlocks{panel};
        combined = vertcat(values{blocks});
        limits = [min(combined),max(combined)];
        if limits(1)==limits(2), limits = limits+[-0.5 0.5]; end
        if viewIndex==1
            [sx,sy,sz] = sphere(30);
            surf(ax,opt.centre(1)+opt.radius*sx, ...
                opt.centre(2)+opt.radius*sy, ...
                opt.centre(3)+opt.radius*sz, ...
                'FaceAlpha',0.06,'EdgeColor','none', ...
                'FaceColor',[0.4 0.7 0.8]);
        else
            theta = linspace(0,pi,300);
            plot(ax,opt.radius*sin(theta), ...
                opt.centre(3)+opt.radius*cos(theta),'k-');
        end
        for b = blocks
            points = basis(b).points;
            marker = 'o';
            if b>1, marker = '^'; end
            if viewIndex==1
                scatter3(ax,points(:,1),points(:,2),points(:,3), ...
                    25,values{b},marker,'filled');
            else
                rho = vecnorm(points(:,1:2)-opt.centre(1:2),2,2);
                scatter(ax,rho,points(:,3),25,values{b}, ...
                    marker,'filled');
            end
        end
        clim(ax,limits); colormap(ax,turbo);
        cb = colorbar(ax);
        if panel==1
            label = '||F_g||';
            title(ax,'Base and cone Stokeslets');
        else
            family = basis(blocks).family;
            if ismember(family,["R","T","T_t","T_n","T_n2"])
                label = '||F_g|| a';
            elseif ismember(family,["D","Q_s"])
                label = '||F_g|| a^2';
            else
                label = '||F_g|| (distance-scaled derivative)';
            end
            title(ax,"Cone "+family,'Interpreter','none');
        end
        cb.Label.String = ['log_{10} max(||coefficient|| / (' ...
            label '), 10^{-16})'];
        axis(ax,'equal'); grid(ax,'on');
        if viewIndex==1
            view(ax,35,25);
            xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
        else
            xlabel(ax,'meridional radius'); ylabel(ax,'z');
        end
    end
end
end

function bubble_self_test()
% Confirm one Stokeslet scale and separate panels for coincident families.
o = bubble_options(struct('sourceRequest',50,'phi',pi/4, ...
    'n_level',1,'n_az_source',6,'n_az_colloc',8));
g = bubble_geometry(o);
s.grid = g;
s.basis = bubble_basis(g,o);
s.lambda = (1:g.coefficientCount).'/g.coefficientCount;
s.Fg = [0;0;-1];
[f3,f2] = bubble_plot_strengths(s,o,'off');
assert(numel(findall(f3,'Type','Scatter'))==2);
assert(numel(findall(f2,'Type','Scatter'))==2);
close(f3); close(f2);
o.refinementFamilies = ["S","T_t","T_n","T_n2","Q_s"];
o.radius = 2;
o = bubble_options(o);
g = bubble_geometry(o);
s.grid = g;
s.basis = bubble_basis(g,o);
s.lambda = ones(g.coefficientCount,1);
[f3,f2] = bubble_plot_strengths(s,o,'off');
assert(numel(findall(f3,'Type','Scatter'))==6);
assert(numel(findall(f2,'Type','Scatter'))==6);
axesHandles = findall(f3,'Type','Axes');
quadrupoleAxes = axesHandles(arrayfun( ...
    @(ax) string(ax.Title.String)=="Cone Q_s",axesHandles));
assert(isscalar(quadrupoleAxes));
quadrupolePoints = findobj(quadrupoleAxes,'Type','Scatter');
physical = image_source_coefficients(s.basis,s.lambda);
expectedColor = log10(norm(physical{end}(1,:))/ ...
    (norm(s.Fg)*o.radius^2));
assert(abs(quadrupolePoints.CData(1)-expectedColor)<1e-12);
close(f3); close(f2);
fprintf('bubble_plot_strengths self-test passed.\n');
end
