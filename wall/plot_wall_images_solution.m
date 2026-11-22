function figures = plot_wall_images_solution(problem,lambda,details)
%PLOT_WALL_IMAGES_SOLUTION Show grids, every source type and node residuals.
% lambda is dimensionless. Plots use physical coordinates and separate
% colorbars for proxy S and every selected wall/pair family. Source panels
% show vector magnitudes, not scalar coefficient infinity norms.
% Collocation residuals are evaluated at the actual merged rout, with
% global rigid-speed normalization (or absolute for zero rigid motion).

a = problem.options.radius; zw = problem.options.wall_height;
q = a*problem.q+[0,0,zw];
figures = gobjects(3,1);
figures(1) = figure('Name',[problem.mode,' wall/pair image discretization'], ...
    'Position',[80,80,1400,520]);
tiledlayout(1,3,'TileSpacing','compact');
for viewIndex = 1:3
    ax = nexttile; draw_scene(ax,q,a,zw);
    nodes = gobjects(4,1);
    nodes(1) = scatter3(ax,details.rout(:,1),details.rout(:,2),details.rout(:,3), ...
        4,[0.2,0.4,0.75],'filled','DisplayName','rout: collocation');
    for k = 1:problem.P
        basis = problem.body{k}.basis;
        for j = 1:numel(basis)
            % Draw shared family locations once; their strengths get separate panels.
            if j>1 && isequal(basis(j).points,basis(j-1).points), continue; end
            rin = a*basis(j).points+[0,0,zw];
            color = [0.1,0.65,0.3]; marker = '.'; count = 15;
            group = 2;
            if startsWith(basis(j).label,"wall")
                color = [0.85,0.2,0.1]; marker = 'o'; count = 28;
                group = 3;
            elseif startsWith(basis(j).label,"pair")
                color = [0.65,0.15,0.8]; marker = 's'; count = 28;
                group = 4;
            end
            nodes(group) = scatter3(ax,rin(:,1),rin(:,2),rin(:,3),count,color,marker, ...
                'DisplayName',char(basis(j).label));
        end
    end
    if viewIndex == 1
        title(ax,sprintf('N=%d proxy; M=[%s]\nWall gap %.3g', ...
            problem.N,num2str(problem.M.'),problem.minWallGap));
        labels = {'rout: collocation','proxy S', ...
            char("wall "+problem.options.wallImageFamilies+" line"), ...
            char("pair "+problem.options.pairImageFamilies+" line")};
        active = isgraphics(nodes);
        legend(ax,nodes(active),labels(active),'Location','southoutside', ...
            'FontSize',8,'NumColumns',2);
    elseif viewIndex == 2
        xlim(ax,q(1,1)+a*[-0.5,0.5]); ylim(ax,q(1,2)+a*[-0.5,0.5]);
        zlim(ax,zw+a*[0,0.5]); view(ax,0,0);
        title(ax,'Wall-facing nodes: close-up');
    else
        if ~isempty(problem.pairs)
            pair = problem.pairs(1,:); mid = mean(q(pair,:),1);
            xlim(ax,mid(1)+a*[-0.4,0.4]); ylim(ax,mid(2)+a*[-0.5,0.5]);
            zlim(ax,mid(3)+a*[-0.6,0.6]); view(ax,0,0);
        end
        title(ax,sprintf('Pair gap %.3g: %s lines', ...
            problem.minBodyGap,problem.options.pairImageFamilies));
    end
end

% Collect plot labels from actual blocks, retaining separate stresslet
% orientations and combining the same pair family across all neighbours.
labels = strings(0,1);
for k = 1:problem.P
    basis = problem.body{k}.basis;
    for j = 1:numel(basis)
        label = basis(j).label;
        if startsWith(label,"pair")
            parts = split(label); label = "pair "+parts(3);
        end
        labels(end+1) = label; %#ok<AGROW>
    end
end
labels = unique(labels,'stable');
columns = min(4,numel(labels)+1); rows = ceil((numel(labels)+1)/columns);
figures(2) = figure('Name',[problem.mode,' source vector magnitudes'], ...
    'Position',[80,80,1400,max(520,350*rows)]);
tiledlayout(rows,columns,'TileSpacing','compact');
for panel = 1:numel(labels)
    ax = nexttile; draw_scene(ax,q,a,zw);
    rin = zeros(0,3); values = zeros(0,1);
    for k = 1:problem.P
        basis = problem.body{k}.basis;
        offset = problem.sourceOffsets(k);
        for j = 1:numel(basis)
            count = size(basis(j).points,1); cols = offset+(1:3*count);
            label = basis(j).label;
            if startsWith(label,"pair")
                parts = split(label); label = "pair "+parts(3);
            end
            if label == labels(panel)
                rin = [rin;a*basis(j).points+[0,0,zw]]; %#ok<AGROW>
                values = [values;vecnorm(reshape(lambda(cols),3,[]).',2,2)]; %#ok<AGROW>
            end
            offset = offset+3*count;
        end
    end
    if isempty(values)
        title(ax,labels(panel)+": inactive"); continue
    end
    % A display floor changes colors only; raw coefficients remain untouched.
    floorValue = max(1e-16,1e-12*max(values));
    scatter3(ax,rin(:,1),rin(:,2),rin(:,3),22,max(values,floorValue),'filled');
    set(ax,'ColorScale','log'); colormap(ax,turbo(256)); colorbar(ax);
    % Background surfaces have their own CData; use only source magnitudes
    % to set color limits, so their constant colors do not dilute the scale.
    lower = max(floorValue,min(values)); upper = max(max(values),10*lower);
    clim(ax,[lower,upper]);
    title(ax,sprintf('%s: max |lambda| = %.3g',labels(panel),max(values)));
end
ax = nexttile; axis(ax,'off');
text(ax,0,1,{'Dimensionless vector magnitudes'; ...
    'Physical S: multiply by mu*a';'Physical T/R: multiply by mu*a^2'; ...
    'Physical D: multiply by mu*a^3';'Every source includes its';'complete wall image'}, ...
    'VerticalAlignment','top','Interpreter','none','FontSize',10);

figures(3) = figure('Name',[problem.mode,' collocation surface residual'], ...
    'Position',[80,80,1400,520]);
tiledlayout(1,3,'TileSpacing','compact');
floorValue = max(1e-16,1e-12*max(details.collocationResidual));
for viewIndex = 1:3
    ax = nexttile; draw_scene(ax,q,a,zw);
    scatter3(ax,details.rout(:,1),details.rout(:,2),details.rout(:,3),12, ...
        max(details.collocationResidual,floorValue),'filled');
    set(ax,'ColorScale','log'); colormap(ax,turbo(256)); colorbar(ax);
    lower = max(floorValue,min(details.collocationResidual));
    upper = max(max(details.collocationResidual),10*lower);
    clim(ax,[lower,upper]);
    if viewIndex == 1
        title(ax,sprintf('Collocation max %.3e; independent %.3e', ...
            details.collocationError,details.surfaceError));
    elseif viewIndex == 2
        xlim(ax,q(1,1)+a*[-0.5,0.5]); ylim(ax,q(1,2)+a*[-0.5,0.5]);
        zlim(ax,zw+a*[0,0.5]); view(ax,0,0); title(ax,'Wall gap residual');
    else
        if ~isempty(problem.pairs)
            pair = problem.pairs(1,:); mid = mean(q(pair,:),1);
            xlim(ax,mid(1)+a*[-0.4,0.4]); ylim(ax,mid(2)+a*[-0.5,0.5]);
            zlim(ax,mid(3)+a*[-0.6,0.6]); view(ax,0,0);
        end
        title(ax,'Pair gap residual');
    end
end
if details.isRelative
    sgtitle(sprintf('%s: pointwise |u-u_rigid| / max |u_rigid| (speed %.3g)', ...
        problem.mode,details.normalizationSpeed),'Interpreter','none');
else
    sgtitle([problem.mode,': absolute pointwise |u-u_rigid| (zero rigid speed)'], ...
        'Interpreter','none');
end
end

function draw_scene(ax,q,a,zw)
% Reuse the same physical wall and translucent bodies in all views.
hold(ax,'on');
[x,y,z] = sphere(20);
for k = 1:size(q,1)
    surf(ax,q(k,1)+a*x,q(k,2)+a*y,q(k,3)+a*z,'FaceAlpha',0.06, ...
        'FaceColor',[0.7,0.75,0.8],'EdgeColor','none','HandleVisibility','off');
end
limits = [min(q(:,1:2),[],1)-a;max(q(:,1:2),[],1)+a];
patch(ax,limits([1,2,2,1],1),limits([1,1,2,2],2),zw*ones(4,1), ...
    [0.45,0.45,0.45],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
axis(ax,'equal'); grid(ax,'on'); view(ax,30,22);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
end
