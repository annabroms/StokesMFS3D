function fig = drop_image_plot_grid(dropGrid,visibility)
%DROP_IMAGE_PLOT_GRID Show one drop's nodes and any stresslet orientations.
if nargin == 0
    self_test(); fig = []; return
end
if nargin < 2, visibility = 'on'; end
i = 1;
if ~isempty(dropGrid.pairs)
    i = dropGrid.pairs(1,1);
end
fig = figure('Name','Drop source and collocation grids', ...
    'Color','w','Visible',visibility);
ax = axes(fig);
body = dropGrid.body(i);
h = gobjects(5+numel(body.basis)+numel(body.interiorBasis),1);
nh = 5;
h(1) = scatter3(ax,body.collocation(:,1),body.collocation(:,2), ...
    body.collocation(:,3),9,[0.25 0.25 0.25],'filled');
hold(ax,'on');
x = body.basis(1).points;
h(2) = scatter3(ax,x(:,1),x(:,2),x(:,3),18,[0 0.35 0.9]);
x = body.imagePoints;
h(3) = scatter3(ax,x(:,1),x(:,2),x(:,3),28,[0.9 0.2 0.1],'filled');
x = body.interiorBasis(1).points;
h(4) = scatter3(ax,x(:,1),x(:,2),x(:,3),18,[0 0.6 0.45]);
x = body.exteriorImagePoints;
h(5) = scatter3(ax,x(:,1),x(:,2),x(:,3),28,[0.7 0.15 0.7],'filled');
% Each stresslet block has its own coefficients and source directions.
names = {'collocation','inner proxy','inner image','outer proxy','outer image'};
arrowLength = 0.13*dropGrid.radii(i);
bases = {body.basis,body.interiorBasis};
sides = ["inner","outer"];
colors = [0.65 0 0;0 0.35 0.85;0.2 0.55 0.1];
for side = 1:2
    basis = bases{side};
    for b = 2:numel(basis)
        family = basis(b).family;
        colorIndex = find(["T","T_n","T_t"]==family,1);
        if isempty(colorIndex), continue; end
        x = basis(b).points;
        d = arrowLength*basis(b).normalDirections;
        nh = nh+1;
        h(nh) = quiver3(ax,x(:,1),x(:,2),x(:,3), ...
            d(:,1),d(:,2),d(:,3),0,'Color',colors(colorIndex,:), ...
            'LineWidth',1.2);
        names{end+1} = char(sides(side)+" "+family+" direction"); %#ok<AGROW>
    end
end
axis(ax,'equal');
set(ax,'XGrid','on','YGrid','on','ZGrid','on');
view(ax,3);
title(ax,sprintf('Drop %d',i));
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
legend(ax,h(1:nh),names,'Location','best');
drawnow;
end

function self_test()
% Confirm both stresslet lines have arrows and the proxy-only case has none.
q = [0 0 0;2.1 0 0]; a = ones(2,1);
opt = drop_image_options(struct('N',30,'N_image',1, ...
    'capPointCounts',[12 12],'mirror_ext',true, ...
    'imageFamilies',"SRTD",'use_orth_dir_images',true),2,a);
g = drop_image_discretization(q,a,opt);
f = drop_image_plot_grid(g,'off');
assert(numel(findall(f,'Type','Scatter'))==5);
assert(isscalar(findall(f,'Type','Axes')));
assert(numel(findall(f,'Type','Quiver'))==2);
close(f);
opt.imageFamilies = "SDRT_nT_t";
g = drop_image_discretization(q,a,opt);
f = drop_image_plot_grid(g,'off');
if numel(findall(f,'Type','Quiver'))~=4
    error('drop_image_plot_grid:SelfTest','Both directions must appear on both lines.');
end
close(f);
opt.useImages = false;
g = drop_image_discretization(q,a,opt);
f = drop_image_plot_grid(g,'off');
assert(numel(findall(f,'Type','Scatter'))==5);
assert(isempty(findall(f,'Type','Quiver')));
close(f);
fprintf('drop_image_plot_grid self-test passed.\n');
end
