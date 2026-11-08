function [figure3D,figureMeridional] = bubble_plot_geometry(geometry,opt,visibility)
%BUBBLE_PLOT_GEOMETRY Show interface, cap edge, sources, and collocation.
if nargin == 0
    bubble_self_test();
    figure3D = [];
    figureMeridional = [];
    return
end
if nargin<3, visibility = 'on'; end
c = opt.centre;
a = opt.radius;
psi = linspace(0,2*pi,120);
figure3D = figure('Name','Bubble geometry: 3D','Color','w', ...
    'Visible',visibility);
ax = axes(figure3D); hold(ax,'on');
for region = 1:2
    if region==1
        theta = linspace(0,opt.phi,35).';
        color = [0.75 0.45 0.7];
    else
        theta = linspace(opt.phi,pi,35).';
        color = [0.3 0.75 0.75];
    end
    x = c(1)+a*sin(theta)*cos(psi);
    y = c(2)+a*sin(theta)*sin(psi);
    z = c(3)+a*cos(theta)*ones(size(psi));
    surf(ax,x,y,z,'FaceColor',color,'FaceAlpha',0.12, ...
        'EdgeColor','none');
end
edge = c+a*[sin(opt.phi)*cos(psi(:)), ...
    sin(opt.phi)*sin(psi(:)),cos(opt.phi)*ones(numel(psi),1)];
plot3(ax,edge(:,1),edge(:,2),edge(:,3),'k-','LineWidth',1.5);
base = geometry.sourcePoints(geometry.sourceRingIndex==0,:);
cone = geometry.sourcePoints(geometry.sourceRingIndex>0,:);
scatter3(ax,base(:,1),base(:,2),base(:,3),16, ...
    [0.1 0.35 0.9],'filled','DisplayName','uniform proxy');
if ~isempty(cone)
    scatter3(ax,cone(:,1),cone(:,2),cone(:,3),22, ...
        [0.95 0.35 0.05],'filled','DisplayName','cone sources');
end
for region = 1:2
    if region==1
        x = geometry.points(geometry.cap,:);
        color = [0.7 0.1 0.5];
        name = 'cap collocation';
    else
        x = geometry.points(geometry.mobile,:);
        color = [0.05 0.55 0.3];
        name = 'mobile collocation';
    end
    scatter3(ax,x(:,1),x(:,2),x(:,3),10,color,'filled', ...
        'DisplayName',name);
end
axis(ax,'equal'); grid(ax,'on'); view(ax,35,25);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
title(ax,sprintf(['%d base + %d cone sites; %d coefficients; ', ...
    '%d collocation points'], ...
    geometry.baseSourceCount,geometry.coneSourceCount, ...
    geometry.coefficientCount,geometry.collocationCount));
legend(ax,'Location','eastoutside');

figureMeridional = figure('Name','Bubble geometry: meridional', ...
    'Color','w','Visible',visibility);
ax = axes(figureMeridional); hold(ax,'on');
theta = linspace(0,pi,400);
plot(ax,a*sin(theta),c(3)+a*cos(theta),'k-','LineWidth',1.2);
plot(ax,a*sin(opt.phi),c(3)+a*cos(opt.phi),'kp', ...
    'MarkerFaceColor','k','MarkerSize',12);
sets = {base,cone,geometry.points(geometry.cap,:), ...
    geometry.points(geometry.mobile,:)};
colors = {[0.1 0.35 0.9],[0.95 0.35 0.05], ...
    [0.7 0.1 0.5],[0.05 0.55 0.3]};
names = {'uniform proxy','cone sources','cap collocation', ...
    'mobile collocation'};
for k = 1:4
    x = sets{k};
    if isempty(x), continue; end
    rho = vecnorm(x(:,1:2)-c(1:2),2,2);
    scatter(ax,rho,x(:,3),17,colors{k},'filled', ...
        'DisplayName',names{k});
end
axis(ax,'equal'); grid(ax,'on');
xlabel(ax,'meridional radius'); ylabel(ax,'z');
title(ax,sprintf('Edge distance: nearest source %.3g a; collocation %.3g a', ...
    geometry.minimumSourceDistance/a,geometry.minimumEdgeDistance/a));
legend(ax,'Location','best');
end

function bubble_self_test()
% Check that both views plot a five-level interior source geometry.
o = bubble_options(struct('sourceRequest',50,'phi',pi/4, ...
    'n_level',5,'n_az_source',8,'n_az_colloc',12));
g = bubble_geometry(o);
[f3,f2] = bubble_plot_geometry(g,o,'off');
assert(isgraphics(f3) && isgraphics(f2));
assert(numel(findall(f3,'Type','Scatter'))==4);
close(f3); close(f2);
fprintf('bubble_plot_geometry self-test passed.\n');
end
