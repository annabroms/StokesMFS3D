function fig = plot_wall_geometry(q,a,rin,opt,validation)
%PLOT_WALL_GEOMETRY Show proxy geometry and surface velocity residuals.
% opt supplies wall_height, delta_body, delta_wall, and visible.
% validation is the output of wall_surface_validation.
if nargin == 0, self_test(); fig=[]; return; end
if nargin < 5, validation = []; end
fig = figure('Name','Wall-bounded proxy geometry','Visible',opt.visible);
layout = tiledlayout(fig,1,2);
[X,Y,Z] = sphere(30);
low = min(q(:,1:2),[],1)-2*a; high = max(q(:,1:2),[],1)+2*a;
ax = nexttile(layout); hold(ax,'on');
surf(ax,[low(1),high(1);low(1),high(1)], ...
    [low(2),low(2);high(2),high(2)], ...
    opt.wall_height*ones(2),'FaceAlpha',0.25,'FaceColor',[0.3,0.4,0.5], ...
    'EdgeColor','none');
for p = 1:size(q,1)
    surf(ax,a*X+q(p,1),a*Y+q(p,2),a*Z+q(p,3), ...
        'FaceAlpha',0.45,'EdgeColor','none','FaceColor',[0.75,0.4,0.2]);
end
n = size(rin,1)/size(q,1);
scatter3(ax,rin(1:n,1),rin(1:n,2),rin(1:n,3),5,'k','filled');
wallPoint = [q(1,1:2),opt.wall_height];
line(ax,[wallPoint(1),wallPoint(1)],[wallPoint(2),wallPoint(2)], ...
    [opt.wall_height,q(1,3)-a],'Color','k','LineWidth',2);
text(ax,wallPoint(1),wallPoint(2),opt.wall_height+opt.delta_wall/2, ...
    sprintf('  delta_{wall}=%.3g',opt.delta_wall));
if size(q,1)>1
    distances = vecnorm(q(2:end,:)-q(1,:),2,2); [~,index] = min(distances);
    other = q(index+1,:); direction = (other-q(1,:))/norm(other-q(1,:));
    endpoints = [q(1,:)+a*direction;other-a*direction];
    plot3(ax,endpoints(:,1),endpoints(:,2),endpoints(:,3),'k-','LineWidth',2);
    middle = mean(endpoints,1);
    text(ax,middle(1),middle(2),middle(3),sprintf(' delta_{body}=%.3g',opt.delta_body));
end
quiver3(ax,low(1),low(2),opt.wall_height,0,0,a,0,'k','LineWidth',2);
text(ax,low(1),low(2),opt.wall_height+a,'wall normal');
axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
view(ax,35,25); title(ax,'Particle geometry and proxy sources');

ax = nexttile(layout); hold(ax,'on');
for p = 1:size(q,1)
    surf(ax,a*X+q(p,1),a*Y+q(p,2),a*Z+q(p,3), ...
        'FaceAlpha',0.12,'EdgeColor','none','FaceColor',[0.75,0.4,0.2]);
end
if isempty(validation)
    title(ax,'Surface residuals not supplied');
else
    assert(size(validation.points,2)==3 && ...
        numel(validation.pointErrors)==size(validation.points,1));
    errors = validation.pointErrors(:);
    assert(all(isfinite(errors)) && all(errors>=0));
    positiveErrors = errors(errors>0);
    if isempty(positiveErrors)
        colorErrors = repmat(realmin,size(errors));
    else
        colorErrors = max(errors,min(positiveErrors));
    end
    scatter3(ax,validation.points(:,1),validation.points(:,2), ...
        validation.points(:,3),24,log10(colorErrors),'filled');  
    cb = colorbar(ax);
    clim([min(log10(colorErrors)),max(log10(colorErrors))]);
    if validation.isRelative
        colorbarLabel = 'log10(Relative velocity error)';
    else
        colorbarLabel = 'log10(Absolute velocity error)';
    end
    title(ax,sprintf('Surface residuals (max %.2e)',max(errors)));
    ylabel(cb,colorbarLabel);
end
axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
view(ax,35,25);
end

function self_test()
% Check both panels with known surface points and residual values.
q = [0,0,1.2;3,0,1.2];
o = struct('wall_height',0,'delta_body',1,'delta_wall',0.2, ...
    'visible',get(groot,'DefaultFigureVisible'));
rin = q+[0.63,0,0];
[X,Y,Z] = sphere(8);
points = [X(:),Y(:),Z(:)]+q(1,:);
validation = struct('points',points,'pointErrors',1e-3*(1+X(:).^2), ...
    'isRelative',true);
f = plot_wall_geometry(q,1,rin,o,validation);
assert(isgraphics(f));
fprintf('plot_wall_geometry: PASS; geometry and surface residual panels show 2 spheres.\n');
end
