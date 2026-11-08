function ring_plot_sources(ax,Yrep,a,centre,radius,proxyRadius,arrowScale,titleText)
%RING_PLOT_SOURCES Shared per-node ring vectors a_r drawn with quiver3.
%
%   ring_plot_sources(ax,Yrep,a,centre,radius,proxyRadius,arrowScale,titleText)
%
% Inputs:  ax          axes handle;
%          Yrep        R-by-3 one representative source per ring (psi = 0, so
%                      in the meridional plane y = centre(2));
%          a           R-by-3 physical per-node Cartesian vectors (row r = a_r',
%                      force units); not completion or nullspace coordinates;
%          centre 1x3, radius (physical sphere), proxyRadius;
%          arrowScale  length per force unit (same value for all panels);
%          titleText   panel title.
% Draws the physical (solid) and proxy (dashed) sphere sections in the plane
% y = centre(2), in an oblique view so the out-of-plane component a_y shows.
if nargin == 0
    self_test();
    return
end
t = linspace(0,2*pi,400);
hold(ax,'on');
plot3(ax,centre(1)+radius*cos(t),centre(2)+0*t,centre(3)+radius*sin(t),'k-');
plot3(ax,centre(1)+proxyRadius*cos(t),centre(2)+0*t, ...
    centre(3)+proxyRadius*sin(t),'k--');
A = arrowScale*a;               % scale 0: draw the vectors exactly as given
quiver3(ax,Yrep(:,1),Yrep(:,2),Yrep(:,3),A(:,1),A(:,2),A(:,3),0, ...
    'LineWidth',1.2,'Color',[0.85 0.2 0.1]);
plot3(ax,Yrep(:,1),Yrep(:,2),Yrep(:,3),'k.','MarkerSize',10);
lim = 1.1*radius;
xlim(ax,centre(1)+[-lim lim]); ylim(ax,centre(2)+[-lim lim]); zlim(ax,centre(3)+[-lim lim]);
daspect(ax,[1 1 1]); grid(ax,'on'); view(ax,35,22);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
title(ax,{char(titleText),sprintf('arrow scale %.3g (length per unit force), max |a_r| = %.3g',arrowScale,max(vecnorm(a,2,2)))},'FontSize',9);
end

function self_test()
%SELF_TEST Draw synthetic arrows in a hidden figure.
fig = figure('Visible','off'); ax = axes(fig);
th = linspace(0.3,2.8,4).';
Y = [0.6*sin(th),zeros(4,1),0.6*cos(th)];
ring_plot_sources(ax,Y,[Y(:,1) Y(:,3) Y(:,2)],[0 0 0],1,0.6,0.2,'test');
assert(numel(findobj(ax,'Type','Quiver'))==1 && all(isfinite(xlim(ax))));
close(fig);
fprintf('ring_plot_sources self-test passed.\n');
end
