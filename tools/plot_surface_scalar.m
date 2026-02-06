function plot_surface_scalar(rvec_out, M, P, values, plot_title, cmap_name)
%PLOT_SURFACE_SCALAR Plot a scalar field on per-body surface nodes.
%
%   plot_surface_scalar(rvec_out, M, P, values, plot_title, cmap_name)
%
%   rvec_out   - (M*P) x 3 array of surface nodes, stacked by body
%   M          - number of nodes per body
%   P          - number of bodies
%   values     - (M*P) x 1 scalar values at surface nodes
%   plot_title - title string
%   cmap_name  - colormap name (e.g. 'parula', 'hot')

figure();
hold on
for k = 1:P
    idx = (k-1)*M+1:k*M;
    verts = rvec_out(idx,:);
    faces = convhull(verts(:,1), verts(:,2), verts(:,3));
    trisurf(faces, verts(:,1), verts(:,2), verts(:,3), values(idx), ...
        'EdgeColor','none', 'FaceColor','interp');
end
axis equal
grid on
view(-30,30);
colorbar
colormap(cmap_name)
title(plot_title)
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');
end
