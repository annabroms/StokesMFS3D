function areas = collocation_patch_area(points,bodyIndices,bodies)
%COLLOCATION_PATCH_AREA Estimate physical area per collocation node.
%   areas = collocation_patch_area(points,bodyIndices,bodies)
%   collocation_patch_area()  % numerical self-test and spheroid figure
%
% POINTS is the complete, deduplicated N-by-3 collocation grid;
% BODYINDICES assigns each row to one entry of BODIES. Each body has
% semi-axes axes (3-by-1), rotation Q (3-by-3), and centre c (3-by-1).
% AREAS is N-by-1 in the same row order as POINTS. This routine estimates
% local surface area; it does not generate collocation points or a cap.
%
% Use the complete grid on each ellipsoid/sphere,
% including any merged contact patch. Map it to the unit parameter sphere,
% triangulate, and share each spherical triangle's solid angle equally among
% its three vertices. Multiply by the ellipsoid map's surface Jacobian.
% Thus denser patches receive smaller area per node without counting their
% overlap with the basic grid twice. These are positive lumped area estimates,
% not the high-order quadrature weights of the original unrefined grid.
% For a sphere of radius R, their sum is exactly 4*pi*R^2 (up to roundoff).
% No source positions or source quadrature weights are involved.

if nargin == 0
    self_test();
    areas = [];
    return
end

bodyIndices = bodyIndices(:);
assert(size(points,2) == 3 && numel(bodyIndices) == size(points,1), ...
    'Surface coordinates and body-index lengths must agree.');
assert(all(ismember(bodyIndices,1:numel(bodies))), ...
    'Every collocation node must belong to one of the supplied bodies.');
areas = zeros(size(points,1),1);
for body = 1:numel(bodies)
    rows = find(bodyIndices == body);
    E = bodies{body};
    axes = E.axes(:);
    directions = ((E.Q.'*(points(rows,:).'-E.c))./axes).';
    lengths = vecnorm(directions,2,2);
    assert(numel(rows) >= 4 && all(abs(lengths-1) < 1e-8), ...
        'Area estimation requires a full grid of nodes on each body surface.');
    directions = directions./lengths; % remove surface-mapping roundoff

    faces = convhull(directions);
    u = directions(faces(:,1),:);
    v = directions(faces(:,2),:);
    w = directions(faces(:,3),:);
    solidAngles = 2*atan2(abs(dot(u,cross(v,w,2),2)), ...
        1+dot(u,v,2)+dot(v,w,2)+dot(w,u,2));
    assert(abs(sum(solidAngles)-4*pi) < 1e-8, ...
        'Collocation nodes must cover the complete parameter sphere.');
    angularAreas = accumarray(faces(:),repmat(solidAngles/3,3,1), ...
        [numel(rows),1],@sum,0);
    surfaceJacobian = prod(axes)*vecnorm(directions./axes.',2,2);
    areas(rows) = angularAreas.*surfaceJacobian;
end
assert(all(isfinite(areas) & areas > 0), ...
    ['Every collocation node must have positive finite area. Check for ', ...
     'duplicate or numerically indistinguishable collocation nodes.']);
end

function self_test()
%SELF_TEST checks areas and shows the finer weighting in a spheroidal cap.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
oldPath = path;
pathCleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'geometry','ellipsoids'));

E = struct('axes',[1.5;0.75;0.75], ...
    'Q',eye(3),'c',zeros(3,1));
n = 700;
k = (0:n-1).';
z = 1-2*(k+0.5)/n;
phi = 2*pi*k/((1+sqrt(5))/2);
dirs = [sqrt(1-z.^2).*cos(phi), ...
    sqrt(1-z.^2).*sin(phi),z];
base = dirs.*E.axes.';

opt = struct('collocationPatchMode',"exponential", ...
    'collocationPatchAngle',pi/5,'collocationPatchPointCount',241, ...
    'collocationPatchLayers',12,'collocationExponent',3, ...
    'collocationMinAzimuth',6,'collocationAzimuthOversampling',1);
cap = collocation_surface_patch(E,[E.axes(1);0;0],opt);
points = [base;cap];
areas = collocation_patch_area(points,ones(size(points,1),1),{E});
assert(all(isfinite(areas) & areas > 0));
assert(mean(areas(n+1:end)) < mean(areas(1:n)));

% Exact prolate-spheroid area checks the Jacobian and solid-angle sum.
a = E.axes(1);
b = E.axes(2);
e = sqrt(1-(b/a)^2);
exactArea = 2*pi*b^2*(1+a*asin(e)/(b*e));
assert(abs(sum(areas)-exactArea)/exactArea < 0.02);

figure('Name','Collocation area on a prolate spheroid');
[sx,sy,sz] = sphere(36);
colour = log10(areas);
limits = [min(colour),max(colour)];
theta = linspace(0,2*pi,240);
beta = opt.collocationPatchAngle;
rim = [a*cos(beta)*ones(size(theta)); ...
    b*sin(beta)*cos(theta);b*sin(beta)*sin(theta)];
for panel = 1:2
    ax = subplot(1,2,panel);
    surf(ax,a*sx,b*sy,b*sz,'FaceAlpha',0.08, ...
        'FaceColor',[0.6,0.7,0.8],'EdgeColor','none');
    hold(ax,'on');
    if panel == 1
        rows = true(size(areas));
    else
        rows = points(:,1) > a*cos(1.2*beta);
    end
    scatter3(ax,points(rows,1),points(rows,2),points(rows,3), ...
        13,colour(rows),'filled');
    plot3(ax,rim(1,:),rim(2,:),rim(3,:),'k-', 'LineWidth',1.3);
    axis(ax,'equal');
    grid(ax,'on');
    view(ax,3);
    clim(ax,limits);
    colorbar(ax);
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    if panel == 1
        title(ax,'Complete spheroid: log_{10} node area');
    else
        xlim(ax,[a-0.5,a+0.04]);
        ylim(ax,[-0.5,0.5]);
        zlim(ax,[-0.5,0.5]);
        title(ax,'Contact cap: smaller areas at dense nodes');
    end
end
fprintf('collocation_patch_area self-test passed: area error %.3g.\n', ...
    abs(sum(areas)-exactArea)/exactArea);
end
