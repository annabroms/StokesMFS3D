function areas = lap3d_collocation_surface_elements(points,bodyIndices,bodies)
%LAP3D_COLLOCATION_SURFACE_ELEMENTS Estimate physical area per surface node.
%
% Use the COMPLETE, deduplicated collocation grid on each ellipsoid/sphere,
% including any merged contact patch. Map it to the unit parameter sphere,
% triangulate, and share each spherical triangle's solid angle equally among
% its three vertices. Multiply by the ellipsoid map's surface Jacobian.
% Thus denser patches receive smaller area per node without counting their
% overlap with the basic grid twice. These are positive lumped area estimates,
% not the high-order quadrature weights of the original unrefined grid.
% For a sphere of radius R, their sum is exactly 4*pi*R^2 (up to roundoff).
% No source positions or source quadrature weights are involved.

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
