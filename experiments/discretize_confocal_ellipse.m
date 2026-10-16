function [points, radialCoordinate, azimuth] = discretize_confocal_ellipse( ...
    a,b,c,nTheta,nRadial,type)
%DISCRETIZE_CONFOCAL_ELLIPSE Place nodes on an ellipsoid's focal ellipse.
%
%   points = DISCRETIZE_CONFOCAL_ELLIPSE(a,b,c,nTheta,nRadial)
%   points = DISCRETIZE_CONFOCAL_ELLIPSE(a,b,c,nTheta,nRadial,type)
%   [points,r,theta] = DISCRETIZE_CONFOCAL_ELLIPSE(...)
%
%   The ellipsoid is aligned with the coordinate axes and has semiaxes
%   (a,b,c), with c smaller than both a and b.  Its filled focal ellipse is
%
%       x^2/(a^2-c^2) + y^2/(b^2-c^2) <= 1,    z = 0.
%
%   Types 1 and 2 use nTheta equispaced angles on each of nRadial rings.
%   Their radial rings are midpoint quantiles u_j = (j-1/2)/nRadial, so no
%   node lies at the duplicated centre or exactly on the focal boundary.
%
%   type = 1 (default) uses
%
%       r = sqrt(1-(1-u)^(2/3)),
%
%   which represents equal-strength nodes with areal density proportional
%   to sqrt(1-r^2).  Type = 2 uses r = sqrt(u), giving a distribution that
%   is uniform per unit area of the ellipse.  Type = 3 uses uniformly
%   spaced midpoint radii and varies the number of angles approximately in
%   proportion to radius.  This gives a quasi-uniform node spacing on the
%   filled ellipse while retaining exactly nTheta*nRadial nodes; adjacent
%   rings are staggered by half an angular step.  Type = 4 uses the same
%   layout as type 3, scaled just enough that the convex hull of its outer
%   ring encloses the complete focal ellipse.
%
%   points is (nTheta*nRadial)-by-3.  Angle varies fastest within each
%   radial ring.  Types 1 and 2 have nTheta rows per ring, whereas types 3
%   and 4 have a radius-dependent number of rows.  r and theta are column
%   vectors in the same order as points.  Type 4 has outer radii greater
%   than one because its outer polygon circumscribes the ellipse.
%
%   Note: a target areal density proportional to r would instead require
%   the radial map r = u^(1/3).

if nargin == 0
    self_test;
    points = [];
    radialCoordinate = [];
    azimuth = [];
    return;
end

if nargin < 6 || isempty(type)
    type = 1;
end

validateattributes(a,{'numeric'}, ...
    {'real','finite','positive','scalar'},mfilename,'a',1);
validateattributes(b,{'numeric'}, ...
    {'real','finite','positive','scalar'},mfilename,'b',2);
validateattributes(c,{'numeric'}, ...
    {'real','finite','positive','scalar'},mfilename,'c',3);
validateattributes(nTheta,{'numeric'}, ...
    {'real','finite','positive','integer','scalar'},mfilename,'nTheta',4);
validateattributes(nRadial,{'numeric'}, ...
    {'real','finite','positive','integer','scalar'},mfilename,'nRadial',5);
validateattributes(type,{'numeric'}, ...
    {'real','finite','integer','scalar'},mfilename,'type',6);

if a <= c || b <= c
    error('discretize_confocal_ellipse:InvalidSemiaxes', ...
        'The x- and y-semiaxes a and b must both be greater than c.');
end
if type ~= 1 && type ~= 2 && type ~= 3 && type ~= 4
    error('discretize_confocal_ellipse:InvalidType', ...
        ['type must be 1 (weighted), 2 (uniform area), ' ...
         '3 (quasi-uniform), or 4 (enclosing quasi-uniform).']);
end

u = ((1:nRadial)-0.5)/nRadial;
if type == 1
    ringRadii = sqrt(1-(1-u).^(2/3));
elseif type == 2
    ringRadii = sqrt(u);
else
    ringRadii = u;
end

if type == 3 || type == 4
    totalNodeCount = nTheta*nRadial;
    idealRingNodeCounts = totalNodeCount*ringRadii/sum(ringRadii);
    ringNodeCounts = max(1,floor(idealRingNodeCounts));

    % Assign rounding remainders without changing the requested total count.
    nodesLeft = totalNodeCount-sum(ringNodeCounts);
    roundingRemainders = idealRingNodeCounts-ringNodeCounts;
    [~,allocationOrder] = sort(roundingRemainders,'descend');
    ringNodeCounts(allocationOrder(1:nodesLeft)) = ...
        ringNodeCounts(allocationOrder(1:nodesLeft))+1;

    radialCoordinate = zeros(totalNodeCount,1);
    azimuth = zeros(totalNodeCount,1);
    nodeOffset = 0;
    for ringIndex = 1:nRadial
        numberOfAngles = ringNodeCounts(ringIndex);
        ringIndices = nodeOffset+(1:numberOfAngles);
        angleOffset = mod(ringIndex-1,2)*pi/numberOfAngles;
        radialCoordinate(ringIndices) = ringRadii(ringIndex);
        azimuth(ringIndices) = angleOffset + ...
            (0:numberOfAngles-1)'*(2*pi/numberOfAngles);
        nodeOffset = nodeOffset+numberOfAngles;
    end

    if type == 4
        outerRingNodeCount = ringNodeCounts(end);
        if outerRingNodeCount < 3
            error('discretize_confocal_ellipse:InsufficientResolution', ...
                ['type 4 requires at least three nodes on the outer ring ' ...
                 'to enclose the focal ellipse.']);
        end

        % A regular polygon of radius sec(pi/N) encloses the unit disc.
        outerRingRadius = ringRadii(end);
        enclosureScale = 1/(outerRingRadius*cos(pi/outerRingNodeCount));
        radialCoordinate = enclosureScale*radialCoordinate;
    end
else
    ringAngles = (0:nTheta-1)*(2*pi/nTheta);
    [angleGrid,radiusGrid] = ndgrid(ringAngles,ringRadii);
    radialCoordinate = radiusGrid(:);
    azimuth = angleGrid(:);
end

focalSemiAxisX = sqrt(a^2-c^2);
focalSemiAxisY = sqrt(b^2-c^2);
points = [focalSemiAxisX*radialCoordinate.*cos(azimuth), ...
          focalSemiAxisY*radialCoordinate.*sin(azimuth), ...
          zeros(numel(radialCoordinate),1)];
end

function self_test
%SELF_TEST Check the four distributions and visualise their focal-ellipse nodes.
a = 1.2;
b = 1.1;
c = 0.7;
nTheta = 12;
nRadial = 1000;

[weightedPoints,weightedRadii] = discretize_confocal_ellipse( ...
    a,b,c,nTheta,nRadial,1);
[uniformPoints,uniformRadii] = discretize_confocal_ellipse( ...
    a,b,c,nTheta,nRadial,2);
[quasiUniformPoints,quasiUniformRadii] = discretize_confocal_ellipse( ...
    a,b,c,nTheta,nRadial,3);
[enclosingPoints,enclosingRadii,enclosingAngles] = ...
    discretize_confocal_ellipse(a,b,c,nTheta,nRadial,4);

assert(isequal(size(weightedPoints),[nTheta*nRadial,3]));
assert(isequal(size(uniformPoints),[nTheta*nRadial,3]));
assert(isequal(size(quasiUniformPoints),[nTheta*nRadial,3]));
assert(isequal(size(enclosingPoints),[nTheta*nRadial,3]));
assert(all(weightedPoints(:,3) == 0));
assert(all(uniformPoints(:,3) == 0));
assert(all(quasiUniformPoints(:,3) == 0));
assert(all(enclosingPoints(:,3) == 0));

weightedLevel = weightedPoints(:,1).^2/(a^2-c^2) + ...
    weightedPoints(:,2).^2/(b^2-c^2);
uniformLevel = uniformPoints(:,1).^2/(a^2-c^2) + ...
    uniformPoints(:,2).^2/(b^2-c^2);
quasiUniformLevel = quasiUniformPoints(:,1).^2/(a^2-c^2) + ...
    quasiUniformPoints(:,2).^2/(b^2-c^2);
enclosingLevel = enclosingPoints(:,1).^2/(a^2-c^2) + ...
    enclosingPoints(:,2).^2/(b^2-c^2);
assert(max(abs(weightedLevel-weightedRadii.^2)) < 20*eps);
assert(max(abs(uniformLevel-uniformRadii.^2)) < 20*eps);
assert(max(abs(quasiUniformLevel-quasiUniformRadii.^2)) < 20*eps);
assert(max(abs(enclosingLevel-enclosingRadii.^2)) < 20*eps);
assert(all(weightedLevel > 0 & weightedLevel < 1));
assert(all(uniformLevel > 0 & uniformLevel < 1));
assert(all(quasiUniformLevel > 0 & quasiUniformLevel < 1));
assert(all(enclosingLevel > 0));

% The exact second radial moments are 2/5, 1/2, and approximately 1/2.
assert(abs(mean(weightedRadii.^2)-2/5) < 1e-5);
assert(abs(mean(uniformRadii.^2)-1/2) < 1e-12);
assert(abs(mean(quasiUniformRadii.^2)-1/2) < 1e-3);

uniqueQuasiUniformRadii = unique(quasiUniformRadii);
ringNodeCounts = arrayfun( ...
    @(radius) sum(quasiUniformRadii == radius),uniqueQuasiUniformRadii);
idealRingNodeCounts = 2*nTheta*uniqueQuasiUniformRadii;
assert(numel(uniqueQuasiUniformRadii) == nRadial);
assert(all(diff(ringNodeCounts) >= 0));
assert(max(abs(ringNodeCounts-idealRingNodeCounts)) <= 1);

outerRingRadius = max(enclosingRadii);
outerRingAngles = enclosingAngles(enclosingRadii == outerRingRadius);
outerRingNodeCount = numel(outerRingAngles);
assert(outerRingNodeCount >= 3);
assert(abs(outerRingRadius*cos(pi/outerRingNodeCount)-1) < 100*eps);

visualNTheta = 36;
visualNRadial = 18;
[weightedPlotPoints,weightedPlotRadii] = discretize_confocal_ellipse( ...
    a,b,c,visualNTheta,visualNRadial,1);
[uniformPlotPoints,uniformPlotRadii] = discretize_confocal_ellipse( ...
    a,b,c,visualNTheta,visualNRadial,2);
[quasiUniformPlotPoints,quasiUniformPlotRadii] = ...
    discretize_confocal_ellipse( ...
        a,b,c,visualNTheta,visualNRadial,3);
[enclosingPlotPoints,enclosingPlotRadii] = discretize_confocal_ellipse( ...
    a,b,c,visualNTheta,visualNRadial,4);

pointSets = {weightedPlotPoints,uniformPlotPoints,quasiUniformPlotPoints, ...
    enclosingPlotPoints};
radiusSets = {weightedPlotRadii,uniformPlotRadii,quasiUniformPlotRadii, ...
    enclosingPlotRadii};
panelTitles = {'type 1: weighted','type 2: uniform area', ...
    'type 3: quasi-uniform','type 4: enclosing quasi-uniform'};
boundaryAngle = linspace(0,2*pi,400);
focalSemiAxisX = sqrt(a^2-c^2);
focalSemiAxisY = sqrt(b^2-c^2);
boundaryX = focalSemiAxisX*cos(boundaryAngle);
boundaryY = focalSemiAxisY*sin(boundaryAngle);

figure('Name','Confocal ellipse discretizations');
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for typeIndex = 1:4
    ax = nexttile(layout);
    pointsToPlot = pointSets{typeIndex};
    radiiToPlot = radiusSets{typeIndex};
    scatter(ax,pointsToPlot(:,1),pointsToPlot(:,2),18, ...
        radiiToPlot,'filled');
    hold(ax,'on');
    plot(ax,boundaryX,boundaryY,'k-','LineWidth',1.5);
    axis(ax,'equal');
    grid(ax,'on');
    xlabel(ax,'x');
    ylabel(ax,'y');
    title(ax,panelTitles{typeIndex});
    clim(ax,[0,max(1,max(radiiToPlot))]);
    colorbar(ax);
end
colormap(turbo);
title(layout,sprintf( ...
    'Focal ellipse for ellipsoid semiaxes (%.1f, %.1f, %.1f)',a,b,c));
fprintf('discretize_confocal_ellipse self-test passed.\n');
end
