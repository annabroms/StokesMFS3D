function [patch,metadata] = collocation_surface_patch(E,contactPoint,opt)
%COLLOCATION_SURFACE_PATCH Contact-facing cap on an ellipsoid or sphere.
%   [patch,metadata] = collocation_surface_patch(E,contactPoint,opt)
%   collocation_surface_patch   % self-test and demonstration figure
%
% E has semi-axes axes (3-by-1), rotation Q (3-by-3), and centre c (3-by-1).
% A sphere is the special case E.axes = radius*ones(3,1). CONTACTPOINT is
% a 3-by-1 surface point. OPT supplies collocationPatchAngle (radians),
% collocationPatchLayers, collocationExponent, collocationMinAzimuth, and
% collocationAzimuthOversampling. Optional collocationPatchMode is
% "exponential" (default) or "fibonacci". Optional
% collocationPatchPointCount requests an exact number of nodes. Exponential
% mode uses graded radial rings; with an exact count, each ring has at least
% collocationMinAzimuth nodes, so request enough points for the chosen layers.
% patch is N-by-3.  N = 1 + sum(metadata.azimuthCount), before the caller
% removes near-duplicates against its basic collocation grid.

if nargin == 0
    self_test;
    patch = zeros(0,3);
    metadata = struct();
    return
end

surface = ellipsoid_param(E.axes(1),E.axes(2),E.axes(3));
contactDirection = (E.Q.'*(contactPoint-E.c))./E.axes;
contactDirection = contactDirection/norm(contactDirection);

% Tangent frame on the parameter sphere; also works for a contact at a pole.
trialDirection = [0;0;1];
tangent1 = trialDirection-contactDirection*(contactDirection.'*trialDirection);
if norm(tangent1) < 1e-12
    trialDirection = [0;1;0];
    tangent1 = trialDirection-contactDirection*(contactDirection.'*trialDirection);
end
tangent1 = tangent1/norm(tangent1);
tangent2 = cross(contactDirection,tangent1);
tangent2 = tangent2/norm(tangent2);

mode = "exponential";
if isfield(opt,'collocationPatchMode')
    mode = string(opt.collocationPatchMode);
end
if mode == "fibonacci"
    assert(isfield(opt,'collocationPatchPointCount'), ...
        'Fibonacci caps require collocationPatchPointCount.');
    count = opt.collocationPatchPointCount;
    validateattributes(count,{'numeric'}, ...
        {'scalar','integer','positive'},mfilename,'collocationPatchPointCount');
    node = (1:count).';
    polarAngle = acos(1-node*(1-cos(opt.collocationPatchAngle))/count);
    goldenRatio = (1+sqrt(5))/2;
    azimuth = 2*pi*mod(node/goldenRatio,1);
    directions = contactDirection*cos(polarAngle).' + ...
        tangent1*(sin(polarAngle).*cos(azimuth)).' + ...
        tangent2*(sin(polarAngle).*sin(azimuth)).';
    localPoints = surface.Z(atan2(directions(2,:),directions(1,:)), ...
        min(1,max(-1,directions(3,:))));
    patch = (E.c+E.Q*localPoints).';
    metadata = struct('mode',mode,'angularRadius',polarAngle, ...
        'generatedCount',count,'contactDirection',contactDirection);
    return
elseif mode ~= "exponential"
    error('collocation_surface_patch:BadMode', ...
        'collocationPatchMode must be exponential or fibonacci.');
end

layer = (1:opt.collocationPatchLayers).';
radius = opt.collocationPatchAngle* ...
    expm1(opt.collocationExponent*layer/opt.collocationPatchLayers)/ ...
    expm1(opt.collocationExponent);
radialSpacing = diff([0;radius]);
if any(diff(radius) <= 0) || any(diff(radialSpacing) <= 0)
    error('collocation_surface_patch:BadPatchSpacing', ...
        'The graded collocation radii and spacings must increase strictly.');
end

circumference = 2*pi*sin(radius);
azimuthCount = max(opt.collocationMinAzimuth, ...
    ceil(opt.collocationAzimuthOversampling*circumference./radialSpacing));
azimuthCount = azimuthCount+mod(azimuthCount,2);
if isfield(opt,'collocationPatchPointCount') && ...
        ~isempty(opt.collocationPatchPointCount)
    count = opt.collocationPatchPointCount;
    minAz = max(3,ceil(opt.collocationMinAzimuth));
    validateattributes(count,{'numeric'}, ...
        {'scalar','integer','>=',1+minAz*numel(layer)},mfilename, ...
        'collocationPatchPointCount');
    azimuthCount = max(minAz, ...
        round(azimuthCount*(count-1)/sum(azimuthCount)));
    while sum(azimuthCount) < count-1
        [~,index] = max(azimuthCount);
        azimuthCount(index) = azimuthCount(index)+1;
    end
    while sum(azimuthCount) > count-1
        eligible = find(azimuthCount > minAz);
        [~,index] = max(azimuthCount(eligible));
        azimuthCount(eligible(index)) = azimuthCount(eligible(index))-1;
    end
end
patch = zeros(1+sum(azimuthCount),3);
patch(1,:) = contactPoint.';
offset = 1;
for k = 1:opt.collocationPatchLayers
    azimuth = 2*pi*(0:azimuthCount(k)-1)/azimuthCount(k);
    ringDirections = contactDirection*cos(radius(k)) + ...
        tangent1*(sin(radius(k))*cos(azimuth)) + ...
        tangent2*(sin(radius(k))*sin(azimuth));
    surfaceAzimuth = atan2(ringDirections(2,:),ringDirections(1,:));
    surfaceHeight = min(1,max(-1,ringDirections(3,:)));
    localPoints = surface.Z(surfaceAzimuth,surfaceHeight);
    rows = offset+(1:azimuthCount(k));
    patch(rows,:) = (E.c+E.Q*localPoints).';
    offset = offset+azimuthCount(k);
end

metadata.angularRadius = radius;
metadata.mode = mode;
metadata.radialSpacing = radialSpacing;
metadata.azimuthCount = azimuthCount;
metadata.spacingRatio = radialSpacing(end)/radialSpacing(1);
metadata.contactDirection = contactDirection;
metadata.generatedCount = size(patch,1);
end


function self_test()
%SELF_TEST checks sphere/rotated-ellipsoid membership, counts, and grading.

projectDirectory = fileparts(fileparts(fileparts(mfilename('fullpath'))));
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(projectDirectory,'geometry','ellipsoids'));
opt.collocationPatchLayers = 18;
opt.collocationExponent = 3;
opt.collocationPatchAngle = 0.5;
opt.collocationMinAzimuth = 8;
opt.collocationAzimuthOversampling = 1;

E.axes = 0.63*ones(3,1);
E.Q = eye(3);
E.c = [0.2;-0.1;0.3];
contact = E.c+E.Q*[E.axes(1);0;0];
[patch,metadata] = collocation_surface_patch(E,contact,opt);
assert(size(patch,1) == 535);
assert(size(patch,1) == 1+sum(metadata.azimuthCount));
assert(all(mod(metadata.azimuthCount,2) == 0));
assert(all(metadata.azimuthCount >= opt.collocationMinAzimuth));
assert(all(diff(metadata.radialSpacing) > 0));
assert(abs(metadata.angularRadius(end)-opt.collocationPatchAngle) < 10*eps);
localPoints = E.Q.'*(patch.'-E.c);
assert(max(abs(sum((localPoints./E.axes).^2,1)-1)) < 100*eps);
assert(all(isfinite(patch),'all'));
assert(size(unique(patch,'rows'),1) == size(patch,1));

figure('Name','Exponentially clustered spherical cap');
ax = subplot(1,2,1);
[sx,sy,sz] = sphere(30);
surf(ax,E.c(1)+E.axes(1)*sx,E.c(2)+E.axes(1)*sy, ...
    E.c(3)+E.axes(1)*sz,'FaceAlpha',0.08,'EdgeColor','none');
hold(ax,'on')
scatter3(ax,patch(:,1),patch(:,2),patch(:,3),9,[0.1,0.4,0.8],'filled');
scatter3(ax,contact(1),contact(2),contact(3),70,'k','filled');
boundaryAngle = linspace(0,2*pi,200);
boundary = E.c+E.axes(1)*(metadata.contactDirection*cos(opt.collocationPatchAngle)+ ...
    [zeros(1,200);cos(boundaryAngle);sin(boundaryAngle)]*sin(opt.collocationPatchAngle));
plot3(ax,boundary(1,:),boundary(2,:),boundary(3,:),'r-','LineWidth',1.5);
axis(ax,'equal')
grid(ax,'on')
view(ax,3)
xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')
title(ax,sprintf('One sphere: %d patch nodes, closest point in black',size(patch,1)))
zoomAxis = subplot(1,2,2);
cosineDistance = sum((patch-E.c.').*metadata.contactDirection.',2)/E.axes(1);
nearContact = acos(max(-1,min(1,cosineDistance))) <= ...
    opt.collocationPatchAngle/4;
scatter(zoomAxis,patch(nearContact,2)-E.c(2), ...
    patch(nearContact,3)-E.c(3),14, ...
    acos(max(-1,min(1,cosineDistance(nearContact)))),'filled');
hold(zoomAxis,'on');
scatter(zoomAxis,0,0,70,'k','filled');
axis(zoomAxis,'equal'); grid(zoomAxis,'on'); colorbar(zoomAxis);
xlabel(zoomAxis,'Tangent y offset');
ylabel(zoomAxis,'Tangent z offset');
title(zoomAxis,'Close-up of clustered region');

E.axes = [1.6;0.7;0.7];
angle = 0.4;
E.Q = [cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
for direction = [1,0;0,0;0,1]
    contact = E.c+E.Q*(E.axes.*direction);
    [patch,metadata] = collocation_surface_patch(E,contact,opt);
    localPoints = E.Q.'*(patch.'-E.c);
    assert(max(abs(sum((localPoints./E.axes).^2,1)-1)) < 100*eps);
    assert(size(patch,1) == metadata.generatedCount);
    assert(norm(patch(1,:).'-contact,inf) == 0);
end
opt.collocationAzimuthOversampling = 2;
[denserPatch,denserMetadata] = collocation_surface_patch(E,contact,opt);
assert(size(denserPatch,1) > size(patch,1));
assert(all(denserMetadata.azimuthCount >= metadata.azimuthCount));
fprintf(['collocation_surface_patch self-test passed: 535 sphere nodes; ', ...
    'rotated ellipsoid, pole, and azimuthal refinement checked.\n']);
opt.collocationPatchMode = "fibonacci";
opt.collocationPatchPointCount = 120;
[uniformCap,uniformMetadata] = collocation_surface_patch(E,contact,opt);
assert(size(uniformCap,1) == 120 && uniformMetadata.mode == "fibonacci");
end
