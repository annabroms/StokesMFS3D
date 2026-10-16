function [patch,metadata] = lap3d_collocation_patch(E,contactPoint,opt)
%LAP3D_COLLOCATION_PATCH Exponentially graded rings on one body surface.
%   [patch,metadata] = lap3d_collocation_patch(E,contactPoint,opt)
%   lap3d_collocation_patch   % self-test and one-body demonstration figure
%
% E has semi-axes axes (3-by-1), rotation Q (3-by-3), and centre c (3-by-1).
% contactPoint is a 3-by-1 surface point.  opt supplies the five controls
% collocationPatchAngle (radians), collocationPatchLayers, collocationExponent,
% collocationMinAzimuth, and collocationAzimuthOversampling, as normalized
% by lap3d_complete_grid_options.  Oversampling 1 preserves the original grid;
% 2 roughly doubles each circumference-based azimuthal count.
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

layer = (1:opt.collocationPatchLayers).';
radius = opt.collocationPatchAngle* ...
    expm1(opt.collocationExponent*layer/opt.collocationPatchLayers)/ ...
    expm1(opt.collocationExponent);
radialSpacing = diff([0;radius]);
if any(diff(radius) <= 0) || any(diff(radialSpacing) <= 0)
    error('lap3d_collocation_patch:BadPatchSpacing', ...
        'The graded collocation radii and spacings must increase strictly.');
end

circumference = 2*pi*sin(radius);
azimuthCount = max(opt.collocationMinAzimuth, ...
    ceil(opt.collocationAzimuthOversampling*circumference./radialSpacing));
azimuthCount = azimuthCount+mod(azimuthCount,2);
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
metadata.radialSpacing = radialSpacing;
metadata.azimuthCount = azimuthCount;
metadata.spacingRatio = radialSpacing(end)/radialSpacing(1);
metadata.contactDirection = contactDirection;
metadata.generatedCount = size(patch,1);
end


function self_test()
%SELF_TEST checks sphere/rotated-ellipsoid membership, counts, and grading.

projectDirectory = fileparts(fileparts(mfilename('fullpath')));
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
[patch,metadata] = lap3d_collocation_patch(E,contact,opt);
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

figure('Name','Enhancing collocation patch self-test');
ax = axes;
scatter3(ax,patch(:,1),patch(:,2),patch(:,3),14,[0.1,0.4,0.8],'filled');
hold(ax,'on')
scatter3(ax,contact(1),contact(2),contact(3),70,'k','filled');
axis(ax,'equal')
grid(ax,'on')
view(ax,3)
xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')
title(ax,sprintf('One sphere: %d patch nodes, closest point in black',size(patch,1)))

E.axes = [1.6;0.7;0.7];
angle = 0.4;
E.Q = [cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
for direction = [1,0;0,0;0,1]
    contact = E.c+E.Q*(E.axes.*direction);
    [patch,metadata] = lap3d_collocation_patch(E,contact,opt);
    localPoints = E.Q.'*(patch.'-E.c);
    assert(max(abs(sum((localPoints./E.axes).^2,1)-1)) < 100*eps);
    assert(size(patch,1) == metadata.generatedCount);
    assert(norm(patch(1,:).'-contact,inf) == 0);
end
opt.collocationAzimuthOversampling = 2;
[denserPatch,denserMetadata] = lap3d_collocation_patch(E,contact,opt);
assert(size(denserPatch,1) > size(patch,1));
assert(all(denserMetadata.azimuthCount >= metadata.azimuthCount));
fprintf(['lap3d_collocation_patch self-test passed: 535 sphere nodes; ', ...
    'rotated ellipsoid, pole, and azimuthal refinement checked.\n']);
end
