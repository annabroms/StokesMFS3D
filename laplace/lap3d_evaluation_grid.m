function [points,bodyIndex,nearContact,metadata] = lap3d_evaluation_grid(pair,opt)
%LAP3D_EVALUATION_GRID Independent full-body and graded near-contact test nodes.
% Shared by the Laplace and Stokes grid builder.  These points never enter
% the least-squares matrices.  The full grid uses repository quadrature;
% optional contact rings use different radii/grading from collocation.
% Near contact means inside collocationPatchAngle on the parameter sphere.

bodies = {pair.E1,pair.E2};
contacts = {pair.p1,pair.p2};
pointBlocks = cell(2,1);
nearBlocks = cell(2,1);
metadata.fullCounts = zeros(2,1);
metadata.patchCounts = zeros(2,1);
metadata.patchMetadata = cell(2,1);
testOpt = opt;
testOpt.collocationPatchLayers = opt.testPatchLayers;
testOpt.collocationExponent = opt.collocationExponent+0.5;
testOpt.collocationPatchAngle = 0.97*opt.collocationPatchAngle;
testOpt.collocationMinAzimuth = opt.collocationMinAzimuth+2;
testOpt.collocationAzimuthOversampling = opt.testPatchAzimuthOversampling;

for body = 1:2
    E = bodies{body};
    reference = ellipsoid_param(E.axes(1),E.axes(2),E.axes(3));
    reference = setupsurfquad(reference,opt.testResolution);
    fullPoints = (E.c+E.Q*reference.x).';
    patch = zeros(0,3);
    metadata.fullCounts(body) = size(fullPoints,1);
    if opt.testPatchLayers > 0
        [patch,metadata.patchMetadata{body}] = ...
            lap3d_collocation_patch(E,contacts{body},testOpt);
    end
    % Include the contact point itself in the independent error maximum,
    % even though this one point also belongs to the fitted patch.
    pointBlocks{body} = [fullPoints;patch];
    metadata.patchCounts(body) = size(patch,1);
    directions = (E.Q.'*(pointBlocks{body}.'-E.c))./E.axes;
    contactDirection = (E.Q.'*(contacts{body}-E.c))./E.axes;
    contactDirection = contactDirection/norm(contactDirection);
    nearBlocks{body} = (contactDirection.'*directions).' >= ...
        cos(opt.collocationPatchAngle)-100*eps;
end
metadata.counts = metadata.fullCounts+metadata.patchCounts;
points = vertcat(pointBlocks{:});
bodyIndex = [ones(metadata.counts(1),1);2*ones(metadata.counts(2),1)];
nearContact = vertcat(nearBlocks{:});
end
