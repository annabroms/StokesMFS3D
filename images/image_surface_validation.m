function [err,points,velocity,relative,gap] = ...
    image_surface_validation(q,grid,lambda,U,opt)
%IMAGE_SURFACE_VALIDATION Check image flow on independent spheres and caps.
%   image_surface_validation() checks a zero-flow one-sphere case.
%
% Each body gets an independent spherical design plus Fibonacci points on
% every contact cap. relative is pointwise error divided by local rigid
% speed, except near zero speed, where the maximum rigid speed is used.

if nargin == 0
    self_test();
    err = []; points = {}; velocity = {}; relative = {}; gap = {};
    return
end
P = size(q,1);
[fullCheck,~] = get_sphdesign(opt.validationCount);
points = cell(P,1);
velocity = cell(P,1);
relative = cell(P,1);
gap = cell(P,1);
maxSpeed = 0;
for k = 1:P
    x = q(k,:)+fullCheck;
    gapMask = false(size(x,1),1);
    for contact = 1:size(grid.body(k).contactDirections,1)
        direction = grid.body(k).contactDirections(contact,:);
        sphere = struct('axes',ones(3,1),'Q',eye(3),'c',q(k,:).');
        capOpt = struct('collocationPatchMode',"fibonacci", ...
            'collocationPatchPointCount',opt.validationGapCount, ...
            'collocationPatchAngle',max(opt.capAngles));
        [cap,~] = collocation_surface_patch( ...
            sphere,q(k,:).'+direction.',capOpt);
        x = [x;cap]; %#ok<AGROW>
        gapMask = [gapMask;true(size(cap,1),1)]; %#ok<AGROW>
    end
    points{k} = x;
    gap{k} = gapMask;
    desired = reshape(getKmat(x,q(k,:))*U(6*(k-1)+(1:6)),3,[]).';
    maxSpeed = max(maxSpeed,max(vecnorm(desired,2,2)));
end
maxSpeed = max(maxSpeed,eps);
err = 0;
coeff = cell(P,1);
for j = 1:P
    coeff{j} = image_source_coefficients(grid.body(j).basis,lambda{j});
end
for k = 1:P
    x = points{k};
    v = zeros(size(x));
    for j = 1:P
        basis = grid.body(j).basis;
        for first = 1:opt.batchSize:size(x,1)
            nodes = first:min(first+opt.batchSize-1,size(x,1));
            v(nodes,:) = v(nodes,:)+ ...
                direct_image_velocity(basis,coeff{j},x(nodes,:));
        end
    end
    desired = reshape(getKmat(x,q(k,:))*U(6*(k-1)+(1:6)),3,[]).';
    denominator = vecnorm(desired,2,2);
    denominator(denominator < 1e-12*maxSpeed) = maxSpeed;
    relative{k} = vecnorm(v-desired,2,2)./denominator;
    velocity{k} = v;
    err = max(err,max(relative{k}));
end
end

function self_test()
%SELF_TEST checks dimensions and zero-flow residual without contacts.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'geometry','spheres'), ...
    fullfile(root,'common','dense'),fullfile(root,'common','matvecs'));
body.basis = struct('family',"S",'points',[0,0,0], ...
    'normalDirections',zeros(0,3),'scale',1);
body.contactDirections = zeros(0,3);
grid.body = body;
opt = struct('validationCount',40,'validationGapCount',21, ...
    'capAngles',[pi/5,pi/60],'batchSize',16);
[err,points,velocity,relative,gap] = image_surface_validation( ...
    [0,0,0],grid,{zeros(3,1)},zeros(6,1),opt);
assert(err == 0 && ~any(gap{1}));
assert(isequal(size(points{1}),size(velocity{1})));
assert(all(relative{1} == 0));
fprintf('image_surface_validation self-test passed.\n');
end
