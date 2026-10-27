function traction = direct_image_traction(basis,strength,targets,normals)
%DIRECT_IMAGE_TRACTION Evaluate -sigma(u,p)*n for physical S/R/T/D sources.
% n points out of the particle. This traction gives force on the fluid,
% consistent with the source strengths and the existing SLP traction map.
% Sources, coefficients, targets and normals use laboratory coordinates.
% Coefficients are interleaved by node and scaled by basis.scale exactly
% once. The viscosity is one, matching stokes_singularity_matrices.

if nargin == 0
    self_test();
    traction = [];
    return
end
validateattributes(targets,{'numeric'}, ...
    {'real','finite','2d','ncols',3});
validateattributes(normals,{'numeric'}, ...
    {'real','finite','size',size(targets)});
coeff = image_source_coefficients(basis,strength);
traction = zeros(size(targets));
A = -6/(8*pi);
B = 1/(4*pi);
for j = 1:size(targets,1)
    m = normals(j,:);
    t = zeros(1,3);
    for blockIndex = 1:numel(basis)
        block = basis(blockIndex);
        c = coeff{blockIndex};
        r = targets(j,:)-block.points;
        r2 = sum(r.^2,2);
        assert(all(r2>0),'Source and traction target coincide.');
        invR3 = r2.^(-3/2);
        invR5 = invR3./r2;
        invR7 = invR5./r2;
        rm = r*m.';
        rc = sum(r.*c,2);
        switch block.family
            case "S"
                local = (6/(8*pi))*(rm.*rc.*invR5).*r;
            case "R"
                w = cross(c,r,2);
                wm = w*m.';
                local = (3/(8*pi))*(rm.*w+wm.*r).*invR5;
            case "T"
                n = block.normalDirections;
                assert(isequal(size(n),size(c)));
                rn = sum(r.*n,2);
                nc = sum(n.*c,2);
                v = rc.*n+rn.*c;
                vm = v*m.';
                pressure = 2*A*(rn.*rc.*invR5- ...
                    nc.*invR3/3);
                local = pressure.*m-A*( ...
                    2*(rn.*rc.*invR5).*m+ ...
                    (vm.*invR5).*r+ ...
                    (rm.*invR5).*v- ...
                    (10*rn.*rc.*rm.*invR7).*r);
            case "D"
                cm = c*m.';
                local = -B*(6*(rm.*c+cm.*r).*invR5+ ...
                    6*(rc.*invR5).*m- ...
                    (30*rc.*rm.*invR7).*r);
            otherwise
                error('direct_image_traction:UnknownFamily', ...
                    'Unknown source family.');
        end
        t = t+sum(local,1);
    end
    traction(j,:) = t;
end
end

function self_test()
%SELF_TEST compares the S branch with the established traction kernel.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(genpath(root));
basis = struct('family',"S", ...
    'points',[0.1,0,0;-0.2,0.1,0.2], ...
    'normalDirections',zeros(0,3),'scale',1);
f = [1;0.2;-0.3;-0.4;0.7;0.5];
x = [1,0.1,0.2;0.2,1,-0.3;-0.4,0.3,1];
n = x./vecnorm(x,2,2);
t = direct_image_traction(basis,f,x,n);
reference = reshape(getTractionDirectSLP( ...
    basis.points,x,n,f),3,[]).';
assert(norm(t-reference,'fro')<1e-12);
% Compare the remaining analytic gradients with central velocity differences.
source = struct('family',"R",'points',[0.1,-0.2,0.05], ...
    'normalDirections',[0.3,0.4,0.5], 'scale',1);
strength = [0.7;-0.4;0.2];
target = [1.2,0.8,-0.7];
normal = target/norm(target);
h = 1e-6;
for family = ["R","T","D"]
    source.family = family;
    gradient = zeros(3,3);
    for coordinate = 1:3
        shift = zeros(1,3); shift(coordinate) = h;
        plus = stokes_mfs_velocity_matrix(source,target+shift)*strength;
        minus = stokes_mfs_velocity_matrix(source,target-shift)*strength;
        gradient(:,coordinate) = (plus-minus)/(2*h);
    end
    pressure = 0;
    if family == "T"
        r = target-source.points;
        nc = source.normalDirections*strength;
        pressure = (-12/(8*pi))*( ...
            (r*source.normalDirections.')*(r*strength)/norm(r)^5- ...
            nc/(3*norm(r)^3));
    end
    finiteDifference = ((gradient+gradient.')- ...
        pressure*eye(3))*normal.';
    actual = direct_image_traction( ...
        source,strength,target,normal);
    assert(norm(finiteDifference+actual.')<1e-7);
end
% The enclosing sphere must recover the physical S force and R torque.
[sphere,~] = get_sphdesign(1000);
area = 4*pi/size(sphere,1);
for family = ["S","R"]
    source.family = family;
    source.points = [0,0,0];
    surfaceTraction = direct_image_traction( ...
        source,strength,sphere,sphere);
    force = area*sum(surfaceTraction,1).';
    torque = area*sum(cross(sphere,surfaceTraction,2),1).';
    if family == "S"
        assert(norm(force-strength)<1e-8);
        assert(norm(torque)<1e-8);
    else
        assert(norm(force)<1e-8);
        assert(norm(torque-strength)<1e-8);
    end
end
fprintf('direct_image_traction self-test passed.\n');
end
