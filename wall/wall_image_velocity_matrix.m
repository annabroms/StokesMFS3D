function A = wall_image_velocity_matrix(basis,targets,kernelOptions)
%WALL_IMAGE_VELOCITY_MATRIX Dense mixed-source resistance velocity map.
%   A = wall_image_velocity_matrix(basis,targets,kernelOptions)
%   wall_image_velocity_matrix() runs a small ordering/scaling regression.
%
% Reuse the image basis fields family, points, normalDirections, scale.
% Rows and each block's columns are interleaved [x1;y1;z1;x2;...].
% All sources include their complete no-slip wall images through
% StokesHalfspaceFMM. No independent source is placed below the wall.
% Preserve the conventions in images/stokes_singularity_matrices:
%   R: tau cross r/(8*pi*mu*r^3),
%   T: -6*r*(r.q)*(r.n)/(8*pi*mu*r^5),
%   D: (-d/r^3+3*r*(r.d)/r^5)/(4*pi*mu).
% Hence the package receives -q for T and -2*d for D. Scales multiply
% fitted columns here and reconstructed coefficients exactly once.

if nargin == 0
    self_test(); A = []; return
end
required = {'blake.velocity','blake.rotlet','blake.stresslet', ...
    'blake.potential_dipole'};
for k = 1:numel(required)
    if isempty(which(required{k}))
        error('wall:MissingImageKernel', ...
            'Add StokesHalfspaceFMM/matlab to the path; missing %s.',required{k});
    end
end
validateattributes(targets,{'double'},{'real','finite','2d','ncols',3});
columns = 3*sum(arrayfun(@(block) size(block.points,1),basis));
A = zeros(3*size(targets,1),columns);
offset = 0;
for b = 1:numel(basis)
    block = basis(b);
    count = size(block.points,1);
    cols = offset+(1:3*count);
    if block.family == "S"
        options = struct('kernel',"blake",'mu',kernelOptions.mu, ...
            'wall_height',kernelOptions.wall_height);
        A(:,cols) = block.scale*stokes_kernel_matrix(block.points,targets,options);
    else
        for j = 1:count
            options = kernelOptions;
            if block.family == "T"
                options.normals = block.normalDirections(j,:).';
            end
            for component = 1:3
                strength = zeros(3,1);
                strength(component) = 1;
                switch block.family
                    case "R"
                        u = blake.rotlet(block.points(j,:).',strength,targets.',options);
                    case "T"
                        u = blake.stresslet(block.points(j,:).',-strength,targets.',options);
                    case "D"
                        u = blake.potential_dipole(block.points(j,:).',-2*strength,targets.',options);
                    otherwise
                        error('wall:ImageFamily','Unknown image family %s.',block.family);
                end
                A(:,offset+3*(j-1)+component) = block.scale*u(:);
            end
        end
    end
    offset = offset+3*count;
end
end

function self_test()
% Check block order, independent contractions, all orientations and scaling.
options = struct('method','blake','mu',1.7,'wall_height',0.2);
source = [0.1,0.2,0.9; -0.2,0.1,1.1];
targets = [0.3,-0.2,0.5; 0.4,0.2,1.8; -0.1,0.5,0.2];
template = struct('family',"S",'points',source, ...
    'normalDirections',zeros(0,3),'scale',1);
basis = repmat(template,6,1);
families = ["S","R","T","T","T","D"];
normals = [0,-1,0;-1,0,0;0,0,-1];
for b = 1:6
    basis(b).family = families(b);
    basis(b).scale = 0.5+b/3;
    if families(b) == "T"
        basis(b).normalDirections = repmat(normals(b-2,:),2,1);
    end
end
A = wall_image_velocity_matrix(basis,targets,options);
z = sin((1:size(A,2)).');
coeff = image_source_coefficients(basis,z);
u = direct_wall_image_velocity(basis,coeff,targets,options);
actionError = norm(A*z-reshape(u.',[],1));
if ~isequal(size(A),[9,36]) || actionError > 1e-12*max(1,norm(A*z))
    error('wall:ImageMatrixTest','Mixed block order/scaling check failed.');
end
% Cartesian unit columns independently expose transpose and sign mistakes.
for b = 1:6
    for j = 1:2
        f = [0.2;-0.4;0.7];
        o = options;
        switch families(b)
            case "S"
                ref = blake.velocity(source(j,:).',f,targets.',o);
            case "R"
                ref = blake.rotlet(source(j,:).',f,targets.',o);
            case "T"
                o.normals = normals(b-2,:).';
                ref = blake.stresslet(source(j,:).',-f,targets.',o);
            case "D"
                ref = blake.potential_dipole(source(j,:).',-2*f,targets.',o);
        end
        cols = (b-1)*6+3*(j-1)+(1:3);
        if norm(A(:,cols)*f-basis(b).scale*ref(:)) > 1e-12
            error('wall:ImageColumnTest','Kernel contraction failed in block %d.',b);
        end
    end
end
fprintf('wall_image_velocity_matrix PASS: 9-by-36, action error %.3e.\n',actionError);
end
