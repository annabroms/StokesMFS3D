function t = image_traction(grid,lambda,x,n,mu)
%IMAGE_TRACTION Exterior stress times supplied target normals for image MFS.
%   t = image_traction(grid,lambda,x,n,mu)
%   image_traction() checks all four families against velocity derivatives.
%
% grid and lambda are details.grid and details.lambda from an image solve.
% lambda contains unit-viscosity velocity coefficients, interleaved by node.
% Thus the physical traction is mu times the unit-viscosity stress.
% On a unit sphere, n=x-centre points from the solid into the fluid.

if nargin == 0
    self_test();
    t = [];
    return
end
validateattributes(x,{'numeric'},{'real','finite','2d','ncols',3});
validateattributes(n,{'numeric'},{'real','finite','2d','ncols',3});
validateattributes(mu,{'numeric'},{'real','finite','scalar','positive'});
if ~isequal(size(x),size(n)) || ...
        any(abs(vecnorm(n,2,2)-1) > 1e-10)
    error('image_traction:BadNormals', ...
        'Supply one unit normal per target.');
end
if ~iscell(lambda) || numel(lambda) ~= numel(grid.body)
    error('image_traction:BadStrengths', ...
        'lambda must contain one coefficient vector per body.');
end

t = zeros(size(x));
batchSize = 128;
for body = 1:numel(grid.body)
    basis = grid.body(body).basis;
    coeff = image_source_coefficients(basis,lambda{body});
    for block = 1:numel(basis)
        y = basis(block).points;
        c = coeff{block};
        if basis(block).family == "T" && ...
                ~isequal(size(basis(block).normalDirections),size(c))
            error('image_traction:BadDirection', ...
                'Each stresslet needs one source direction.');
        end
        for first = 1:batchSize:size(x,1)
            rows = first:min(first+batchSize-1,size(x,1));
            z = x(rows,:);
            nn = n(rows,:);
            r = cell(1,3);
            for j = 1:3
                r{j} = z(:,j)-y(:,j).';
            end
            r2 = r{1}.^2+r{2}.^2+r{3}.^2;
            if any(r2 == 0,'all')
                error('image_traction:CoincidentPoints', ...
                    'A target coincides with a source.');
            end
            ir3 = r2.^(-3/2);
            ir5 = ir3./r2;
            rn = r{1}.*nn(:,1)+r{2}.*nn(:,2)+r{3}.*nn(:,3);
            rc = r{1}.*c(:,1).'+r{2}.*c(:,2).'+r{3}.*c(:,3).';
            switch basis(block).family
                case "S"
                    scalar = (-3/(4*pi))*rc.*rn.*ir5;
                    for j = 1:3
                        t(rows,j) = t(rows,j)+sum(scalar.*r{j},2);
                    end
                case "R"
                    w = {c(:,2).'.*r{3}-c(:,3).'.*r{2}, ...
                        c(:,3).'.*r{1}-c(:,1).'.*r{3}, ...
                        c(:,1).'.*r{2}-c(:,2).'.*r{1}};
                    wn = w{1}.*nn(:,1)+w{2}.*nn(:,2)+ ...
                        w{3}.*nn(:,3);
                    for j = 1:3
                        t(rows,j) = t(rows,j)-3/(8*pi)*sum( ...
                            (w{j}.*rn+r{j}.*wn).*ir5,2);
                    end
                case "T"
                    q = basis(block).normalDirections;
                    a = r{1}.*q(:,1).'+r{2}.*q(:,2).'+ ...
                        r{3}.*q(:,3).';
                    qc = sum(q.*c,2).';
                    v = cell(1,3);
                    for j = 1:3
                        v{j} = q(:,j).'.*rc+c(:,j).'.*a;
                    end
                    vn = v{1}.*nn(:,1)+v{2}.*nn(:,2)+ ...
                        v{3}.*nn(:,3);
                    ir7 = ir5./r2;
                    for j = 1:3
                        value = -qc.*nn(:,j).*ir3/(2*pi) ...
                            -3*(r{j}.*vn+v{j}.*rn).*ir5/(4*pi) ...
                            +15*a.*rc.*r{j}.*rn.*ir7/(2*pi);
                        t(rows,j) = t(rows,j)+sum(value,2);
                    end
                case "D"
                    cn = nn*c.';
                    ir7 = ir5./r2;
                    for j = 1:3
                        value = 3*(c(:,j).'.*rn+r{j}.*cn+ ...
                            rc.*nn(:,j)).*ir5/(2*pi) ...
                            -15*rc.*r{j}.*rn.*ir7/(2*pi);
                        t(rows,j) = t(rows,j)+sum(value,2);
                    end
                otherwise
                    error('image_traction:BadFamily', ...
                        'Unknown source family "%s".',basis(block).family);
            end
        end
    end
end
t = mu*t;
end

function self_test()
%SELF_TEST checks families, mixtures, normal sign, viscosity, and PDE data.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','eval'),fullfile(root,'common','dense'), ...
    fullfile(root,'drops'),fullfile(root,'images'));
base = struct('family',"S",'points',[0,0,0], ...
    'normalDirections',zeros(0,3),'scale',1);
families = ["S","R","T","D"];
x = [1.1,0.4,-0.3;0.5,-0.8,1.2];
n = x./vecnorm(x,2,2);
strength = [0.7;-0.2;0.4];
muValues = [0.6,1.7,2.4];
maxError = 0;
maxPde = 0;
maxStokes = 0;
for family = families
    block = base;
    block.family = family;
    if family == "T"
        block.normalDirections = [0.3,0.4,sqrt(0.75)];
    end
    grid.body.basis = block;
    for mu = muValues
        actual = image_traction(grid,{strength},x,n,mu);
        flipped = image_traction(grid,{strength},x,-n,mu);
        maxError = max(maxError,max(abs(actual+flipped),[],'all'));
        for target = 1:size(x,1)
            point = x(target,:);
            step = 1e-5;
            gradient = zeros(3);
            for j = 1:3
                shift = zeros(1,3); shift(j) = step;
                up = stokes_mfs_velocity_matrix(block,point+shift)*strength;
                um = stokes_mfs_velocity_matrix(block,point-shift)*strength;
                gradient(:,j) = (up-um)/(2*step);
            end
            p = image_pressure(grid,{strength},point,mu);
            expected = (-p*eye(3)+mu*(gradient+gradient.'))*n(target,:).';
            maxError = max(maxError,norm(actual(target,:).'-expected));
            maxPde = max(maxPde,abs(trace(gradient)));
            step = 2e-4;
            lap = zeros(1,3);
            gp = zeros(1,3);
            u0 = stokes_mfs_velocity_matrix(block,point)*strength;
            for j = 1:3
                shift = zeros(1,3); shift(j) = step;
                up = stokes_mfs_velocity_matrix(block,point+shift)*strength;
                um = stokes_mfs_velocity_matrix(block,point-shift)*strength;
                lap = lap+(up-2*u0+um).' /step^2;
                gp(j) = (image_pressure(grid,{strength},point+shift,mu)- ...
                    image_pressure(grid,{strength},point-shift,mu))/(2*step);
            end
            maxStokes = max(maxStokes,norm(-mu*lap+gp));
        end
    end
end

% A scaled mixture on two bodies must equal the sum of its ordered blocks.
mix = repmat(base,4,1);
raw = cell(4,1);
for j = 1:4
    mix(j).family = families(j);
    mix(j).points = [0.1*j,-0.03*j,0.02*j];
    mix(j).scale = 0.5+0.2*j;
    if j == 3
        mix(j).normalDirections = [0,1,0];
    end
    raw{j} = sin((1:3).'+j);
end
grid.body = repmat(struct('basis',base),2,1);
grid.body(1).basis = mix(1:2);
grid.body(2).basis = mix(3:4);
lambda = {cell2mat(raw(1:2)),cell2mat(raw(3:4))};
actual = image_traction(grid,lambda,x,n,1.7);
expected = zeros(size(x));
for j = 1:4
    one.body.basis = mix(j);
    expected = expected+image_traction(one,raw(j),x,n,1.7);
end
eMix = max(abs(actual-expected),[],'all');
one.body.basis = base;
reference = reshape(drop_traction_block(base.points,x,n)*strength,3,[]).';
eS = max(abs(image_traction(one,{strength},x,n,1)-reference),[],'all');
fprintf(['image_traction self-test: derivative %.3e, mixture %.3e, ', ...
    'S sign %.3e, divergence %.3e, PDE %.3e; traction scale %.3e\n'], ...
    maxError,eMix,eS,maxPde,maxStokes,max(abs(actual),[],'all'));
if max([maxError,eMix,eS,maxPde]) > 1e-7*max(1,max(abs(actual),[],'all')) || ...
        maxStokes > 1e-6*max(1,max(abs(actual),[],'all'))
    error('image_traction:SelfTestFailed','Traction self-test failed.');
end
end
