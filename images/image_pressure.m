function p = image_pressure(grid,lambda,x,mu)
%IMAGE_PRESSURE Exterior pressure from the complete mixed image representation.
%   p = image_pressure(grid,lambda,x,mu)
%   image_pressure() checks individual families and a scaled mixture.
%
% grid and lambda are details.grid and details.lambda from an image solve.
% lambda contains unit-viscosity velocity coefficients, interleaved by node.
% S/R/T/Q_s/S_y physical strengths are mu times these coefficients;
% D is unchanged.
% The pressure gauge is p -> 0 at infinity. R and D have zero pressure.

if nargin == 0
    self_test();
    p = [];
    return
end
validateattributes(x,{'numeric'},{'real','finite','2d','ncols',3});
validateattributes(mu,{'numeric'},{'real','finite','scalar','positive'});
if ~iscell(lambda) || numel(lambda) ~= numel(grid.body)
    error('image_pressure:BadStrengths', ...
        'lambda must contain one coefficient vector per body.');
end

p = zeros(size(x,1),1);
batchSize = 128;
for body = 1:numel(grid.body)
    basis = grid.body(body).basis;
    coeff = image_source_coefficients(basis,lambda{body});
    for block = 1:numel(basis)
        y = basis(block).points;
        c = coeff{block};
        for first = 1:batchSize:size(x,1)
            rows = first:min(first+batchSize-1,size(x,1));
            z = x(rows,:);
            switch basis(block).family
                case "S"
                    p(rows) = p(rows)+getPressure( ...
                        reshape(c.',[],1),y,z);
                case "S_y"
                    K = stokes_directional_stokeslet_matrix( ...
                        y,z,basis(block).normalDirections, ...
                        basis(block).derivativeDistances,"pressure");
                    p(rows) = p(rows)+K*reshape(c.',[],1);
                case {"T","T_n","T_n2","T_t"}
                    q = basis(block).normalDirections;
                    if ~isequal(size(q),size(c))
                        error('image_pressure:BadDirection', ...
                            'Each stresslet needs one source direction.');
                    end
                    rx = z(:,1)-y(:,1).';
                    ry = z(:,2)-y(:,2).';
                    rz = z(:,3)-y(:,3).';
                    r2 = rx.^2+ry.^2+rz.^2;
                    if any(r2 == 0,'all')
                        error('image_pressure:CoincidentPoints', ...
                            'A target coincides with a source.');
                    end
                    a = rx.*q(:,1).'+ry.*q(:,2).'+rz.*q(:,3).';
                    b = rx.*c(:,1).'+ry.*c(:,2).'+rz.*c(:,3).';
                    qc = sum(q.*c,2).';
                    p(rows) = p(rows)+sum( ...
                        (qc./r2.^(3/2)-3*a.*b./r2.^(5/2))/(2*pi),2);
                case "Q_s"
                    p(rows) = p(rows)+stokes_axial_quadrupole_field( ...
                        y,z,basis(block).normalDirections,c,"pressure");
                case {"R","D"}
                    % Both velocity fields are harmonic outside their sources.
                    if any((z(:,1)-y(:,1).').^2+ ...
                            (z(:,2)-y(:,2).').^2+ ...
                            (z(:,3)-y(:,3).').^2 == 0,'all')
                        error('image_pressure:CoincidentPoints', ...
                            'A target coincides with a source.');
                    end
                otherwise
                    error('image_pressure:BadFamily', ...
                        'Unknown source family "%s".',basis(block).family);
            end
        end
    end
end
p = mu*p;
end

function self_test()
%SELF_TEST checks pressure values, zero families, scaling, and Stokes PDE.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','eval'),fullfile(root,'images'));
base = struct('family',"S",'points',[0,0,0], ...
    'normalDirections',zeros(0,3),'scale',1);
x = [1.2,0.4,-0.3;0.5,-0.8,1.1];
f = [0.7;-0.2;0.4];
grid.body.basis = base;
p = image_pressure(grid,{f},x,1.7);
expected = 1.7*getPressure(f,base.points,x);
eS = max(abs(p-expected));
base.family = "R";
grid.body.basis = base;
eR = max(abs(image_pressure(grid,{f},x,1.7)));
base.family = "D";
grid.body.basis = base;
eD = max(abs(image_pressure(grid,{f},x,1.7)));
base.family = "T";
base.normalDirections = [0.3,0.4,sqrt(0.75)];
grid.body.basis = base;
q = base.normalDirections;
r = x-base.points;
a = r*q.';
b = r*f;
rho = vecnorm(r,2,2);
expected = 1.7*((q*f)./rho.^3-3*a.*b./rho.^5)/(2*pi);
eT = 0;
for family = ["T","T_n","T_t","T_n2"]
    grid.body.basis.family = family;
    eT = max(eT,max(abs(image_pressure(grid,{f},x,1.7)-expected)));
end
base.family = "Q_s";
base.normalDirections = [1 0 0];
grid.body.basis = base;
step = 1e-3;
shift = step*base.normalDirections;
pPlus = getPressure(f,base.points+shift,x);
pMinus = getPressure(f,base.points-shift,x);
pZero = getPressure(f,base.points,x);
qExpected = 1.7*(pPlus-2*pZero+pMinus)/step^2;
eQ = max(abs(image_pressure(grid,{f},x,1.7)-qExpected));
base.family = "T";
base.normalDirections = q;
grid.body.basis = base;

% The normalized velocity and physical pressure satisfy Stokes' equation.
coeff = image_source_coefficients(base,f);
step = 2e-4;
point = x(1,:);
lap = zeros(1,3);
gp = zeros(1,3);
for j = 1:3
    shift = zeros(1,3); shift(j) = step;
    up = direct_image_velocity(base,coeff,point+shift);
    um = direct_image_velocity(base,coeff,point-shift);
    u0 = direct_image_velocity(base,coeff,point);
    lap = lap+(up-2*u0+um)/step^2;
    gp(j) = (image_pressure(grid,{f},point+shift,1.7)- ...
        image_pressure(grid,{f},point-shift,1.7))/(2*step);
end
ePde = norm(-1.7*lap+gp);
mix = repmat(base,4,1);
families = ["S","R","T","D"];
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
actual = image_pressure(grid,lambda,x,1.7);
separate = zeros(size(actual));
for j = 1:4
    one.body.basis = mix(j);
    separate = separate+image_pressure(one,raw(j),x,1.7);
end
eMix = max(abs(actual-separate));
eMu = max(abs(image_pressure(grid,lambda,x,0.6)-0.6/1.7*actual));
fprintf(['image_pressure self-test: S %.3e, R %.3e, D %.3e, ', ...
    'T %.3e, Q_s %.3e, mixture %.3e, viscosity %.3e, PDE %.3e; ', ...
    'pressure scale %.3e\n'], ...
    eS,eR,eD,eT,eQ,eMix,eMu,ePde,max(abs(expected)));
if max([eS,eR,eD,eT,eMix,eMu]) > ...
        1e-12*max(1,max(abs(expected))) || ...
        eQ > 1e-5*max(1,max(abs(qExpected))) || ...
        ePde > 1e-6*max(1,norm(gp))
    error('image_pressure:SelfTestFailed','Pressure self-test failed.');
end
derivative = base;
derivative.family = "S_y";
derivative.points = [0.1 0 0];
derivative.normalDirections = [0 1 0];
derivative.derivativeDistances = 0.2;
grid.body = struct('basis',derivative);
actual = image_pressure(grid,{f},x,1.7);
expected = 1.7*stokes_directional_stokeslet_matrix( ...
    derivative.points,x,derivative.normalDirections, ...
    derivative.derivativeDistances,"pressure")*f;
assert(norm(actual-expected,inf)<1e-13);
end
