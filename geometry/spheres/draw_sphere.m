figure; clf
set(gcf,'Color','w')
ax = axes;
hold(ax,'on')
axis(ax,'equal')
axis(ax,'off')
view(ax,35,20)

%% Sphere
R = 1;
[xs,ys,zs] = sphere(100);

surf(R*xs,R*ys,R*zs, ...
    'FaceColor',[0.85 0.8 0.85], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.8);

lighting gouraud
%camlight headlight
camlight left
material dull

%% Parameters
arrowLW = 2.2;
fs = 20;

%% Net force: straight arrow
F0 = [0, 0, 0];
F  = [1.45, 0.35, 0.15];

drawArrow3(F0,F, ...
    'Color',[0.85 0.15 0.10], ...
    'LineWidth',arrowLW, ...
    'HeadLength',0.18, ...
    'HeadRadius',0.07);

text(F(1)+0.28,F(2),F(3), ...
    '$\mathbf{f}$', ...
    'Interpreter','latex', ...
    'FontSize',fs, ...
    'Color',[0.85 0.15 0.10]);

%% Translational velocity: straight arrow, displaced below sphere
U0 = [0,0,0];
U  = [ 1.25, -1.55, -0.75];

drawArrow3(U0,U, ...
    'Color',[0.10 0.35 0.85], ...
    'LineWidth',arrowLW, ...
    'HeadLength',0.18, ...
    'HeadRadius',0.07);

text(U(1)+0.18,U(2),U(3), ...
    '$\mathbf{v}$', ...
    'Interpreter','latex', ...
    'FontSize',fs, ...
    'Color',[0.10 0.35 0.85]);

%%  normal direction
U0 = [0,0,0];
U  = -[ 1.2 1.2 1.2];%/sqrt(3);

drawArrow3(U0,U, ...
    'Color','k', ...
    'LineWidth',arrowLW, ...
    'HeadLength',0.18, ...
    'HeadRadius',0.07);

text(U(1)+0.18,U(2)+0.2,U(3), ...
    '$\mathbf{n}$', ...
    'Interpreter','latex', ...
    'FontSize',fs, ...
    'Color','k');

%% Net torque: curved arrow around sphere
drawCircularArrow3( ...
    [0,0,0], ...          % centre
    [0,0,1], ...          % axis of rotation
    1.25, ...             % radius
    deg2rad(20), ...      % start angle
    deg2rad(150), ...     % end angle
    'Color',[0.15 0.55 0.20], ...
    'LineWidth',arrowLW, ...
    'HeadLength',0.16, ...
    'HeadRadius',0.06);

text(-1.10,0.95,0.10, ...
    '$\mathbf{t}$', ...
    'Interpreter','latex', ...
    'FontSize',fs, ...
    'Color',[0.15 0.55 0.20]);

%% Angular velocity: visible curved arrow outside sphere
drawCircularArrow3( ...
    [0,0,0], ...              % centre
    [1,1,1], ...              % rotation axis
    1.35, ...                 % radius: larger than sphere radius
    deg2rad(-60), ...
    deg2rad(20), ...
    'Color',[0.55 0.15 0.75], ...
    'LineWidth',arrowLW, ...
    'HeadLength',0.18, ...
    'HeadRadius',0.07);

text(0.95,-0.55,0.95, ...
    '$\mathbf{\omega}$', ...
    'Interpreter','latex', ...
    'FontSize',fs, ...
    'Color',[0.55 0.15 0.75]);

%% Optional: coordinate axes limits
xlim([-1.7 1.8])
ylim([-1.7 1.5])
zlim([-1.2 1.4])

function drawArrow3(p0,p1,varargin)
% Draw a 3D arrow from p0 to p1 using a line and a cone arrowhead.

opts.Color      = [0 0 0];
opts.LineWidth  = 2;
opts.HeadLength = 0.15;
opts.HeadRadius = 0.05;

opts = parseOpts(opts,varargin{:});

p0 = p0(:).';
p1 = p1(:).';
v = p1 - p0;
L = norm(v);

if L == 0
    return
end

e = v/L;

% Shaft ends before the cone
pShaftEnd = p1 - opts.HeadLength*e;

plot3([p0(1),pShaftEnd(1)], ...
      [p0(2),pShaftEnd(2)], ...
      [p0(3),pShaftEnd(3)], ...
      'Color',opts.Color, ...
      'LineWidth',opts.LineWidth);

% Cone in local z-direction
[nx,ny,nz] = cylinder([opts.HeadRadius 0],32);
nz = opts.HeadLength*nz;

% Move cone so local tip is at z = HeadLength
pts = [nx(:), ny(:), nz(:)];

% Rotate local z-axis to direction e
R = rotationFromZ(e);
pts = pts*R.';

% Translate so cone base is at pShaftEnd and tip at p1
pts = pts + pShaftEnd;

nx = reshape(pts(:,1),size(nx));
ny = reshape(pts(:,2),size(ny));
nz = reshape(pts(:,3),size(nz));

surf(nx,ny,nz, ...
    'FaceColor',opts.Color, ...
    'EdgeColor','none');
end

function drawCircularArrow3(c,axisVec,r,theta0,theta1,varargin)
% Draw a circular 3D arrow around axisVec, centred at c.

opts.Color      = [0 0 0];
opts.LineWidth  = 2;
opts.HeadLength = 0.15;
opts.HeadRadius = 0.05;

opts = parseOpts(opts,varargin{:});

c = c(:).';
a = axisVec(:).';
a = a/norm(a);

% Build two orthonormal vectors perpendicular to a
tmp = [1 0 0];
if abs(dot(tmp,a)) > 0.9
    tmp = [0 1 0];
end

e1 = tmp - dot(tmp,a)*a;
e1 = e1/norm(e1);
e2 = cross(a,e1);

theta = linspace(theta0,theta1,120);

pts = c + r*cos(theta(:))*e1 + r*sin(theta(:))*e2;

plot3(pts(:,1),pts(:,2),pts(:,3), ...
    'Color',opts.Color, ...
    'LineWidth',opts.LineWidth);

% Arrowhead tangent direction at end
pEnd = pts(end,:);
pPrev = pts(end-4,:);
drawArrow3(pPrev,pEnd, ...
    'Color',opts.Color, ...
    'LineWidth',opts.LineWidth, ...
    'HeadLength',opts.HeadLength, ...
    'HeadRadius',opts.HeadRadius);
end

function R = rotationFromZ(e)
% Rotation matrix mapping [0 0 1] to unit vector e.

e = e(:)/norm(e);
z = [0;0;1];

if norm(cross(z,e)) < 1e-14
    if dot(z,e) > 0
        R = eye(3);
    else
        R = diag([1 -1 -1]);
    end
    return
end

v = cross(z,e);
s = norm(v);
c = dot(z,e);

vx = [  0   -v(3)  v(2);
       v(3)   0   -v(1);
      -v(2) v(1)   0  ];

R = eye(3) + vx + vx^2*((1-c)/s^2);
end

function opts = parseOpts(opts,varargin)
% Small name-value parser.

for k = 1:2:numel(varargin)
    name = varargin{k};
    val  = varargin{k+1};

    if isfield(opts,name)
        opts.(name) = val;
    else
        error('Unknown option "%s".',name)
    end
end
end