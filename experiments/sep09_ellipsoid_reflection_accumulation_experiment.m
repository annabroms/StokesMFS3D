function results = sep09_ellipsoid_reflection_accumulation_experiment(caseName, delta, visualise, semiAxes)
%ELLIPSOID_REFLECTION_ACCUMULATION_EXPERIMENT
% Directional osculating-sphere model for two nearby ellipsoids.
%
%   results = sep09_ellipsoid_reflection_accumulation_experiment
%   results = sep09_ellipsoid_reflection_accumulation_experiment("all", 0.05)
%   results = sep09_ellipsoid_reflection_accumulation_experiment("tip", 0.02)
%   results = sep09_ellipsoid_reflection_accumulation_experiment("side", 0.02)
%   results = sep09_ellipsoid_reflection_accumulation_experiment("oblique", 0.02)
%   results = sep09_ellipsoid_reflection_accumulation_experiment("random", 0.02)
%   randomizes both rotations and the gap normal, with an exact support gap.
%   results = sep09_ellipsoid_reflection_accumulation_experiment(...,false)
%   suppresses figures while retaining the printed numerical summary.
%   results = sep09_ellipsoid_reflection_accumulation_experiment( ...
%       "tip",0.05,false,[1,0.7,0.7])
%   overrides the editable semiaxes without changing standalone defaults.
%
% The experiment:
%   1. Finds the closest points p1, p2 of the two ellipsoids.
%   2. Builds a common tangent basis at the gap.
%   3. Computes each body's normal curvature kappa_j(theta).
%   4. Replaces each normal section by its osculating circle/sphere,
%      R_j(theta) = 1/kappa_j(theta).
%   5. Computes the two fixed points of the corresponding two-sphere
%      inversion problem for every theta in [0,pi).
%
% The resulting accumulation points all lie on the common normal line.
% Thus the theta sweep predicts two axial intervals, not two off-axis curves.
%
% Three distinct positive semiaxes are supported by the curvature model.
% The first local axis sets the tip geometry and the theta reference, not
% a required symmetry axis.  These are osculating-sphere predictions,
% not an exact ellipsoid reflection or a proof of its singularity locations.
%
% No Optimization Toolbox is required.

if nargin < 1 || isempty(caseName)
    caseName = "all";
end
if nargin < 2 || isempty(delta)
    delta = 0.05;
end
if nargin < 3 || isempty(visualise)
    visualise = true;
end
if ~isscalar(visualise)
    error('visualise must be a scalar logical value.');
end
visualise = logical(visualise);

caseName = lower(string(caseName));

% Standalone geometry defaults are preserved.  A caller such as the
% Laplace experiment can supply its own geometry through the fourth input.
if nargin < 4 || isempty(semiAxes)
    a = 1;                            % major semiaxis
    b = 1; %0.70;                         % minor semiaxis
    semiAxes = [a; b; b];
end
validateattributes(semiAxes,{'numeric'}, ...
    {'real','finite','vector','numel',3,'positive'},mfilename,'semiAxes');
semiAxes = semiAxes(:);


switch caseName
    case "all"
        names = ["tip", "side", "oblique"];
    case {"tip", "side", "oblique", "random"}
        names = caseName;
    otherwise
        error('caseName must be "all", "tip", "side", "oblique", or "random".');
end

for k = 1:numel(names)
    [E1,E2,label,p1,p2] = make_case(names(k), semiAxes, delta);

    % For the oblique configuration, enforce the requested closest gap.
    if names(k) == "oblique"
        E2 = adjust_gap(E1,E2,delta);
    end

    if names(k) == "random"
        % The facing support points attain the separating-plane distance.
        % Hence their distance is the global closest gap, without a search.
        gap = norm(p2-p1);
    else
        [p1,p2,gap] = closest_points(E1,E2);
    end

    % Recenter for cleaner plots without changing the geometry.
    mid = 0.5*(p1+p2);
    E1.c = E1.c - mid;
    E2.c = E2.c - mid;
    p1 = p1 - mid;
    p2 = p2 - mid;

    out = analyse_pair(E1,E2,p1,p2,gap,label);
    if k == 1
        results = repmat(out,1,numel(names));
    else
        results(k) = out;
    end

    if visualise
        plot_pair(out, k);
        plot_theta_sweep(out, k);
    end

    fprintf('\n%s\n', repmat('=',1,68));
    fprintf('%s\n', label);
    fprintf('%s\n', repmat('=',1,68));
    fprintf('Closest gap = %.12g\n', out.gap);
    fprintf('Common normal n = [% .6f % .6f % .6f]^T\n', out.n);
    fprintf('\nCurvature matrix K1 in common tangent coordinates:\n');
    disp(out.K1)
    fprintf('Curvature matrix K2 in common tangent coordinates:\n');
    disp(out.K2)

    fprintf('Principal curvatures body 1: %.8g, %.8g\n', out.kappa1_principal);
    fprintf('Principal curvatures body 2: %.8g, %.8g\n', out.kappa2_principal);
    fprintf('Predicted z_- range: [% .8g, % .8g]\n', min(out.zminus), max(out.zminus));
    fprintf('Predicted z_+ range: [% .8g, % .8g]\n', min(out.zplus),  max(out.zplus));

    if max(out.zplus)-min(out.zplus) < 1e-8
        fprintf('The theta sweep is essentially degenerate: one point per body.\n');
    else
        fprintf('The theta sweep gives a nontrivial axial interval in each body.\n');
    end
end

end


function [E1,E2,label,p1,p2] = make_case(name, semiAxes, delta)
%MAKE_CASE supplies the fixed cases and exact-gap random ellipsoid pairs.

a = semiAxes(1);
b = semiAxes(2);

E1.axes = semiAxes;
E2.axes = semiAxes;
E1.Q = eye(3);
E2.Q = eye(3);
p1 = [];
p2 = [];

% Reuse the Z-Y-Z Euler rotation convention from geometry/ellipsoids.
Rz = @(t) [cos(t),-sin(t),0; sin(t),cos(t),0; 0,0,1];
Ry = @(t) [cos(t),0,-sin(t); 0,1,0; sin(t),0,cos(t)];

switch name
    case "tip"
        % Major axes and line of centers both along x.
        E1.c = [-(a + delta/2); 0; 0];
        E2.c = [ +(a + delta/2); 0; 0];
        label = 'Tip-to-tip prolate spheroids';

    case "side"
        % Major axes along x; line of centers along y.
        E1.c = [0; -(b + delta/2); 0];
        E2.c = [0; +(b + delta/2); 0];
        label = 'Side-by-side prolate spheroids';

    case "oblique"
        % Two differently oriented spheroids.  Start comfortably separated;
        % adjust_gap subsequently translates body 2 until the closest gap is delta.
        degrees = pi/180;
        E1.Q = Rz(18*degrees) * Ry(14*degrees);
        E2.Q = Rz(94.178379866655348*degrees) * ...
               Ry(32.598849049231667*degrees) * ...
               Rz(-122.57709199192384*degrees);

        E1.c = [0;0;0];
        d = [1.0; 0.62; 0.37];
        d = d/norm(d);
        E2.c = 3.25*d;
        label = 'Oblique, non-symmetric prolate spheroids';

    case "random"
        normal = randn(3,1);
        normal = normal/norm(normal);
        contacts = [-delta/2*normal,delta/2*normal];
        bodies = {E1,E2};
        for body = 1:2
            % Uniform pole direction and uniform twist, as in ellipsoid_cluster.
            pole = randn(3,1);
            pole = pole/norm(pole);
            alpha = atan2(pole(2),pole(1));
            beta = acos(pole(3));
            gamma = 2*pi*rand;
            rotation = Rz(alpha)*Ry(beta)*Rz(gamma);
            direction = (3-2*body)*normal;
            localDirection = rotation.'*direction;
            supportOffset = rotation*(semiAxes.^2.*localDirection) / ...
                norm(semiAxes.*localDirection);
            bodies{body}.Q = rotation;
            bodies{body}.c = contacts(:,body)-supportOffset;
        end
        E1 = bodies{1};
        E2 = bodies{2};
        p1 = contacts(:,1);
        p2 = contacts(:,2);
        label = 'Randomly positioned and oriented ellipsoids';
end
end


function out = analyse_pair(E1,E2,p1,p2,gap,label)
%ANALYSE_PAIR computes directional curvatures and surrogate fixed points.

n = (p2-p1)/gap;

% Choose theta=0 using the projection of body 1's first local axis if possible.
e1 = E1.Q(:,1) - n*(n.'*E1.Q(:,1));
if norm(e1) < 1e-10
    trial = [0;0;1];
    if abs(trial.'*n) > 0.9
        trial = [0;1;0];
    end
    e1 = trial - n*(n.'*trial);
end
e1 = e1/norm(e1);
e2 = cross(n,e1);
e2 = e2/norm(e2);
T = [e1,e2];

K1 = curvature_matrix(E1,p1,T);
K2 = curvature_matrix(E2,p2,T);

% theta and theta+pi give the same unoriented section plane and the same
% normal curvature, so [0,pi) contains each section exactly once.
Ntheta = 721;
theta = linspace(0,pi,Ntheta+1);
theta(end) = [];

v = [cos(theta); sin(theta)];
kappa1 = sum(v.*(K1*v),1);
kappa2 = sum(v.*(K2*v),1);

R1 = 1./kappa1;
R2 = 1./kappa2;

% Fixed points of alternating Kelvin inversions in two spheres of radii
% R1(theta), R2(theta), separated by the physical gap.
disc = gap .* (2*R1 + gap) .* (2*R2 + gap) .* ...
              (2*R1 + 2*R2 + gap);

zminus = (-gap.*(R1-R2) - sqrt(disc)) ./ ...
         (2*(R1+R2+gap));
zplus  = (-gap.*(R1-R2) + sqrt(disc)) ./ ...
         (2*(R1+R2+gap));

mid = 0.5*(p1+p2);
Pminus = mid*ones(1,Ntheta) + n*zminus;
Pplus  = mid*ones(1,Ntheta) + n*zplus;

lam1 = sort(eig(K1)).';
lam2 = sort(eig(K2)).';

out.label = label;
out.E1 = E1;
out.E2 = E2;
out.p1 = p1;
out.p2 = p2;
out.gap = gap;
out.mid = mid;
out.n = n;
out.T = T;
out.K1 = K1;
out.K2 = K2;
out.kappa1_principal = lam1;
out.kappa2_principal = lam2;
out.theta = theta;
out.kappa1 = kappa1;
out.kappa2 = kappa2;
out.R1 = R1;
out.R2 = R2;
out.zminus = zminus;
out.zplus = zplus;
out.Pminus = Pminus;
out.Pplus = Pplus;
end


function K = curvature_matrix(E,p,T)
%CURVATURE_MATRIX evaluates the second fundamental form in tangent basis T.

A = E.Q * diag(1./(E.axes.^2)) * E.Q.';
r = p-E.c;

gradF = 2*A*r;
H = 2*A;

K = T.'*H*T / norm(gradF);
K = 0.5*(K+K.');
end


function E2 = adjust_gap(E1,E2,targetGap)
%ADJUST_GAP translates body 2 along the current closest-point normal.

for it = 1:30
    [p1,p2,gap] = closest_points(E1,E2);

    if abs(gap-targetGap) < 2e-9
        return
    end

    n = (p2-p1)/gap;
    step = targetGap-gap;

    % When moving the bodies closer, avoid an overly aggressive translation.
    if step < 0
        step = max(step,-0.70*gap);
    end

    E2.c = E2.c + step*n;
end

[p1,p2,gap] = closest_points(E1,E2); %#ok<ASGLU>
warning('Gap adjustment stopped at %.4e instead of %.4e.', gap, targetGap);
end


function [p1,p2,gap] = closest_points(E1,E2)
%CLOSEST_POINTS finds the closest surface points using toolbox-free fminsearch.

d = E2.c-E1.c;
d = d/norm(d);

u1 = support_parameter(E1, d);
u2 = support_parameter(E2,-d);

q0 = [unit_to_angles(u1); unit_to_angles(u2)];

% A few deterministic multistarts make the search robust for the oblique case.
h = 0.45;
offsets = [zeros(4,1), ...
           h*eye(4), -h*eye(4), ...
           [h;h;0;0], [-h;h;0;0], ...
           [0;0;h;h], [0;0;-h;h]];

opts = optimset('Display','off', ...
                'TolX',1e-11, ...
                'TolFun',1e-15, ...
                'MaxIter',2500, ...
                'MaxFunEvals',10000);

best = inf;
qbest = q0;

for k = 1:size(offsets,2)
    qstart = q0 + offsets(:,k);
    [q,fval] = fminsearch(@(qq) distance_sq(qq,E1,E2), qstart, opts);
    if fval < best
        best = fval;
        qbest = q;
    end
end

p1 = surface_point(E1,qbest(1:2));
p2 = surface_point(E2,qbest(3:4));
gap = norm(p2-p1);
end


function val = distance_sq(q,E1,E2)
%DISTANCE_SQ is the squared distance between two parameterized surface points.

x1 = surface_point(E1,q(1:2));
x2 = surface_point(E2,q(3:4));
r = x2-x1;
val = r.'*r;
end


function x = surface_point(E,q)
%SURFACE_POINT maps azimuth/elevation parameters to an ellipsoid point.

az = q(1);
el = q(2);
u = [cos(el)*cos(az); cos(el)*sin(az); sin(el)];
x = E.c + E.Q*(E.axes.*u);
end


function u = support_parameter(E,d)
%SUPPORT_PARAMETER gives the unit-sphere parameter of the support point.

w = E.axes .* (E.Q.'*d);
u = w/norm(w);
end


function q = unit_to_angles(u)
%UNIT_TO_ANGLES converts a unit vector to azimuth/elevation.

u = u/norm(u);
q = [atan2(u(2),u(1)); atan2(u(3),hypot(u(1),u(2)))];
end


function plot_pair(out, figIndex)
%PLOT_PAIR shows bodies, closest points, normal, and predicted accumulation set.

figure(10+figIndex); clf
hold on

plot_ellipsoid(out.E1,48);
plot_ellipsoid(out.E2,48);

plot3(out.p1(1),out.p1(2),out.p1(3),'ko','MarkerFaceColor','k','MarkerSize',5);
plot3(out.p2(1),out.p2(2),out.p2(3),'ko','MarkerFaceColor','k','MarkerSize',5);
plot3([out.p1(1),out.p2(1)], ...
      [out.p1(2),out.p2(2)], ...
      [out.p1(3),out.p2(3)], 'k-', 'LineWidth',1.2);

% The theta-sweep points are collinear; the dense markers reveal the interval.
scatter3(out.Pminus(1,:),out.Pminus(2,:),out.Pminus(3,:),12,out.theta,'filled');
scatter3(out.Pplus(1,:), out.Pplus(2,:), out.Pplus(3,:), 12,out.theta,'filled');

% Draw theta=0 and theta=pi/2 directions in the common tangent plane.
L = 0.8*min(out.E1.axes);
e1 = out.T(:,1);
e2 = out.T(:,2);
m = out.mid;
plot3(m(1)+[-L,L]*e1(1),m(2)+[-L,L]*e1(2),m(3)+[-L,L]*e1(3), ...
      'k--','LineWidth',0.9);
plot3(m(1)+[-L,L]*e2(1),m(2)+[-L,L]*e2(2),m(3)+[-L,L]*e2(3), ...
      'k:','LineWidth',0.9);

axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
title({out.label, ...
       sprintf('gap = %.3g; accumulation points from \\theta \\in [0,\\pi)',out.gap)})
cb = colorbar;
cb.Label.String = '\theta';
end


function plot_theta_sweep(out, figIndex)
%PLOT_THETA_SWEEP shows radii and accumulation depth versus theta.

figure(20+figIndex); clf
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(out.theta*180/pi,out.R1,'LineWidth',1.4); hold on
plot(out.theta*180/pi,out.R2,'LineWidth',1.4);
grid on
xlim([0,180])
ylabel('R_j(\theta)')
legend('R_1','R_2','Location','best')
title(out.label)

nexttile
plot(out.theta*180/pi,out.zminus,'LineWidth',1.4); hold on
plot(out.theta*180/pi,out.zplus,'LineWidth',1.4);
grid on
xlim([0,180])
xlabel('\theta (degrees)')
ylabel('z_\pm(\theta)')
legend('z_-','z_+','Location','best')

title(tl,'Directional osculating-sphere prediction')
end


function plot_ellipsoid(E,N)
%PLOT_ELLIPSOID draws an oriented ellipsoid surface.

[x,y,z] = sphere(N);
P = [E.axes(1)*x(:), E.axes(2)*y(:), E.axes(3)*z(:)].';
P = E.Q*P + E.c;

X = reshape(P(1,:),size(x));
Y = reshape(P(2,:),size(y));
Z = reshape(P(3,:),size(z));

surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.18);
end
