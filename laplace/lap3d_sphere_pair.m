function pair = lap3d_sphere_pair(radius,delta)
%LAP3D_SPHERE_PAIR Shared equal-sphere geometry and reflection accumulation.
% Inputs are the positive radius and closest surface gap, validated by the
% experiment entry points.  This geometry is reused by the single run and sweep.

semiAxes = radius*ones(3,1);
pair.E1.axes = semiAxes;
pair.E1.Q = eye(3);
pair.E1.c = [-(radius+delta/2);0;0];
pair.E2.axes = semiAxes;
pair.E2.Q = eye(3);
pair.E2.c = [radius+delta/2;0;0];
pair.p1 = [-delta/2;0;0];
pair.p2 = [delta/2;0;0];
pair.gap = delta;
pair.mid = zeros(3,1);
pair.n = [1;0;0];
pair.T = [0,0;1,0;0,1];
pair.K1 = eye(2)/radius;
pair.K2 = eye(2)/radius;
pair.kappa1_principal = [1,1]/radius;
pair.kappa2_principal = [1,1]/radius;
pair.label = sprintf('Two equal spheres of radius %.4g',radius);

numberOfAngles = 721;
pair.theta = linspace(0,pi,numberOfAngles+1);
pair.theta(end) = [];
pair.kappa1 = ones(1,numberOfAngles)/radius;
pair.kappa2 = pair.kappa1;
pair.R1 = radius*ones(1,numberOfAngles);
pair.R2 = pair.R1;

fixedPointOffset = sqrt(delta*(4*radius+delta))/2;
pair.zminus = -fixedPointOffset*ones(1,numberOfAngles);
pair.zplus = fixedPointOffset*ones(1,numberOfAngles);
pair.Pminus = pair.mid+pair.n*pair.zminus;
pair.Pplus = pair.mid+pair.n*pair.zplus;
end
