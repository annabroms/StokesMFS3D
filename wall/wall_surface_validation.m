function validation = wall_surface_validation(problem,lambda,U)
%WALL_SURFACE_VALIDATION Check rigid motion on an independent surface grid.
% No extra wall-facing targets or patches are added. The error is maximum
% pointwise velocity discrepancy / global maximum rigid speed; zero data
% instead report an absolute error. Inputs and reported arrays are physical.
if nargin == 0, self_test(); validation = []; return; end
assert(problem.options.computeResidual,'Prepare with computeResidual=true.');
assert(numel(lambda)==3*problem.N*problem.P && numel(U)==6*problem.P);
v = problem.validation;
velocity = reshape(apply_stokes_kernel(lambda/problem.strengthUnit,v.kernel),3,[]).';
desired = zeros(size(velocity));
u = U./problem.velocityUnits;
for p = 1:problem.P
    rows = (p-1)*v.count+(1:v.count); rb = (p-1)*6+(1:6);
    desired(rows,:) = reshape(getKmat(v.points(rows,:),problem.q(p,:))*u(rb),3,[]).';
end
absolute = vecnorm(velocity-desired,2,2);
speed = max(vecnorm(desired,2,2));
relative = absolute;
if speed > 0, relative = absolute/speed; end
a = problem.options.radius(1);
points = a*v.points+[0,0,problem.options.wall_height];
validation = struct('error',max(relative),'absoluteError',max(absolute), ...
    'normalizationSpeed',speed,'isRelative',speed>0,'points',points, ...
    'velocity',velocity,'desiredVelocity',desired,'pointErrors',relative);
end

function self_test()
% Verify zero flow and that all check points lie on the physical sphere.
q = [0,0,3]; o = wall_solver_options(struct('N',18, ...
    'backend',"blake",'validationResolution',[12,8])); o.radius = 2;
[s,t] = init_spheres(q,o.proxyRadius, ...
    o.N,o.collocationRatio,2);
p = prepare_wall_problem(q,s,t,o);
v = wall_surface_validation(p,zeros(3*size(s,1),1),zeros(6,1));
assert(v.error==0 && ~v.isRelative);
assert(max(abs(vecnorm(v.points-q,2,2)-2))<1e-12);
fprintf('wall_surface_validation: PASS; %d independent surface points, zero-flow absolute error %.3e.\n', ...
    size(v.points,1),v.absoluteError);
fprintf('  Surface-radius error %.3e; no extra wall-facing sampling.\n', ...
    max(abs(vecnorm(v.points-q,2,2)-2)));
end
