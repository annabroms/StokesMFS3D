function check = wall_force_torque_check(problem,lambda)
%WALL_FORCE_TORQUE_CHECK Compare physical source moments with full Blake traction.
% blake.fields returns sigma*n for outward particle normals; negate this
% traction to recover force/torque ON THE FLUID, matching K_in'*lambda.
% Error norms scale force and torque by problem.forceUnits before comparison.
if nargin == 0, self_test(); check = []; return; end
o = problem.options; a = o.radius(1); P = problem.P;
s = ellipsoid_param(a,a,a); s = setupsurfquad(s,o.tractionResolution);
integrated = zeros(6*P,1); moments = integrated;
for p = 1:P
    q = problem.physical.q(p,:); x = q+s.x.';
    opts = struct('method','blake','mu',o.mu,'wall_height',o.wall_height, ...
        'requested',{{'traction'}},'normals',s.nx);
    fields = blake.fields(problem.physical.rin.',reshape(lambda,3,[]),x.',opts);
    traction = -fields.traction;
    force = sum(traction.*s.w,2);
    torque = sum(cross(s.x,traction,1).*s.w,2);
    rb = (p-1)*6+(1:6); ns = (p-1)*3*problem.N+(1:3*problem.N);
    integrated(rb) = [force;torque];
    moments(rb) = getKmat(problem.physical.rin((p-1)*problem.N+(1:problem.N),:),q).'*lambda(ns);
end
scaledError = (integrated-moments)./problem.forceUnits;
scaledMoment = moments./problem.forceUnits;
check = struct('integrated',integrated,'sourceMoments',moments, ...
    'absoluteError',norm(scaledError,inf),'relativeError',NaN);
if norm(scaledMoment,inf)>0
    check.relativeError = norm(scaledError,inf)/norm(scaledMoment,inf);
end
end

function self_test()
% Use a single off-centre force inside an enclosing sphere near the wall.
q = [0,0,1.1]; o = wall_solver_options(...
    struct('N',18,'backend',"blake",'computeResidual',false)); o.radius = 1;
[s,t] = init_spheres(q,o.proxyRadius, ...
    o.N,o.collocationRatio,1);
p = prepare_wall_problem(q,s,t,o);
p.N = 1; p.physical.rin = q+[0.1,-0.08,0.15];
c = wall_force_torque_check(p,[0.7;-0.4;0.2]);
assert(c.relativeError<1e-10);
fprintf('wall_force_torque_check: PASS; complete Blake traction/moment relative error %.3e.\n', ...
    c.relativeError);
figure('Name','Blake traction versus source moments: self-test');
bar([c.sourceMoments,c.integrated]);
xticks(1:6); xticklabels({'Fx','Fy','Fz','Tx','Ty','Tz'});
ylabel('Force / torque (a=1, mu=1)'); grid on;
legend('Proxy-source moments','Integrated Blake traction','Location','best');
title(sprintf('Complete traction check: relative error %.2e',c.relativeError));
end
