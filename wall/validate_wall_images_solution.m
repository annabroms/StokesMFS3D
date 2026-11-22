function details = validate_wall_images_solution(problem,lambda,U)
%VALIDATE_WALL_IMAGES_SOLUTION Report physical surface errors independently.
% lambda and U are normalized; source coefficients are body/family-grouped.
% All surface fields here use DIRECT complete-wall kernels, independently of
% GMRES's hybrid evaluation. The collocation and independent errors divide
% by the same largest rigid speed over the sampled surface sets. Zero rigid
% speed selects absolute error. No accuracy goal triggers an assertion.
% Independent targets include full-sphere quadrature, broad/inner Fibonacci
% caps at EVERY wall and pair contact and the exact closest surface points.
% Points in returned details are physical; velocities use velocity unit one.

timer = tic;
opt = problem.options; a = opt.radius;
uvec = apply_wall_images_velocity(lambda,problem,problem.rout,true);
u_bndry = zeros(size(uvec));
Fvec = zeros(6*problem.P,1);
for k = 1:problem.P
    rows = problem.targetOffsets(k)+(1:3*problem.M(k));
    cols = problem.sourceOffsets(k)+(1:3*problem.N_total(k));
    u_bndry(rows) = problem.body{k}.Kout*U(6*(k-1)+(1:6));
    Fvec(6*(k-1)+(1:6)) = problem.body{k}.Kin.'*lambda(cols);
end
speed = max(vecnorm(reshape(u_bndry,3,[]),2,1));
rcheck = zeros(0,3); ucheck = zeros(0,1); computed = zeros(0,1);
if opt.compute_residual
    sphere = setupsurfquad(ellipsoid_param(1,1,1),opt.validationResolution);
    cap = struct('collocationPatchMode',"fibonacci", ...
        'collocationPatchPointCount',opt.validationCount);
    for k = 1:problem.P
        points = problem.q(k,:)+sphere.x.';
        E = struct('axes',ones(3,1),'Q',eye(3),'c',problem.q(k,:).');
        % Check wall and neighbour directions even if their image lines are off.
        contacts = problem.q(k,:)-[0,0,1];
        pairs = problem.pairs;
        neighbours = [pairs(pairs(:,1)==k,2);pairs(pairs(:,2)==k,1)];
        for j = neighbours.'
            contacts(end+1,:) = problem.q(k,:)+ ...
                (problem.q(j,:)-problem.q(k,:))/norm(problem.q(j,:)-problem.q(k,:)); %#ok<AGROW>
        end
        for j = 1:size(contacts,1)
            for angle = [opt.capAngle,opt.capAngle/6]
                cap.collocationPatchAngle = angle;
                points = [points;collocation_surface_patch(E,contacts(j,:).',cap)]; %#ok<AGROW>
            end
        end
        points = uniquetol([points;contacts],1e-12,'ByRows',true,'DataScale',1);
        rcheck = [rcheck;points]; %#ok<AGROW>
        ucheck = [ucheck;getKmat(points,problem.q(k,:))*U(6*(k-1)+(1:6))]; %#ok<AGROW>
    end
    computed = apply_wall_images_velocity(lambda,problem,rcheck,true);
    speed = max(speed,max(vecnorm(reshape(ucheck,3,[]),2,1)));
end
absolute = vecnorm(reshape(uvec-u_bndry,3,[]).',2,2);
relative = absolute;
pointErrors = vecnorm(reshape(computed-ucheck,3,[]).',2,2);
if speed > 0
    relative = relative/speed; pointErrors = pointErrors/speed;
end
surfaceError = NaN;
if opt.compute_residual, surfaceError = max(pointErrors); end
% Sample the wall beneath and around each particle; wall points are never fitted.
wall = zeros(0,3);
for k = 1:problem.P
    [x,y] = meshgrid(problem.q(k,1)+[-1,0,1],problem.q(k,2)+[-1,0,1]);
    wall = [wall;x(:),y(:),zeros(numel(x),1)]; %#ok<AGROW>
end
wallVelocity = apply_wall_images_velocity(lambda,problem,wall,true);
forceUnits = repmat([opt.mu*a*ones(3,1);opt.mu*a^2*ones(3,1)],problem.P,1);
details = struct('problem',problem,'options',opt,'N',problem.N, ...
    'M',problem.M,'N_total',problem.N_total,'lambda',lambda, ...
    'sourceUnits',problem.sourceUnits,'rin',a*problem.rin+[0,0,opt.wall_height], ...
    'rout',a*problem.rout+[0,0,opt.wall_height], ...
    'collocationResidual',relative,'collocationAbsoluteResidual',absolute, ...
    'collocationError',max(relative),'normalizationSpeed',speed, ...
    'isRelative',speed>0,'validationPoints',a*rcheck+[0,0,opt.wall_height], ...
    'validationResidual',pointErrors,'surfaceError',surfaceError, ...
    'realizedForceTorque',Fvec.*forceUnits,'rigidVelocity', ...
    U./repmat([1;1;1;a;a;a],problem.P,1),'wallError',norm(wallVelocity,inf), ...
    'retainedRanks',problem.ranks,'timing',problem.timing,'figures',[]);
details.timing.validationSeconds = toc(timer);
end
