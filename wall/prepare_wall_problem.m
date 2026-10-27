function problem = prepare_wall_problem(q,rin,rout,opt)
%PREPARE_WALL_PROBLEM Build reusable geometry, kernel, and one-body factors.
%   problem = PREPARE_WALL_PROBLEM(q,rin,rout,opt) prepares the data shared by
%   solve_resistance_wall and solve_mobility_wall for equal spheres above a
%   stationary no-slip plane. Call it once for a fixed geometry and pass the
%   result back through opt.preparedProblem when solving either or both modes.
%
%   Inputs:
%     q     P-by-3 sphere centres in physical coordinates.
%     rin   (N*P)-by-3 proxy-source positions, grouped by sphere.
%     rout  (M*P)-by-3 surface collocation positions, grouped by sphere.
%     opt   Wall solver options; see WALL_SOLVER_OPTIONS. The radius may be a
%           scalar or a length-P vector of equal radii. wall_height locates
%           the plane z=wall_height. backend selects the physical Blake
%           interaction evaluator; preconditioner is "wall" or "freeSpace".
%
%   Geometry requirements: all spheres must be strictly above the plane and
%   mutually nonoverlapping. Each sphere must have N>=3 sources and M>=N
%   targets. Sources must lie on the proxy sphere of radius proxyRadius*a;
%   targets must lie on the physical sphere of radius a. Components of
%   vectorized fields are interleaved by node: [x1;y1;z1;x2;...].
%
%   Returned fields:
%     options       Validated options, with preparedProblem cleared.
%     P,N,M         Sphere count and source/target counts per sphere.
%     physical      Original dimensional q, rin, and rout.
%     q,rin,rout    Coordinates translated to the wall plane and divided by
%                   a. These normalized arrays are used by the kernel.
%     kernel        Prepared global Blake kernel for the requested backend.
%     body{p}       Per-sphere H=K_in', K=K_out, C=H'/(H*H'), local offsets,
%                   and normalized height. In oneBody mode it also stores
%                   Gpre and truncated-SVD resistance/mobility factors.
%     velocityUnits, forceUnits, strengthUnit
%                   Dimensional conversion factors used by the solvers.
%     minWallGap, minBodyGap
%                   Minimum physical clearances; minBodyGap is NaN for P=1.
%     validation    Independent surface-check points and a prepared kernel
%                   when computeResidual=true; otherwise this field is absent.
%     dense         Global truncated-SVD factors only when solverMode is
%                   "denseSVD". This diagnostic mode is restricted to P<=2
%                   and explicitly assembles dense global matrices.
%     timing        Setup, assembly, factorization times and uniqueBlocks,
%                   the number of distinct one-body factorizations made.
%
%   The one-body resistance factor is based on Gpre. The mobility factor is
%   based on Gpre+(K-Gpre*C)*H. A "wall" preconditioner includes the wall in
%   Gpre and can be reused only between identical body grids at the same
%   height. A "freeSpace" preconditioner omits the wall only in Gpre and can
%   be reused across heights; the global physical kernel remains Blake in
%   both cases. Identical eligible blocks are factored once and shared.
%
%   In the default oneBody mode, no global dense interaction matrix is formed:
%   the global kernel handle is applied matrix-free during GMRES. Dense local
%   one-body blocks are formed for their truncated-SVD preconditioners.
%   Invalid geometry or mismatched proxy/target radii raises a wall:* error.
if nargin == 0, self_test(); problem = []; return; end
setupTimer = tic;
opt = wall_solver_options(opt);
opt.preparedProblem = [];
validateattributes(q,{'double'},{'real','finite','2d','ncols',3,'nonempty'});
P = size(q,1); a = opt.radius(1);
assert(isscalar(opt.radius) || numel(opt.radius) == P);
assert(size(rin,2) == 3 && size(rout,2) == 3);
N = size(rin,1)/P; M = size(rout,1)/P;
assert(N == floor(N) && M == floor(M) && N >= 3 && M >= N);
assert(all(isfinite(rin(:))) && all(isfinite(rout(:))));
assert(isreal(rin) && isreal(rout));
clearance = q(:,3)-a-opt.wall_height;
assert(all(clearance > 0),'wall:WallIntersection','Particles intersect the wall.');
nearest = inf(P,1);
for p = 1:P
    others = [1:p-1,p+1:P];
    if ~isempty(others)
        nearest(p) = min(vecnorm(q(others,:)-q(p,:),2,2)-2*a);
    end
end
assert(all(nearest > 0),'wall:Overlap','Particles must have positive separation.');
shift = [0,0,opt.wall_height];
problem = struct('options',opt,'P',P,'N',N,'M',M, ...
    'q',(q-shift)/a,'rin',(rin-shift)/a,'rout',(rout-shift)/a, ...
    'physical',struct('q',q,'rin',rin,'rout',rout), ...
    'velocityUnits',repmat([1;1;1;1/a;1/a;1/a],P,1), ...
    'forceUnits',repmat([opt.mu*a*ones(3,1);opt.mu*a^2*ones(3,1)],P,1), ...
    'strengthUnit',opt.mu*a,'minWallGap',min(clearance),'minBodyGap',NaN);
if P > 1, problem.minBodyGap = min(nearest); end
kernelOpt = struct('kernel',"blake",'backend',opt.backend,'mu',1, ...
    'wall_height',0,'fmmTolerance',opt.fmmTolerance);
problem.kernel = prepare_stokes_kernel(problem.rin,problem.rout,kernelOpt);
problem.body = cell(P,1);
assemblySeconds = 0; factorizationSeconds = 0; uniqueBlocks = 0;
for p = 1:P
    ns = (p-1)*N+(1:N); ms = (p-1)*M+(1:M);
    source = problem.rin(ns,:); target = problem.rout(ms,:);
    offsetSource = source-problem.q(p,:);
    offsetTarget = target-problem.q(p,:);
    assert(max(abs(vecnorm(offsetSource,2,2)-opt.proxyRadius)) < 1e-9, ...
        'wall:ProxyGeometry','Sources must be on the configured proxy sphere.');
    assert(max(abs(vecnorm(offsetTarget,2,2)-1)) < 1e-9, ...
        'wall:SurfaceGeometry','Targets must be on the particle sphere.');
    body = struct('H',getKmat(source,problem.q(p,:)).', ...
        'K',getKmat(target,problem.q(p,:)), ...
        'offsetSource',offsetSource,'offsetTarget',offsetTarget, ...
        'height',problem.q(p,3));
    body.C = body.H.'/(body.H*body.H.');
    reused = 0;
    for j = 1:p-1
        other = problem.body{j};
        sameHeight = opt.preconditioner == "freeSpace" || body.height==other.height;
        if sameHeight && norm(offsetSource-other.offsetSource,inf)<1e-12 && ...
                norm(offsetTarget-other.offsetTarget,inf)<1e-12
            reused = j; break
        end
    end
    if opt.solverMode == "oneBody"
        if reused
            fprintf('Reusing one-body preconditioner for body %d/%d (body %d).\n', ...
                p,P,reused);
            other = problem.body{reused};
            body.Gpre = other.Gpre;
            body.resistance = other.resistance;
            body.mobility = other.mobility;
        else
            fprintf('Computing one-body preconditioner for body %d/%d...\n',p,P);
            timer = tic;
            blockOpt = kernelOpt; blockOpt.backend = "blake";
            if opt.preconditioner == "freeSpace", blockOpt.kernel = "freeSpace"; end
            body.Gpre = stokes_kernel_matrix(source,target,blockOpt);
            A = body.Gpre+(body.K-body.Gpre*body.C)*body.H;
            assemblySeconds = assemblySeconds+toc(timer);
            timer = tic;
            [Y,U,s] = getPseudoFactors(body.Gpre,opt.svdCutoff,0);
            body.resistance = struct('Y',Y,'U',U,'singularValues',s,'rank',numel(s));
            [Y,U,s] = getPseudoFactors(A,opt.svdCutoff,0);
            body.mobility = struct('Y',Y,'U',U,'singularValues',s,'rank',numel(s));
            factorizationSeconds = factorizationSeconds+toc(timer);
            uniqueBlocks = uniqueBlocks+1;
        end
    end
    problem.body{p} = body;
end
if opt.solverMode == "denseSVD"
    assert(P <= 2,'wall:DenseSize','Dense reference mode is restricted to one or two bodies.');
    timer = tic;
    G = stokes_kernel_matrix(problem.rin,problem.rout,kernelOpt);
    A = G;
    for p = 1:P
        cols = (p-1)*3*N+(1:3*N); rows = (p-1)*3*M+(1:3*M);
        b = problem.body{p};
        A(:,cols) = G(:,cols)-(G(:,cols)*b.C)*b.H;
        A(rows,cols) = A(rows,cols)+b.K*b.H;
    end
    assemblySeconds = toc(timer);
    timer = tic;
    [Y,U,s] = getPseudoFactors(G,opt.svdCutoff,0);
    problem.dense.resistance = struct('Y',Y,'U',U,'singularValues',s,'rank',numel(s));
    [Y,U,s] = getPseudoFactors(A,opt.svdCutoff,0);
    problem.dense.mobility = struct('Y',Y,'U',U,'singularValues',s,'rank',numel(s));
    factorizationSeconds = toc(timer);
end
if opt.computeResidual
    sphere = ellipsoid_param(1,1,1);
    sphere = setupsurfquad(sphere,opt.validationResolution);
    checks = zeros(P*sphere.N,3);
    for p = 1:P
        checks((p-1)*sphere.N+(1:sphere.N),:) = sphere.x.'+problem.q(p,:);
    end
    problem.validation = struct('points',checks,'count',sphere.N, ...
        'kernel',prepare_stokes_kernel(problem.rin,checks,kernelOpt));
end
totalSeconds = toc(setupTimer);
problem.timing = struct('setupSeconds',totalSeconds-assemblySeconds-factorizationSeconds, ...
    'assemblySeconds',assemblySeconds,'factorizationSeconds',factorizationSeconds, ...
    'totalPreparationSeconds',totalSeconds,'uniqueBlocks',uniqueBlocks);
end

function self_test()
% Verify scaling, lateral factor reuse, and distinct-height wall blocks.
o = struct('N',30,'backend',"blake",'computeResidual',false,'radius',2,'mu',3);
q = [0,0,3;6,0,3]; [s,t,o] = init_spheres_wall(q,2,o);
p = prepare_wall_problem(q,s,t,o);
assert(p.timing.uniqueBlocks == 1 && p.strengthUnit == 6);
fprintf('prepare_wall_problem: equal heights reuse %d block; strength unit mu*a=%.3g.\n', ...
    p.timing.uniqueBlocks,p.strengthUnit);
q(2,3) = 4; [s,t,o] = init_spheres_wall(q,2,o);
p = prepare_wall_problem(q,s,t,o);
assert(p.timing.uniqueBlocks == 2);
fprintf('prepare_wall_problem: PASS; unequal heights require %d separate wall blocks.\n', ...
    p.timing.uniqueBlocks);
end
