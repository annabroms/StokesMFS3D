function problem = prepare_wall_images_problem(q,rin,rout,opt,mode)
%PREPARE_WALL_IMAGES_PROBLEM Build wall/pair grids and one-body SVD factors.
%   problem = prepare_wall_images_problem(q,rin,rout,opt,mode)
%   mode is 'resistance' or 'mobility'. No global dense matrix is built.
%
% q is P-by-3; rin and rout are particle-grouped ordinary proxy/surface
% grids with the same initial N and M on each equal sphere. Coordinates
% are physical. Internally lengths are divided by opt.radius and the wall
% is moved to z=0. Point vector components are interleaved [x1;y1;z1;...].
%
% Wall lines default to S,T_t,D; close-pair lines default to S,D,R.
% opt.wallImageFamilies and opt.pairImageFamilies select them independently.
% Tokens: S, R, D, T_t (along the line), T_n/T_n2 (transverse), or T (all
% three stresslet directions). "none" disables one type of line. Proxy S
% is always present. All families use
% complete no-slip half-space kernels, even on particle-pair lines.
% N is the proxy count. N_image counts geometric points on ONE image line.
% M(k) counts final surface targets; N_total(k) counts source vectors,
% counting each family separately at a shared line point.
%
% Flat options (missing fields are completed here, once):
%   radius=1, viscosity=1, wall_height=0  physical length/viscosity controls
%   SVDsolve=false                    true: global dense truncated SVD,
%                                     false: one-body Helsing GMRES
%   N_image=30, nearThreshold=0.15     line count; pair surface gap/radius
%   wallImageFamilies="ST_tD"         wall line families
%   pairImageFamilies="SDR"           particle-pair line families
%   imageOffsetFraction=0.05          line starts outside the proxy sphere
%   capAngle=pi/5, capCount=901        exponential fitting cap per contact
%   capLayers=30, capExponent=4        cap grading (includes closest pole)
%   tol=1e-12                         relative one-body SVD cutoff
%   maxit=400, gmres_tol=1e-9, gmres_verbose=1
%   fmm=1, fmm_tol=1e-12              FMM for S; T_t,D,R remain direct
%   compute_residual=1, plot=0         independent validation and figures
%   validationResolution=[32,40]       independent full-sphere quadrature
%   validationCount=2001              Fibonacci points per contact cap
%   batchSize=256                     direct-kernel target block size
% N_image=0 disables wall AND pair lines and their additional fitting caps.
% Unknown fields from init_spheres are retained, not interpreted as kernels.
%
% body{k}.Kin' maps coefficients to force/torque on the fluid. S contributes
% force and moment-arm torque, R contributes torque, T_t/D neither.
% body{k}.Kout maps [translation;rotation] to surface velocity.
% Y{k}*UU{k} applies the truncated one-body pseudoinverse. UU has retained
% left singular vectors in ROWS; Y contains right vectors divided by their
% singular values. U is reserved for rigid velocities in the public solvers.

if nargin < 4, opt = struct(); end
% Keep the old mu spelling usable; viscosity is the public option name.
if ~isfield(opt,'viscosity') && isfield(opt,'mu'), opt.viscosity = opt.mu; end
defaults = struct('radius',1,'viscosity',1,'wall_height',0,'N_image',30, ...
    'SVDsolve',false, ...
    'wallImageFamilies',"ST_tD",'pairImageFamilies',"SDR", ...
    'nearThreshold',0.15,'imageOffsetFraction',0.05, ...
    'capAngle',pi/5,'capCount',901,'capLayers',30,'capExponent',4, ...
    'tol',1e-12,'maxit',400,'gmres_tol',1e-9,'gmres_verbose',1, ...
    'fmm',1,'fmm_tol',1e-12,'compute_residual',1,'plot',0, ...
    'validationResolution',[32,40],'validationCount',2001,'batchSize',256);
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opt,names{k}), opt.(names{k}) = defaults.(names{k}); end
end
% Decode both combinations once. Masks use [S,R,T_n,T_n2,T_t,D].
families = false(2,6);
names = {'wallImageFamilies','pairImageFamilies'};
for k = 1:2
    combination = strtrim(string(opt.(names{k})));
    if ~isscalar(combination) || ismissing(combination)
        error('wall:ImageCombination','%s must be one combination string.',names{k});
    end
    if strcmpi(combination,"none")
        combination = "none";
    else
        families(k,:) = wall_image_combination(combination);
    end
    opt.(names{k}) = combination;
end
validateattributes(opt,{'struct'},{'scalar'});
validateattributes(q,{'double'},{'real','finite','nonempty','ncols',3});
validateattributes(rin,{'double'},{'real','finite','nonempty','ncols',3});
validateattributes(rout,{'double'},{'real','finite','nonempty','ncols',3});
validateattributes(opt.radius,{'double'},{'scalar','positive','finite'});
validateattributes(opt.viscosity,{'double'},{'scalar','real','positive','finite'});
opt.mu = opt.viscosity;
validateattributes(opt.wall_height,{'double'},{'scalar','real','finite'});
validateattributes(opt.N_image,{'double'},{'scalar','integer','nonnegative'});
validateattributes(opt.nearThreshold,{'double'},{'scalar','positive','finite'});
validateattributes(opt.imageOffsetFraction,{'double'},{'scalar','nonnegative','finite'});
validateattributes(opt.capAngle,{'double'},{'scalar','>',0,'<',pi});
validateattributes(opt.capLayers,{'double'},{'scalar','integer','positive'});
validateattributes(opt.capCount,{'double'},{'scalar','integer','>=',1+8*opt.capLayers});
validateattributes(opt.capExponent,{'double'},{'scalar','positive','finite'});
validateattributes(opt.tol,{'double'},{'scalar','>',0,'<',1});
validateattributes(opt.gmres_tol,{'double'},{'scalar','>',0,'<',1});
validateattributes(opt.fmm_tol,{'double'},{'scalar','>',0,'<',1});
validateattributes(opt.maxit,{'double'},{'scalar','integer','positive'});
validateattributes(opt.validationCount,{'double'},{'scalar','integer','>=',4});
validateattributes(opt.validationResolution,{'double'},{'numel',2,'integer','>=',4});
validateattributes(opt.batchSize,{'double'},{'scalar','integer','positive'});
for name = {'fmm','compute_residual','plot','gmres_verbose','SVDsolve'}
    validateattributes(opt.(name{1}),{'numeric','logical'},{'scalar','binary'});
end
opt.SVDsolve = logical(opt.SVDsolve);
mode = validatestring(mode,{'resistance','mobility'});
timer = tic;
P = size(q,1); N = size(rin,1)/P; M = size(rout,1)/P;
if N ~= floor(N) || M ~= floor(M) || N < 3 || M < N
    error('wall:ImagesGrid','Initial rin/rout must have equal N/M per particle, M>=N>=3.');
end
% Normalize once; the physical source conversion is returned in sourceUnits.
a = opt.radius;
q = (q-[0,0,opt.wall_height])/a;
rin = (rin-[0,0,opt.wall_height])/a;
rout = (rout-[0,0,opt.wall_height])/a;
if any(q(:,3) <= 1)
    error('wall:WallIntersection','Every sphere must have positive wall clearance.');
end
pairs = zeros(0,2); minGap = inf;
for k = 1:P
    for j = k+1:P
        gap = norm(q(j,:)-q(k,:))-2;
        if gap <= 0, error('wall:Overlap','Particles must have positive separation.'); end
        minGap = min(minGap,gap);
        if gap < opt.nearThreshold, pairs(end+1,:) = [k,j]; end %#ok<AGROW>
    end
end
problem = struct('options',opt,'mode',mode,'P',P,'N',N,'q',q, ...
    'pairs',pairs,'minWallGap',a*min(q(:,3)-1),'minBodyGap',a*minGap);
if P == 1, problem.minBodyGap = NaN; end
problem.M = zeros(P,1); problem.N_total = zeros(P,1);
problem.body = cell(P,1); problem.Y = cell(P,1); problem.UU = cell(P,1);
problem.rin = zeros(0,3); problem.rout = zeros(0,3);
problem.sourceUnits = zeros(0,1); problem.powers = zeros(0,1);
problem.sourceOffsets = zeros(P+1,1); problem.targetOffsets = zeros(P+1,1);
problem.ranks = zeros(P,1); problem.singularValues = cell(P,1);
% This template matches the established wall matrix and coefficient adapters.
block = struct('family',"S",'points',zeros(0,3), ...
    'normalDirections',zeros(0,3),'scale',1,'label',"proxy S");
cap = struct('collocationPatchMode',"exponential", ...
    'collocationPatchAngle',opt.capAngle,'collocationPatchPointCount',opt.capCount, ...
    'collocationPatchLayers',opt.capLayers,'collocationExponent',opt.capExponent, ...
    'collocationMinAzimuth',8,'collocationAzimuthOversampling',1);
kernelOpt = struct('method','blake','mu',1,'wall_height',0);
assemblySeconds = 0; factorSeconds = 0;
tokens = ["S","R","T_n","T_n2","T_t","D"];
for k = 1:P
    rin_block = rin((k-1)*N+(1:N),:);
    rout_block = rout((k-1)*M+(1:M),:);
    Rp = mean(vecnorm(rin_block-q(k,:),2,2));
    if Rp >= 1 || Rp <= 0 || ...
            max(abs(vecnorm(rin_block-q(k,:),2,2)-Rp)) > 1e-9 || ...
            max(abs(vecnorm(rout_block-q(k,:),2,2)-1)) > 1e-9
        error('wall:ImagesGrid','rin must be a proxy sphere; rout must lie on radius a.');
    end
    sphere = struct('axes',ones(3,1),'Q',eye(3),'c',q(k,:).');
    basis = block; basis.points = rin_block;
    contacts = zeros(0,3);
    % Reflect the centre: the equal reflected-pair gap is twice the wall gap.
    h = q(k,3); gap = h-1;
    rstar = 1/(h+sqrt(gap*(2+gap)));
    if opt.N_image > 0 && any(families(1,:)) && ...
            Rp*(1+opt.imageOffsetFraction) < rstar
        mirror = q(k,:); mirror(3) = -mirror(3);
        line = image_line_points(q(k,:),mirror,1,Rp,opt.N_image, ...
            opt.imageOffsetFraction,1);
        % This order preserves the default [wall S;wall T_t;wall D] layout.
        normals = [zeros(2,3);0,-1,0;-1,0,0;0,0,-1;zeros(1,3)];
        for f = [1,5,6,2,3,4]
            if ~families(1,f), continue; end
            b = block; b.family = tokens(f); b.points = line;
            b.label = "wall "+tokens(f);
            if startsWith(tokens(f),"T")
                b.family = "T";
                b.normalDirections = repmat(normals(f,:),size(line,1),1);
            end
            basis(end+1) = b; %#ok<AGROW>
        end
        contacts(end+1,:) = q(k,:)-[0,0,1]; %#ok<AGROW>
    end
    neighbours = [pairs(pairs(:,1)==k,2);pairs(pairs(:,2)==k,1)];
    neighbours = sort(neighbours);
    for j = neighbours.'
        if opt.N_image == 0 || ~any(families(2,:)), continue; end
        line = image_line_points(q(k,:),q(j,:),1,Rp,opt.N_image, ...
            opt.imageOffsetFraction,1);
        % T_t points toward this neighbour; n/n2 form a deterministic frame
        % perpendicular to that line, also for nonhorizontal particle pairs.
        t = (q(j,:)-q(k,:))/norm(q(j,:)-q(k,:));
        [~,axisIndex] = min(abs(t));
        axis = zeros(1,3); axis(axisIndex) = 1;
        n = cross(t,axis); n = n/norm(n);
        n2 = cross(t,n);
        normals = [zeros(2,3);n;n2;t;zeros(1,3)];
        % Preserve the default [pair S;pair D;pair R] coefficient order.
        for f = [1,6,2,5,3,4]
            if ~families(2,f), continue; end
            b = block; b.family = tokens(f); b.points = line;
            b.label = "pair "+j+" "+tokens(f);
            if startsWith(tokens(f),"T")
                b.family = "T";
                b.normalDirections = repmat(normals(f,:),size(line,1),1);
            end
            basis(end+1) = b; %#ok<AGROW>
        end
        contacts(end+1,:) = q(k,:)+(q(j,:)-q(k,:))/norm(q(j,:)-q(k,:)); %#ok<AGROW>
    end
    for j = 1:size(contacts,1)
        rout_block = [rout_block;collocation_surface_patch( ...
            sphere,contacts(j,:).',cap)]; %#ok<AGROW>
    end
    rout_block = uniquetol(rout_block,1e-12,'ByRows',true,'DataScale',1);
    rin_block = vertcat(basis.points);
    if any(rin_block(:,3) <= 0) || any(vecnorm(rin_block-q(k,:),2,2) >= 1)
        error('wall:ImagesGrid','All sources must be inside their sphere and above the wall.');
    end
    problem.M(k) = size(rout_block,1);
    problem.N_total(k) = size(rin_block,1);
    if problem.M(k) < 2*problem.N_total(k)
        error('wall:ImagesOversampling', ...
            'Body %d: M=%d needs at least 2*N_total=%d; increase the initial rout or capCount.', ...
            k,problem.M(k),2*problem.N_total(k));
    end
    % Kin' maps mixed coefficients to moments. Reflected wall sources are
    % prescribed by the kernels and never enter this particle moment sum.
    Kin = zeros(3*problem.N_total(k),6); powers = zeros(size(Kin,1),1);
    offset = 0;
    for j = 1:numel(basis)
        count = size(basis(j).points,1); cols = offset+(1:3*count);
        if basis(j).family == "S"
            Kin(cols,:) = getKmat(basis(j).points,q(k,:));
        elseif basis(j).family == "R"
            Kin(cols,4:6) = repmat(eye(3),count,1);
        end
        if basis(j).family == "T" || basis(j).family == "R", powers(cols) = 1; end
        if basis(j).family == "D", powers(cols) = 2; end
        offset = offset+3*count;
    end
    Kout = getKmat(rout_block,q(k,:));
    % Six columns suffice for both completion and the LL projection action.
    completion = Kin/(Kin.'*Kin);
    body = struct('basis',basis,'contacts',contacts,'Kin',Kin,'Kout',Kout, ...
        'completion',completion,'rin',rin_block,'rout',rout_block);
    problem.body{k} = body;
    problem.rin = [problem.rin;rin_block];
    problem.rout = [problem.rout;rout_block];
    problem.powers = [problem.powers;powers];
    problem.sourceUnits = [problem.sourceUnits;opt.mu*a.^(1+powers)];
    problem.sourceOffsets(k+1) = problem.sourceOffsets(k)+size(Kin,1);
    problem.targetOffsets(k+1) = problem.targetOffsets(k)+size(Kout,1);
    % The dense path needs these physical maps, but no local factorizations.
    if opt.SVDsolve, continue; end
    % Reuse only genuinely identical complete blocks, including pair lines.
    reuse = 0;
    for j = 1:k-1
        other = problem.body{j};
        if q(k,3)==q(j,3) && isequal(size(rin_block),size(other.rin)) && ...
                isequal(size(rout_block),size(other.rout)) && ...
                isequal([basis.family],[other.basis.family]) && ...
                isequal({basis.normalDirections},{other.basis.normalDirections}) && ...
                norm(rin_block-q(k,:)-(other.rin-q(j,:)),inf)<1e-12 && ...
                norm(rout_block-q(k,:)-(other.rout-q(j,:)),inf)<1e-12
            reuse = j; break
        end
    end
    if reuse
        problem.Y{k} = problem.Y{reuse}; problem.UU{k} = problem.UU{reuse};
        problem.ranks(k) = problem.ranks(reuse);
        problem.singularValues{k} = problem.singularValues{reuse};
        continue
    end
    if opt.gmres_verbose
        fprintf('%s body %d: one-body SVD %d x %d (%.1f MiB per matrix).\n', ...
            mode,k,size(Kout,1),size(Kin,1),8*size(Kout,1)*size(Kin,1)/2^20);
    end
    t = tic;
    Sblock = wall_image_velocity_matrix(basis,rout_block,kernelOpt);
    if strcmp(mode,'mobility')
        % Remove the moment-carrying field, then store rigid motion in it.
        Sblock = Sblock+(Kout-Sblock*completion)*Kin.';
    end
    assemblySeconds = assemblySeconds+toc(t);
    t = tic;
    [Y,UU,Si] = getPseudoFactors(Sblock,opt.tol,0);
    if isempty(Si), error('wall:ImagesRank','SVD cutoff removed every one-body mode.'); end
    problem.Y{k} = Y; problem.UU{k} = UU.';
    problem.ranks(k) = numel(Si); problem.singularValues{k} = Si;
    factorSeconds = factorSeconds+toc(t);
end
% Gather each family once, keeping its coefficient columns in body ordering.
problem.kernel = struct('family',{},'cols',{},'rin',{},'geometry',{}, ...
    'normalDirections',{});
for family = ["S","T","D","R"]
    cols = zeros(0,1); points = zeros(0,3);
    normals = zeros(0,3);
    for k = 1:P
        offset = problem.sourceOffsets(k);
        basis = problem.body{k}.basis;
        for j = 1:numel(basis)
            count = 3*size(basis(j).points,1);
            if basis(j).family == family
                cols = [cols;offset+(1:count).']; %#ok<AGROW>
                points = [points;basis(j).points]; %#ok<AGROW>
                if family == "T"
                    normals = [normals;basis(j).normalDirections]; %#ok<AGROW>
                end
            end
            offset = offset+count;
        end
    end
    if isempty(cols), continue; end
    kernel = struct('family',family,'cols',cols,'rin',points,'geometry',[], ...
        'normalDirections',normals);
    if family == "S"
        kernel.geometry = blake.prepare(points.',problem.rout.',0, ...
            logical(opt.fmm && ~opt.SVDsolve));
    end
    problem.kernel(end+1) = kernel;
end
problem.timing = struct('assemblySeconds',assemblySeconds, ...
    'factorizationSeconds',factorSeconds,'preparationSeconds',toc(timer));
end
