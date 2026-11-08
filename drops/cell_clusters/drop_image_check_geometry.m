function report = drop_image_check_geometry(grid,opt)
%DROP_IMAGE_CHECK_GEOMETRY Audit paired image rays, radii and source ownership.
% Radii/clearances are in physical units; KelvinError and RayError are
% dimensionless. NaN denotes an absent outer line or a one-point spacing.
if nargin == 0
    self_test(); report = table(); return
end
names = {'Drop','Neighbor','NInner','NOuter','RInnerFirst','RInnerEnd', ...
    'ROuterFirst','ROuterEnd','InnerClearance','OuterClearance', ...
    'MinOuterSpacing','KelvinError','RayError'};
data = nan(2*size(grid.pairs,1),numel(names));
offset = zeros(numel(grid.body),1);
row = 0;
tol = 1e-11;

% Check every source block against its owner's physical boundary.
for i = 1:numel(grid.body)
    q = grid.centres(i,:); a = grid.radii(i);
    body = grid.body(i);
    rin = vecnorm(vertcat(body.basis.points)-q,2,2);
    rout = vecnorm(vertcat(body.interiorBasis.points)-q,2,2);
    if any(rin>=a) || any(rout<=a)
        error('drop_image_check_geometry:WrongSide','A source crosses its owner boundary.');
    end
    for b = 2:numel(body.basis)
        if ~isequal(body.basis(b).points,body.imagePoints)
            error('drop_image_check_geometry:Ownership','Inner basis and image nodes differ.');
        end
    end
    for b = 2:numel(body.interiorBasis)
        if ~isequal(body.interiorBasis(b).points,body.exteriorImagePoints)
            error('drop_image_check_geometry:Ownership','Outer basis and image nodes differ.');
        end
    end
end

% Read lines in the same pair-major order, but test their geometry directly.
for pair = 1:size(grid.pairs,1)
    ids = grid.pairs(pair,1:2); N = grid.pairs(pair,4);
    if ~opt.useImages || N==0, continue; end
    pairOffset = offset;
    for side = 1:2
        i = ids(side); j = ids(3-side);
        q = grid.centres(i,:); a = grid.radii(i);
        ray = grid.centres(j,:)-q; ray = ray/norm(ray);
        nodes = offset(i)+(1:N);
        x = grid.body(i).imagePoints(nodes,:)-q;
        rin = x*ray.';
        rayError = max(vecnorm(x-rin*ray,2,2))/a;
        if any(rin<=0 | rin>=a) || any(diff(rin)>=0)
            error('drop_image_check_geometry:InnerLine','Invalid inner radial ordering.');
        end
        % The two accumulation points are mutually related by inversion.
        rj = norm(grid.body(j).imagePoints(pairOffset(j)+1,:)-grid.centres(j,:));
        d = norm(grid.centres(j,:)-q);
        if abs(rin(1)*(d-rj)/a^2-1)>tol
            error('drop_image_check_geometry:Accumulation','Incorrect accumulation point.');
        end
        row = row+1;
        data(row,:) = [i,j,N,0,rin(1),rin(end),NaN,NaN, ...
            a-max(rin),NaN,NaN,NaN,rayError];
        if opt.mirror_ext
            x = grid.body(i).exteriorImagePoints(nodes,:)-q;
            rout = x*ray.';
            rayError = max(rayError,max(vecnorm(x-rout*ray,2,2))/a);
            kelvinError = abs(rin(1)*rout(1)/a^2-1);
            R = a*opt.outerProxyRatio;
            endpoint = R*(1-opt.outerImageOffsetFraction);
            if any(rout<=a | rout>=R) || any(diff(rout)<=0) || ...
                    kelvinError>tol || (N>1 && abs(rout(end)-endpoint)>tol*a)
                error('drop_image_check_geometry:OuterLine', ...
                    'Invalid Kelvin point, exterior radial bounds, ordering or endpoint.');
            end
            spacing = NaN;
            if N>1, spacing = min(diff(rout)); end
            data(row,[4,7,8,10,11,12,13]) = ...
                [N,rout(1),rout(end),min(rout)-a,spacing,kelvinError,rayError];
        end
        if rayError>tol
            error('drop_image_check_geometry:Ray','An image node is off its contact ray.');
        end
        offset(i) = offset(i)+N;
    end
end
for i = 1:numel(grid.body)
    if size(grid.body(i).imagePoints,1)~=offset(i) || ...
            size(grid.body(i).exteriorImagePoints,1)~=opt.mirror_ext*offset(i)
        error('drop_image_check_geometry:Count','Paired image counts are inconsistent.');
    end
end
report = array2table(data(1:row,:),'VariableNames',names);
end

function self_test()
% Exercise unequal spheres, zero/one/many nodes and a corrupted outer ray.
q = [0 0 0;3.05 0 0]; a = [1;2];
opt = drop_image_options(struct('N',30,'N_image',3, ...
    'collocationRatio',2,'capPointCounts',[20 20], ...
    'mirror_ext',true,'imageFamilies',"SDRT_nT_t"),2,a);
for N = [0 1 3]
    opt.N_image = N;
    g = drop_image_discretization(q,a,opt);
    report = drop_image_check_geometry(g,opt);
    if height(report)~=2*(N>0) || any(report.NOuter~=N) || ...
            any(report.KelvinError>1e-13)
        error('drop_image_check_geometry:SelfTest','Incorrect geometry report.');
    end
end
g.body(1).exteriorImagePoints(1,2) = 0.01;
for b = 2:numel(g.body(1).interiorBasis)
    g.body(1).interiorBasis(b).points = g.body(1).exteriorImagePoints;
end
rejected = false;
try
    drop_image_check_geometry(g,opt);
catch err
    rejected = strcmp(err.identifier,'drop_image_check_geometry:Ray');
end
if ~rejected, error('drop_image_check_geometry:SelfTest','Off-ray node accepted.'); end
fprintf('drop_image_check_geometry self-test passed.\n');
end
