function grid = drop_image_discretization(q,a,opt)
%DROP_IMAGE_DISCRETIZATION Add outer proxy and optional Kelvin image lines.
if nargin == 0
    self_test(); grid = struct(); return
end
grid = build_image_discretization(q,opt,a);
P = size(q,1);
families = image_family_list(opt.imageFamilies);
template = struct('family',"S",'points',zeros(0,3), ...
    'normalDirections',zeros(0,3),'scale',1);
outerPoints = cell(P,1);
outerDirections = cell(P,1);
for body = 1:P
    outerPoints{body} = zeros(0,3);
    outerDirections{body} = zeros(0,3);
end
lineOffset = zeros(P,1);
for pair = 1:size(grid.pairs,1)
    i = grid.pairs(pair,1); j = grid.pairs(pair,2);
    N = grid.pairs(pair,4);
    if ~opt.useImages || ~opt.mirror_ext || N==0, continue; end
    for body = [i,j]
        nodes = lineOffset(body)+(1:N);
        inner = grid.body(body).imagePoints(nodes,:);
        direction = grid.body(body).imageDirections(nodes(1),:);
        rstar = norm(inner(1,:)-q(body,:));
        [line,~] = exterior_image_line_points(q(body,:),direction, ...
            a(body),a(body)*opt.outerProxyRatio,rstar,N, ...
            opt.outerImageOffsetFraction);
        outerPoints{body} = [outerPoints{body};line];
        outerDirections{body} = [outerDirections{body}; ...
            repmat(direction,N,1)];
        lineOffset(body) = lineOffset(body)+N;
    end
end
for body = 1:P
    grid.body(body).normals = ...
        (grid.body(body).collocation-q(body,:))/a(body);
    innerProxy = grid.body(body).basis(1).points;
    unit = (innerProxy-q(body,:))/(a(body)*opt.proxyRadius);
    b = template;
    b.points = q(body,:)+a(body)*opt.outerProxyRatio*unit;
    basis = b;
    if opt.useImages && opt.mirror_ext && ~isempty(outerPoints{body})
        for family = families
            b = template;
            b.family = string(family);
            b.points = outerPoints{body};
            if ismember(family,["T","T_n","T_n2","T_t"])
                b.normalDirections = image_stresslet_directions( ...
                    outerDirections{body},opt.use_orth_dir_images,family);
            elseif family=="Q_s"
                b.normalDirections = outerDirections{body};
                b.scale = a(body)^2;
            end
            basis(end+1,1) = b; %#ok<AGROW>
        end
    end
    grid.body(body).interiorBasis = basis;
    grid.body(body).exteriorImagePoints = outerPoints{body};
    rIn = vecnorm(vertcat(grid.body(body).basis.points)-q(body,:),2,2);
    rOut = vecnorm(vertcat(basis.points)-q(body,:),2,2);
    if any(rIn>=a(body)) || any(rOut<=a(body))
        error('drop_image_discretization:WrongSide', ...
            'An exterior-flow or interior-flow source is on the wrong side.');
    end
    ncol = 3*sum(arrayfun(@(s) size(s.points,1), ...
        [grid.body(body).basis;basis]));
    if 6*size(grid.body(body).collocation,1)<=ncol
        error('drop_image_discretization:Undersampled', ...
            'Increase collocation nodes or decrease source count.');
    end
end
end

function self_test()
% Check unequal radii, contact counts, and both source-side inequalities.
q = [0 0 0;3.05 0 0]; a = [1;2];
opt = drop_image_options(struct('N',30,'collocationRatio',2, ...
    'N_image',3,'mirror_ext',true,'capPointCounts',[20 20], ...
    'imageFamilies',"SRTD",'use_orth_dir_images',true),2,a);
g = drop_image_discretization(q,a,opt);
assert(size(g.body(1).exteriorImagePoints,1)==3);
assert(size(g.body(2).exteriorImagePoints,1)==3);
for k = 1:2
    assert(all(vecnorm(g.body(k).exteriorImagePoints-q(k,:),2,2)>a(k)));
    assert(abs(norm(g.body(k).exteriorImagePoints(end,:)-q(k,:)) - ...
        a(k)*opt.outerProxyRatio*(1-opt.outerImageOffsetFraction))<1e-12);
    innerT = find(arrayfun(@(b) b.family=="T",g.body(k).basis),1);
    outerT = find(arrayfun( ...
        @(b) b.family=="T",g.body(k).interiorBasis),1);
    ray = g.body(k).imageDirections;
    innerDirection = g.body(k).basis(innerT).normalDirections;
    outerDirection = g.body(k).interiorBasis(outerT).normalDirections;
    assert(max(abs(sum(innerDirection.*ray,2)))<1e-14);
    assert(norm(innerDirection-outerDirection,'fro')<1e-14);
end
opt.use_orth_dir_images = false;
parallel = drop_image_discretization(q,a,opt);
for k = 1:2
    innerT = find(arrayfun(@(b) b.family=="T",parallel.body(k).basis),1);
    assert(norm(parallel.body(k).basis(innerT).normalDirections- ...
        parallel.body(k).imageDirections,'fro')<1e-14);
end
opt.imageFamilies = "SDRT_nT_t";
joint = drop_image_discretization(q,a,opt);
for k = 1:2
    bases = {joint.body(k).basis,joint.body(k).interiorBasis};
    ray = joint.body(k).imageDirections;
    for side = 1:2
        basis = bases{side};
        if numel(basis)~=6 || ...
                max(abs(sum(basis(5).normalDirections.*ray,2)))>1e-14 || ...
                norm(basis(6).normalDirections-ray,'fro')>1e-14
            error('drop_image_discretization:SelfTest','Joint stresslet directions failed.');
        end
    end
end
opt.imageFamilies = "SDRT_tT_nT_n2Q_s";
full = drop_image_discretization(q,a,opt);
for k = 1:2
    bases = {full.body(k).basis,full.body(k).interiorBasis};
    ray = full.body(k).imageDirections;
    for side = 1:2
        basis = bases{side};
        first = basis(5).normalDirections;
        second = basis(6).normalDirections;
        assert(max(abs(sum(first.*second,2)))<1e-14);
        assert(max(abs(sum(second.*ray,2)))<1e-14);
        assert(norm(basis(7).normalDirections-cross(ray,second,2),'fro')<1e-14);
        assert(norm(basis(8).normalDirections-ray,'fro')<1e-14);
        assert(basis(8).scale==a(k)^2);
    end
end
opt.use_orth_dir_images = true;
opt.outerImageOffsetFraction = 0.005;
closer = drop_image_discretization(q,a,opt);
assert(norm(closer.body(1).exteriorImagePoints(end,:)-q(1,:)) > ...
    norm(g.body(1).exteriorImagePoints(end,:)-q(1,:)));
opt.mirror_ext = false;
g = drop_image_discretization(q,a,opt);
assert(all(arrayfun(@(b) isscalar(b.interiorBasis),g.body)));
fprintf('drop_image_discretization self-test passed.\n');
end
