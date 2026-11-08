function y = drop_image_interactions(lambda,grid,factors,opt)
%DROP_IMAGE_INTERACTIONS Apply all off-diagonal mixed-source drop blocks.
if nargin == 0
    self_test(); y = []; return
end
P = numel(grid.body);
M = arrayfun(@(b) size(b.collocation,1),grid.body);
y = zeros(6*sum(M),1);
coeff = cell(P,1);
for j = 1:P
    coeff{j} = image_source_coefficients(grid.body(j).basis,lambda{j});
end
row0 = 0;
for i = 1:P
    body = grid.body(i);
    n = body.normals;
    normal = zeros(M(i),1);
    tangent = zeros(M(i),2);
    for first = 1:opt.batchSize:M(i)
        rows = first:min(first+opt.batchSize-1,M(i));
        x = body.collocation(rows,:);
        v = zeros(numel(rows),3);
        t = zeros(numel(rows),3);
        for j = 1:P
            if j==i, continue; end
            v = v+direct_image_velocity(grid.body(j).basis,coeff{j},x);
            one.body.basis = grid.body(j).basis;
            t = t+image_traction(one,lambda(j),x,n(rows,:),1);
        end
        normal(rows) = sum(v.*n(rows,:),2);
        tangent(rows,1) = sum(t.*factors{i}.tangentOne(rows,:),2);
        tangent(rows,2) = sum(t.*factors{i}.tangentTwo(rows,:),2);
    end
    w = body.rowScale;
    y(row0+3*M(i)+(1:M(i))) = w.*normal/opt.velocityScale;
    traction = (1-opt.interiorViscosity(i))* ...
        opt.tractionRowWeight/opt.tractionScale* ...
        (w.*tangent);
    y(row0+4*M(i)+(1:2*M(i))) = reshape(traction.',[],1);
    row0 = row0+6*M(i);
end
end

function self_test()
% Match explicit pair blocks with all four image-source families.
q = [-1.025 0 0;1.025 0 0]; a = [1;1];
opt = drop_image_options(struct('N',30,'collocationRatio',2, ...
    'N_image',1,'imageFamilies',"SRTD",'mirror_ext',true, ...
    'capPointCounts',[12 12]),2,a);
g = drop_image_discretization(q,a,opt);
f = cell(2,1); lambda = cell(2,1);
for j = 1:2
    [~,H,t1,t2] = drop_image_self_block(g.body(j),2,1/(8*pi));
    f{j}.tangentOne = t1; f{j}.tangentTwo = t2;
    z = sin((1:size(H,2)).'); C = H'/(H*H');
    lambda{j} = z-C*(H*z);
end
y = drop_image_interactions(lambda,g,f,opt);
o = 0; expected = zeros(size(y));
for i = 1:2
    j = 3-i;
    block = drop_image_interaction_block(g.body(j),g.body(i), ...
        f{i}.tangentOne,f{i}.tangentTwo,2);
    local = block*[lambda{j};zeros(3*sum(arrayfun( ...
        @(b) size(b.points,1),g.body(j).interiorBasis)),1)];
    m = size(g.body(i).collocation,1);
    local(1:4*m) = local(1:4*m)/opt.velocityScale;
    local(4*m+1:end) = ...
        local(4*m+1:end)*opt.tractionRowWeight/opt.tractionScale;
    expected(o+(1:6*m)) = local; o = o+6*m;
end
assert(norm(y-expected,inf)<1e-10);
fprintf('drop_image_interactions self-test passed.\n');
end
