function y = drop_image_apply_A(z,grid,factors,opt)
%DROP_IMAGE_APPLY_A Apply the original scaled rectangular coefficient map.
if nargin == 0
    self_test(); y = []; return
end
P = numel(grid.body);
M = arrayfun(@(b) size(b.collocation,1),grid.body);
y = zeros(6*sum(M),1);
lambda = cell(P,1);
row0 = 0; col0 = 0;
for i = 1:P
    ne = size(factors{i}.H,2);
    ni = size(factors{i}.selfMatrix,2)-ne;
    zi = z(col0+(1:ne+ni));
    y(row0+(1:6*M(i))) = factors{i}.selfMatrix*zi;
    alpha = zi(1:ne);
    lambda{i} = alpha-factors{i}.C*(factors{i}.H*alpha);
    row0 = row0+6*M(i); col0 = col0+ne+ni;
end
if col0~=numel(z), error('drop_image_apply_A:BadLength','Bad coefficient length.'); end
y = y+drop_image_interactions(lambda,grid,factors,opt);
end

function self_test()
% Compare a complete mixed-image action with explicit global assembly.
q = [-1.025 0 0;1.025 0 0]; a = ones(2,1);
opt = drop_image_options(struct('N',30,'collocationRatio',2, ...
    'N_image',1,'imageFamilies',"SDRT_nT_t",'mirror_ext',true, ...
    'capPointCounts',[12 12]),2,a);
g = drop_image_discretization(q,a,opt);
f = cell(2,1); blocks = cell(2,2);
for i = 1:2
    [blocks{i,i},H,t1,t2] = drop_image_self_block( ...
        g.body(i),2,opt.completionScale(i));
    f{i}.selfMatrix = blocks{i,i}; f{i}.H = H;
    f{i}.C = H'/(H*H');
    f{i}.tangentOne = t1; f{i}.tangentTwo = t2;
end
for i = 1:2
    j = 3-i;
    blocks{i,j} = drop_image_interaction_block( ...
        g.body(j),g.body(i),f{i}.tangentOne,f{i}.tangentTwo,2);
end
A = [blocks{1,1},blocks{1,2};blocks{2,1},blocks{2,2}];
z = sin((1:size(A,2)).');
y = drop_image_apply_A(z,g,f,opt);
assert(norm(y-A*z,inf)<1e-10);
fprintf('drop_image_apply_A self-test passed.\n');
end
