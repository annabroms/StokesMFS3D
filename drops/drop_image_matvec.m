function y = drop_image_matvec(eta,grid,factors,opt)
%DROP_IMAGE_MATVEC Apply identity-diagonal one-body preconditioned operator.
if nargin == 0
    self_test(); y = []; return
end
P = numel(grid.body);
M = arrayfun(@(b) size(b.collocation,1),grid.body);
if numel(eta)~=6*sum(M), error('drop_image_matvec:BadLength','Bad data length.'); end
lambda = cell(P,1);
o = 0;
for i = 1:P
    zi = applyOneBodyPinv(eta(o+(1:6*M(i))),factors{i});
    ne = size(factors{i}.H,2);
    alpha = zi(1:ne);
    lambda{i} = alpha-factors{i}.C*(factors{i}.H*alpha);
    o = o+6*M(i);
end
y = eta+drop_image_interactions(lambda,grid,factors,opt);
end

function self_test()
% Verify the range-discard identity with both image representations.
q = [-1.025 0 0;1.025 0 0]; a = ones(2,1);
opt = drop_image_options(struct('N',30,'collocationRatio',2, ...
    'N_image',1,'imageFamilies',"SDRT_nT_t",'mirror_ext',true, ...
    'capPointCounts',[12 12]),2,a);
g = drop_image_discretization(q,a,opt);
f = cell(2,1); M = arrayfun(@(b) size(b.collocation,1),g.body);
for i = 1:2
    f{i} = drop_image_one_body_preconditioner( ...
        g.body(i),2,opt.completionScale(i),opt);
end
eta = sin((1:6*sum(M)).');
z = []; discarded = zeros(size(eta)); o = 0;
for i = 1:2
    rows = o+(1:6*M(i));
    z = [z;applyOneBodyPinv(eta(rows),f{i})]; %#ok<AGROW>
    L = f{i}.leftTranspose;
    discarded(rows) = eta(rows)-L.'*(L*eta(rows));
    o = o+6*M(i);
end
assert(norm(drop_image_apply_A(z,g,f,opt)- ...
    (drop_image_matvec(eta,g,f,opt)-discarded),inf)<1e-9);
fprintf('drop_image_matvec self-test passed.\n');
end
