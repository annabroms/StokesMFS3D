function [maxError,pointError,points,normals] = ...
    drop_image_validate_surface_conditions(grid,lambda,gamma,U,S1,e,opt,points,normals)
%DROP_IMAGE_VALIDATE_SURFACE_CONDITIONS Check all three physical BCs off grid.
if nargin == 0
    self_test(); maxError=[]; pointError={}; points={}; normals={}; return
end
P = numel(grid.body);
if nargin<8 || isempty(points)
    [points,normals] = drop_image_validation_grid( ...
        grid,opt.validationCount,opt.validationGapCount,max(opt.capAngles));
end
if isscalar(S1), S1 = repmat(S1,P,1); end
if isequal(size(e),[1 3]), e = repmat(e,P,1); end
e = e./vecnorm(e,2,2);
maxError = zeros(P,3); pointError = cell(P,1);
for i = 1:P
    x = points{i}; n = normals{i};
    pointError{i} = zeros(size(x,1),3);
    for first = 1:opt.batchSize:size(x,1)
        rows = first:min(first+opt.batchSize-1,size(x,1));
        xx = x(rows,:); nn = n(rows,:);
        [u0,~,t0] = drop_image_field(xx,nn,0,grid,lambda,gamma,opt);
        [ui,~,ti] = drop_image_field(xx,nn,i,grid,lambda,gamma,opt);
        dv = u0-ui;
        normalError = abs(sum((u0-U(i,:)).*nn,2));
        jump = t0-ti;
        tangent = jump-sum(jump.*nn,2).*nn;
        d = nn*e(i,:).';
        target = S1(i)*d.*(e(i,:)-d.*nn);
        pointError{i}(rows,:) = [vecnorm(dv,2,2),normalError, ...
            vecnorm(tangent-target,2,2)];
    end
    maxError(i,:) = max(pointError{i},[],1);
end
end

function self_test()
% Zero flow and zero active stress satisfy all surface conditions.
q = [0 0 0]; a = 1;
opt = drop_image_options(struct('N',30,'useImages',false, ...
    'validationCount',50,'validationGapCount',21),1,a);
g = drop_image_discretization(q,a,opt);
lam = {zeros(3*size(g.body(1).basis(1).points,1),1)};
gam = {zeros(3*size(g.body(1).interiorBasis(1).points,1),1)};
[err,point,x] = drop_image_validate_surface_conditions( ...
    g,lam,gam,zeros(1,3),0,[0 0 1],opt);
assert(all(err==0) && all(point{1}==0,'all') && ~isempty(x{1}));
fprintf('drop_image_validate_surface_conditions self-test passed.\n');
end
