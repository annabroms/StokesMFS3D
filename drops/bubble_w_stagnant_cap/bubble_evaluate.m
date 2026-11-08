function field = bubble_evaluate(solution,points,opt,normals)
%BUBBLE_EVALUATE Batched physical exterior fields in the laboratory frame.
% Pass normals to obtain traction and its tangential projection as well.
% Only evaluate exterior/surface targets; no interior bubble flow is defined.
if nargin==0
    self_test();
    field = [];
    return
end
validateattributes(points,{'numeric'},{'2d','ncols',3,'real','finite'});
count = size(points,1);
withTraction = nargin>=4;
if withTraction
    assert(isequal(size(normals),size(points)) && ...
        all(abs(vecnorm(normals,2,2)-1)<1e-11));
    field.traction = zeros(count,3);
    field.tangentTraction = zeros(count,3);
end
assert(all(vecnorm(points-opt.centre,2,2)>=opt.radius*(1-1e-11)), ...
    'bubble_evaluate:InteriorTarget','Bubble field targets must be exterior or on the surface.');
field.velocity = zeros(count,3);
for first = 1:opt.batchSize:count
    rows = first:min(first+opt.batchSize-1,count);
    field.velocity(rows,:) = reshape(stokes_mfs_velocity_matrix( ...
        solution.basis,points(rows,:))*(solution.lambda/opt.muOut),3,[]).';
    if withTraction
        t = reshape(stokes_mfs_traction_matrix(solution.basis, ...
            points(rows,:),normals(rows,:))*solution.lambda,3,[]).';
        field.traction(rows,:) = t;
        field.tangentTraction(rows,:) = t-sum(t.*normals(rows,:),2).*normals(rows,:);
    end
end
field.relativeVelocity = field.velocity-solution.U(:).';
if withTraction
    field.normalVelocity = sum(field.relativeVelocity.*normals,2);
    field.tangentVelocity = field.relativeVelocity-field.normalVelocity.*normals;
end
end

function self_test()
% Verify viscosity, interleaving, batching, force sign and far-field decay.
o = bubble_options(struct('radius',1.2,'centre',[0.2 -0.3 0.4], ...
    'muOut',2.3,'batchSize',2));
unit = [eye(3);-eye(3)];
points = o.centre+o.radius*unit;
s.basis = struct('family',"S",'points',o.centre+[0.1 0.1 0], ...
    'normalDirections',zeros(0,3),'scale',1);
s.lambda = [0.3;-0.4;1];
s.U = [0.1;-0.2;0.3];
f = bubble_evaluate(s,points,o,unit);
u = reshape(drop_velocity_block(s.basis.points,points,o.muOut)*s.lambda,3,[]).';
t = reshape(drop_traction_block(s.basis.points,points,unit)*s.lambda,3,[]).';
assert(norm(f.velocity-u,'fro')<1e-13 && norm(f.traction-t,'fro')<1e-13);
assert(norm(f.relativeVelocity-(u-s.U.'),'fro')<1e-13);
assert(max(abs(sum(f.tangentTraction.*unit,2)))<1e-13);
near = bubble_evaluate(s,o.centre+[10 0 0],o);
far = bubble_evaluate(s,o.centre+[100 0 0],o);
assert(norm(far.velocity)<0.11*norm(near.velocity));
try
    bubble_evaluate(s,o.centre,o);
    error('bubble_evaluate:SelfTest','Interior target accepted.');
catch ME
    assert(strcmp(ME.identifier,'bubble_evaluate:InteriorTarget'));
end
fprintf('bubble_evaluate self-test passed.\n');
end
