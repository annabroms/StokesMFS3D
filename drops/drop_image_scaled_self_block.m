function block = drop_image_scaled_self_block(body,mu,c,opt)
%DROP_IMAGE_SCALED_SELF_BLOCK Build one row-scaled drop matrix without SVD.
% Shared by the iterative preconditioner and the dense global reference.
if nargin == 0
    self_test(); block = struct(); return
end
[A,H,t1,t2] = drop_image_self_block(body,mu,c);
M = size(body.collocation,1);
w = body.rowScale;
row = [repelem(w,3);w;repelem(w,2)];
row(1:4*M) = row(1:4*M)/opt.velocityScale;
row(4*M+1:end) = ...
    row(4*M+1:end)*opt.tractionRowWeight/opt.tractionScale;
block.selfMatrix = row.*A;
block.H = H;
block.C = H'/(H*H');
block.tangentOne = t1;
block.tangentTwo = t2;
block.rowScale = row;
end

function self_test()
% Match the existing pure-Stokeslet block after row scaling.
v = [eye(3);-eye(3)];
x = [1 1 1;1 -1 -1;-1 1 -1;-1 -1 1]/sqrt(3);
b = struct('family',"S",'points',0.6*v, ...
    'normalDirections',zeros(0,3),'scale',1);
body.basis = b;
b.points = 2*v; body.interiorBasis = b;
body.collocation = x; body.normals = x; body.rowScale = ones(4,1);
opt = drop_image_options(struct('N',30,'velocityScale',2, ...
    'tractionScale',3,'tractionRowWeight',1.5),1,1);
block = drop_image_scaled_self_block(body,2,1/(8*pi),opt);
A = drop_recompleted_matrix(0.6*v,2*v,x,x,2,1/(8*pi));
A(1:16,:) = A(1:16,:)/2;
A(17:end,:) = A(17:end,:)/2;
assert(norm(block.selfMatrix-A,'fro')<1e-12);
assert(isequal(size(block.selfMatrix),[24 36]));
fprintf('drop_image_scaled_self_block self-test passed.\n');
end
