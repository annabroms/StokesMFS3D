function G = stokes_kernel_matrix(rin,rout,opt)
%STOKES_KERNEL_MATRIX Dense direct reference or one-body velocity block.
% Global production interactions must use apply_stokes_kernel instead.
if nargin == 0, self_test(); G = []; return; end
if string(opt.kernel) == "freeSpace"
    assert(~any(ismember(rin,rout,'rows')),'Source and target coincide.');
    G = stokes_SLP_mat(rin,rout)/opt.mu;
    return
end
G = zeros(3*size(rout,1),3*size(rin,1));
opts = struct('method','blake','mu',opt.mu,'wall_height',opt.wall_height);
for j = 1:size(rin,1)
    for component = 1:3
        force = zeros(3,1); force(component) = 1;
        u = blake.velocity(rin(j,:).',force,rout.',opts);
        G(:,3*(j-1)+component) = u(:);
    end
end
end

function self_test()
% Compare dense assembly against the independent GGV direct representation.
o = struct('kernel',"blake",'backend',"ggv-direct",'mu',1, ...
    'wall_height',0,'fmmTolerance',1e-12);
s = [0,0,1;0.2,0.1,0.7]; t = [0.2,-0.1,0;0.4,0.3,2];
f = (1:6).'; k = prepare_stokes_kernel(s,t,o);
G = stokes_kernel_matrix(s,t,o);
assert(norm(G*f-apply_stokes_kernel(f,k)) < 1e-13);
fprintf('stokes_kernel_matrix: PASS; %d-by-%d block, direct/GGV action error %.3e.\n', ...
    size(G,1),size(G,2),norm(G*f-apply_stokes_kernel(f,k)));
end
