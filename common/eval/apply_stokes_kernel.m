function [u,info] = apply_stokes_kernel(lambda,kernel)
%APPLY_STOKES_KERNEL Map interleaved point forces to interleaved velocities.
if nargin == 0, self_test(); u = []; info = []; return; end
assert(iscolumn(lambda) && numel(lambda) == 3*size(kernel.rin,1));
assert(isreal(lambda) && all(isfinite(lambda)));
o = kernel.options;
if string(o.kernel) == "blake"
    opts = struct('method',char(o.backend),'mu',o.mu, ...
        'wall_height',o.wall_height,'tol',o.fmmTolerance);
    [velocity,info] = blake.apply(kernel.geometry,reshape(lambda,3,[]),opts);
    u = velocity(:);
elseif string(o.backend) == "fmm"
    vars = struct('fmm',true,'fmm_tol',o.fmmTolerance);
    u = getStokesletFlow(lambda,kernel.rin,kernel.rout,vars)/o.mu;
    info = struct();
else
    % Bounded target batches avoid a global free-space interaction matrix.
    u = zeros(3*size(kernel.rout,1),1);
    for first = 1:128:size(kernel.rout,1)
        nodes = first:min(first+127,size(kernel.rout,1));
        rows = reshape(3*nodes-2+(0:2).',[],1);
        u(rows) = stokes_SLP_mat(kernel.rin,kernel.rout(nodes,:))*lambda/o.mu;
    end
    info = struct();
end
assert(numel(u) == 3*size(kernel.rout,1) && all(isfinite(u)));
end

function self_test()
% Compare all force directions with the independent direct Blake evaluator.
o = struct('kernel',"blake",'backend',"blake",'mu',2, ...
    'wall_height',0,'fmmTolerance',1e-12);
s = [0,0,1;0.2,-0.1,0.8]; t = [0.1,0.4,0;0.3,0.2,2];
k = prepare_stokes_kernel(s,t,o); f = (1:6).';
u = apply_stokes_kernel(f,k);
ref = blake.velocity(s.',reshape(f,3,[]),t.', ...
    struct('method','blake','mu',2));
assert(norm(u-ref(:)) < 1e-13);
assert(norm(apply_stokes_kernel(2*f,k)-2*u) < 1e-13);
assert(norm(u(1:3)) < 1e-13);
fprintf('apply_stokes_kernel: PASS; reference error %.3e, linearity error %.3e, wall speed %.3e.\n', ...
    norm(u-ref(:)),norm(apply_stokes_kernel(2*f,k)-2*u),norm(u(1:3)));
end
