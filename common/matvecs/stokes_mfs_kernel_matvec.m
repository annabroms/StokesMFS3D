function result = stokes_mfs_kernel_matvec(x,problem,mode)
%STOKES_MFS_KERNEL_MATVEC Preserve physical Blake interactions and matched correction.
% x is target-space data, length 3MP. z=Bpre^+ x is auxiliary source data.
% The correction uses Gpre, retaining the physical wall self field when
% free-space preconditioning is selected. Mobility K*H terms cancel here.
if nargin == 0 
    self_test(); result = []; 
    return;
end
P = problem.P; 
N = problem.N; 
M = problem.M;
assert(numel(x) == 3*M*P && iscolumn(x));
lambda = zeros(3*N*P,1); self = zeros(size(x));
for p = 1:P
    ns = (p-1)*3*N+(1:3*N); 
    ms = (p-1)*3*M+(1:3*M);
    b = problem.body{p}; 
    factor = b.(mode);
    z = factor.Y*(factor.U.'*x(ms));
    if strcmp(mode,'mobility') 
        z = z-b.C*(b.H*z); 
    end
    lambda(ns) = z;
    self(ms) = b.Gpre*z;
end

result = x+apply_stokes_kernel(lambda,problem.kernel)-self;

end

function self_test()
% Independently assemble both square operators for unequal-height bodies.
q = [0,0,1.4;3,0,2];
for pre = ["wall","freeSpace"]
    [s,t,o] = init_spheres_wall(q,1,struct('N',18,'backend',"blake", ...
        'preconditioner',pre,'computeResidual',false));
    p = prepare_wall_problem(q,s,t,o);
    G = stokes_kernel_matrix(p.rin,p.rout,p.kernel.options);
    x = sin((1:3*p.M*p.P).');
    for mode = {'resistance','mobility'}
        expected = x;
        for j = 1:p.P
            cols = (j-1)*3*p.N+(1:3*p.N);
            rows = (j-1)*3*p.M+(1:3*p.M);
            b = p.body{j}; f = b.(mode{1}); z = f.Y*(f.U.'*x(rows));
            if strcmp(mode{1},'mobility'), z = z-b.C*(b.H*z); end
            expected = expected+G(:,cols)*z;
            expected(rows) = expected(rows)-b.Gpre*z;
        end
        actual = stokes_mfs_kernel_matvec(x,p,mode{1});
        assert(norm(actual-expected)/norm(expected)<1e-12);
        fprintf('stokes_mfs_kernel_matvec: PASS; %s, %s, relative action error %.3e.\n', ...
            pre,mode{1},norm(actual-expected)/norm(expected));
    end
end
end
