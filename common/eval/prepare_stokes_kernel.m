function kernel = prepare_stokes_kernel(rin,rout,opt)
%PREPARE_STOKES_KERNEL Prepare fixed NP-by-3 sources and MP-by-3 targets.
% opt.kernel is blake or freeSpace; backend, mu, wall_height and
% fmmTolerance are explicit. Blake preparation caches geometry, not trees.
if nargin == 0, self_test(); kernel = []; return; end
validateattributes(rin,{'double'},{'real','finite','2d','ncols',3});
validateattributes(rout,{'double'},{'real','finite','2d','ncols',3});
kernel = struct('rin',rin,'rout',rout,'options',opt,'geometry',[]);
assert(ismember(string(opt.kernel),["blake","freeSpace"]));
if string(opt.kernel) == "blake"
    assert(~isempty(which('blake.prepare')),'wall:MissingBlake', ...
        'Add StokesHalfspaceFMM/matlab to the MATLAB path.');
    kernel.geometry = blake.prepare(rin.',rout.',opt.wall_height, ...
        string(opt.backend) == "fmm");
else
    assert(~any(ismember(rin,rout,'rows')), ...
        'wall:Coincidence','A source coincides with a target.');
end
end

function self_test()
% Verify unequal source/target counts and prepared Blake geometry.
o = struct('kernel',"blake",'backend',"blake",'mu',1, ...
    'wall_height',0,'fmmTolerance',1e-12);
k = prepare_stokes_kernel([0,0,1],[0.2,0.3,0;0.1,0,2],o);
assert(k.geometry.ns == 1 && k.geometry.nt == 2);
fprintf('prepare_stokes_kernel: PASS; %d source, %d targets, wall height %.3g.\n', ...
    k.geometry.ns,k.geometry.nt,k.options.wall_height);
end
