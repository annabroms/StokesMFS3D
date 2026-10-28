function uvec = apply_wall_images_velocity(lambda,problem,rout,direct)
%APPLY_WALL_IMAGES_VELOCITY Batched complete-wall mixed-source velocity action.
%   uvec = apply_wall_images_velocity(lambda,problem)
%   uvec = apply_wall_images_velocity(lambda,problem,rout,direct)
% lambda is dimensionless, interleaved and body/family-grouped as documented
% in prepare_wall_images_problem. Optional rout is dimensionless; omitting
% it reuses prepared collocation geometry. direct=true selects direct S as
% well as T_t,D,R, for independent small references and surface diagnostics.
%
% Hybrid backend: ONE global half-space S FMM call; ONE blocked direct call
% per optional family. All S blocks, including wall/pair line S, are batched.
% No pair loop, dense interaction matrix or column assembly occurs here.
% T_t,D,R are direct-only in the installed StokesHalfspaceFMM package.
% Package strengths are -lambda for T and -2*lambda for D. Each stresslet
% keeps its prepared source orientation: axial or transverse to its own
% wall/pair line. Geometry is normalized to wall=0 and viscosity=1.

prepared = nargin < 3;
if prepared, rout = problem.rout; end
if nargin < 4, direct = false; end
opt = problem.options;
uvec = zeros(3*size(rout,1),1);
if ~any(lambda), return; end
options = struct('mu',1,'wall_height',0,'method','blake', ...
    'block_size',[opt.batchSize,256]);
for j = 1:numel(problem.kernel)
    kernel = problem.kernel(j);
    strength = reshape(lambda(kernel.cols),3,[]);
    if ~any(strength,'all'), continue; end
    switch kernel.family
        case "S"
            o = options;
            if opt.fmm && ~direct && ~opt.SVDsolve
                o.method = 'fmm'; o.tol = opt.fmm_tol;
            end
            if prepared
                u = blake.apply(kernel.geometry,strength,o);
            else
                u = blake.velocity(kernel.rin.',strength,rout.',o);
            end
        case "T"
            o = options; o.normals = kernel.normalDirections.';
            u = blake.stresslet(kernel.rin.',-strength,rout.',o);
        case "D"
            u = blake.potential_dipole(kernel.rin.',-2*strength,rout.',options);
        case "R"
            u = blake.rotlet(kernel.rin.',strength,rout.',options);
    end
    uvec = uvec+u(:);
end
end
