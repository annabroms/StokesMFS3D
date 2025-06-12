function opt = init_MFS(N)
%INIT_MFS(N) initializes parameters for MFS spheres with N source points

opt.fib = 0; % if using the fibonacci grid on both proxy and collocation surfaces
if nargin<1
    opt.des_n = 700;
else
    opt.des_n = N;
end
% if nargin < 2
%     opt.fib_n = nph;
%     nph = [];
%     nth = [];
%     opt.fib = 1; % if using the fibonacci grid on both surfaces
% else
%     opt.fib_n = nph*nth;
%     
% end
% 
% %Default values for MFS for spheres
% opt.nph = nph;
% opt.nth = nth;
opt.stokes = 2;
opt.theta = 0.5; 
opt.ver = 1; % Set outer radius to the radius of the bc, i.e. 1
opt.inout = 1; %always true for MFS
%ss = read_mobility(1);
%opt.ss = ss; %if we want to compare to true solution for a pair of particles
opt.vel_check = 0; 
opt.quad = 0; %use quadrature weights to compute integrals
opt.trans = 0; 
opt.tol = 1e-14; %default tolerance for the truncated SVD
opt.plot = 0; % debug parameter
opt.beta = 3; %number of singularities per image point
opt.sing_vec = [1 1 0 1]; %source types at image points when images are used
opt.maxit = 400;  %max number of GMRES iterations
opt.gmres_tol = 1e-6; %GMRES tol

opt.profile = 0; 

opt.dense = 0; 
opt.mapped = 0; 
opt.precond = 0; 
opt.split_precond = 0; 

opt.remove_last = 0; %remove last singular value in SVD?
opt.image = 0; %default not to use images
opt.design = 1; %use spherical designs
opt.fib = 0; %fibonacci grid instead?

opt.eps = 1e-8; %used for FMM
%opt.eps = 1e-10; %change tol in fmm to see effect
opt.ellipsoid = 0; 
opt.fmm = 0; 
opt.disc = 0;
opt.get_traction = 0; 
opt.a_glob = 1.2; %M (number of collocation point) > N (number of sources on proxy surface)

%% discretisation with images
opt.lines = 1; 
opt.cone = 0;
opt.cheb = 1; 
opt.np = 1; 
opt.theta_max_out = [pi/5; pi/50];
opt.double_patch = 1; 



end