function opt = init_MFS(des_n)
%INIT_MFS(des_n) initializes parameters for MFS spheres with close to des_n source points
%
% %   NOTE: Spherical designs do not exist for all integers, so the true
%         number of proxy points N will be close to des_n
if nargin<1    
    opt.des_n = 700;
else
    opt.des_n = des_n;
end
opt.tol = 1e-14; %default tolerance for the truncated SVD
opt.plot = 0; % debug parameter

opt.maxit = 400;  %max number of GMRES iterations
opt.gmres_tol = 1e-6; %GMRES tol

opt.profile = 0; %activate profiling


opt.remove_last = 0; %remove last singular value in SVD?
opt.image = 0; %default not to use images

opt.eps = 1e-8; %used for FMM
%opt.eps = 1e-10; %change tol in fmm to see effect

opt.fmm = 0; 
opt.ellipsoid = 0; %spheres if false 

opt.get_traction = 0; 
opt.a_glob = 1.2; %M (number of collocation point) > N (number of sources on proxy surface)

%% discretisation with images: 
%NB: to be updated
opt.beta = 3; %number of singularities per image point
opt.sing_vec = [1 1 0 1]; %source types at image points when images are used
opt.lines = 1; 
opt.cone = 0;
opt.cheb = 1; 
opt.np = 1; 
opt.theta_max_out = [pi/5; pi/50];
opt.double_patch = 1; 


end