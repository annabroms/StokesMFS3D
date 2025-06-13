function [rvec_in,rvec_out,opt] = init_spheres(q,Rp,N)
%INIT_SPHERES(q,Rp,N) 

% Set default values if Rg and N are not provided
if nargin < 2
    Rp = 0.68; %proxy radius
    N = 700; % approximate number of proxy sources on every particle
elseif nargin < 3
    N = 700;
end

P = size(q,1); 

% initialize a bunch of parameters. Do not change if you don't really want.
opt = init_MFS(N);
%opt = init_MFS(1000);
opt.Rp = Rp;
opt.fmm = 0; 
opt.maxit = 200; %max number of gmres iterations
opt.gmres_tol = 1e-7;
opt.plot = 0; %visualise?

%% Discretize one body
% Inner proxy surface, outer collocation grid for single sphere using spherical design nodes
[rin,rout] =  getDesignGrid(opt.Rp,opt); 


rvec_in = [];
rvec_out = [];
 

for k = 1:P
    rvec_in = [rvec_in; rin+q(k,:)];
    rvec_out = [rvec_out; rout+q(k,:)];
     
end

end