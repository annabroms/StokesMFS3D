function [rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a_glob,a)
%INIT_SPHERES Ordinary proxy and surface grids for equal spheres.
% [rin,rout,opt] = init_spheres(q,Rp,N,a_glob,a)
% q is P-by-3 physical centres; a is the physical radius (default 1).
% Rp is the proxy radius divided by a; N requests proxy nodes and a_glob
% requests the ratio M/N of surface to proxy nodes. Spherical design counts
% can round down. No wall-specific geometry initializer is needed: wall
% clearance is checked by the wall solver, not by this grid generator.
% The third output retains legacy solver defaults and records opt.radius=a.

% Set default values if Rg and N are not provided
if nargin < 2
    Rp = 0.68; %proxy radius
    N = 700; % approximate number of proxy sources on every particle
    a_glob = 1.2;
elseif nargin < 3
    N = 700;
    a_glob = 1.2; 
elseif nargin < 4
    a_glob = 1.2;
end
if nargin < 5, a = 1; end
validateattributes(a,{'double'},{'scalar','positive','finite'});

P = size(q,1); 

% initialize a bunch of parameters. Do not change if you don't really want.
opt = init_MFS(N);
opt.a_glob = a_glob;
opt.lr = 0; %no deflation based long-range preconditioning
opt.Rp = Rp;
opt.radius = a;
opt.fmm = 0; 
opt.maxit = 200; %max number of gmres iterations
opt.gmres_tol = 1e-7;
opt.plot = 0; %visualise?
opt.debug = 0; % determine system matrix to look at eigenvalues?

%% Discretize one body
% Inner proxy surface, outer collocation grid for single sphere using spherical design nodes
[rin,rout] =  getDesignGrid(opt.Rp,opt); 


rvec_in = zeros(P*size(rin,1),3);
rvec_out = zeros(P*size(rout,1),3);
 

for k = 1:P
    rvec_in((k-1)*size(rin,1)+(1:size(rin,1)),:) = a*rin+q(k,:);
    rvec_out((k-1)*size(rout,1)+(1:size(rout,1)),:) = a*rout+q(k,:);
     
end

end
