function [C,TW,A] = getSymTWS(rvec_out,rvec_in,n_out,w,varargin)
%GETSYMTWS Form the symmetrized weighted T*S matrix needed for the fluctuating 
% velocity field for Brownian motion with MFS.
%
%   C = GETSYMTWS(rvec_out,rvec_in,n_out,w)
%   C = GETSYMTWS(rvec_out,rvec_in,n_out,w,n_in,opt)
%
%   INPUTS:
%       rvec_out - M x 3 proxy/source points.
%       rvec_in  - N x 3 collocation/target points.
%       n_out    - M x 3 normal or dipole direction vectors attached to rvec_out.
%       w        - Quadrature weights for rvec_out, length M.
%       n_in     - Optional N x 3 inner proxy normals for opt.add_rank1.
%       opt      - Optional struct with add_rank1 and rank1_scale.
%
%   OUTPUTS: 
%       C
%       TW       - Transpose of W'T', with W quadrature weights and T' the 
%                  traction operator - for debugging
%       A        - C = (A+A')/2

if nargin == 0
    self_test;
    return;
end

nin = [];
opt = struct();
for iarg = 1:numel(varargin)
    arg = varargin{iarg};
    if isstruct(arg)
        opt = arg;
    elseif ~isempty(arg)
        nin = arg;
    end
end
if ~isfield(opt,'add_rank1')
    opt.add_rank1 = false;
end
if ~isfield(opt,'rank1_scale')
    opt.rank1_scale = 1;
end

T = stokes_DLP_mat(rvec_out,rvec_in,n_out);
S = stokes_SLP_mat(rvec_in,rvec_out);
w = w(:).';
ww = repmat(w,3,1);
W = diag(ww(:));
TW = T*W;
A = TW*S;
if opt.add_rank1
    if isempty(nin)
        if isfield(opt,'nin') && ~isempty(opt.nin)
            nin = opt.nin;
        elseif isfield(opt,'n_in') && ~isempty(opt.n_in)
            nin = opt.n_in;
        end
    end
    if isempty(nin)
        error('getSymTWS:missingInnerNormals', ...
            'opt.add_rank1=true requires inner proxy normals.');
    end
    A = A + normal_rank1_block(nin,opt.rank1_scale);
end
C = (A+A')/2;

end

function self_test()
rng(2);
close all; 

E0 = [0.7 0.45 1.1];
%E0 = [1 1 1];
P = 1;
delta = 0.12;
N1 = 24;
N2 = 16;
sep = 0.15;

[~, R, qvec, ~] = ellipsoid_cluster(E0,P,delta);
R{1} = eye(3); 
[rvec_in, rvec_out, ~, ~, ~, w, n_out] = getEllipsoidGrids(E0,P,delta,N1,N2,sep,R,qvec);

figure()
scatter3(rvec_in(:,1),rvec_in(:,2),rvec_in(:,3));
hold on
scatter3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3));

% q = [0 0 0];
% [rvec_in,rvec_out,] = init_spheres(q,0.3,100,1.2);
% n_out = rvec_out; 
% M = size(rvec_out,1); 
% w = 4*pi/M*ones(M,1);

[~,TW] = getSymTWS(rvec_out,rvec_in,n_out,w);

Kout = getKmat(rvec_out,qvec{1}');
Kin = getKmat(rvec_in,qvec{1}');
diff_mat = Kout'*TW'-Kin';
N = size(rvec_in,1);
lambda = rand(3*N,1); 
F1= Kout'*TW'*lambda;
F2 = Kin'*lambda;
ft_rel_err = norm(F1-F2)/norm(F1);

sym_err = norm(diff_mat,'fro') / max(norm(Kout,'fro'),eps);
fprintf('getSymTWS self-test:\n');
fprintf('  symmetry error = %.3e\n', sym_err);
fprintf('  force/torque relative error = %.3e\n', ft_rel_err);


end
