function [Y,UU,LL,Kin,Kout] = getSL(rin,rout,q)
%getSL(rin,rout) Compute factors for the pseudoinverse for the mobility 
% matrix of a single particle, given discretisation 
if nargin<3   
    q = [0 0 0]; %The particle is assumed to not be rotated
end

S = generate_stokes_matrix_vec_inout(rin,rout);
Kout = getKmat(rout,q);
Kin = getKmat(rin,q);

%Get projection matrix: 
% tic
% for k = 1:200
% LL = getProjection(Kin,rin,q);
% end
% toc
% 
% %KtK = getKtK(rin,q); 
% 
% %Lr = Kout*Kin';
% 
% tic
% for k = 1:200
% LL2 = Kin*((Kin'*Kin)\Kin');
% end
% toc
Lr = Kout*Kin';
LL = Kin*((Kin'*Kin)\Kin'); 

%Form 1-body matrix explicitly
A = (S*(eye(size(LL))-LL)+Lr);

%cond(A) %or skeel(A) for debug.


[U,SS,V] = svd(A,"econ");
SS = diag(SS);
%tol = 1e-14; 
%ra = sum(SS>max(SS)*tol); 
%iS = 1./SS(1:ra); % rank
%Y = V(:,1:ra)*diag(iS); 
%UU = U(:,1:ra)'

iS = 1./SS; % rank

Y = V*diag(iS); 
UU = U';


end