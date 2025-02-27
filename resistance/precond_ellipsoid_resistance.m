function [Y,UU] = precond_ellipsoid_resistance(rin,rout)
%PRECOND_ELLIPSOID_RESISTANCE Compute factors for the pseudoinverse for a single particle in a Stokes resistance problem, given source distribution rin and distribution of collocation point rout. 

S = generate_stokes_matrix_vec_inout(rin,rout);

% [L,U] = lu(S);
% return

%Look at singular values
% [U,SS,V] = svd(S);
% %[Q,R] = qr(S); 
% 
% SS = diag(SS);
% tol = 1e-14; 
% ra = sum(SS>max(SS)*tol); 
% iS = 1./SS(1:ra); % rank
% 
% figure(55)
% % clf;
% semilogy(SS,'o-');
% hold on
% semilogy(ra*ones(1,2),logspace(-15,5,2),'r--')
% c = max(SS)/min(SS);  %Condition number
% str = sprintf('Self condition number %1.3e, max sing %1.3e',c,max(SS));
% title(str,'interpreter','latex');
% grid on
% xlabel('$j$','interpreter','latex')
% ylabel('$\sigma_j$','interpreter','latex')
% 
% 
% Y = V(:,1:ra)*diag(iS); 
% UU = U(:,1:ra)';

%UU = U';
[U,SS,V] = svd(S,"econ");
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
