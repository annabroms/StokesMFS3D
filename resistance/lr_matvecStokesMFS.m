function res = lr_matvecStokesMFS(x,rvec_in,rvec_out,q,UU,Y,opt,R,Rinv,AN,AM)
%LR_MATVECSTOKESMFS applies longrange preconditioning to the matrix 
% corresponding to one body preconditioning for the 2D Stokes resistance problem 

%Apply the preconditioner P, see applyPmat.m
vel = matvecStokesMFS(x,rvec_in,rvec_out,q,UU,Y,opt,1,R);

res = applyPmat(vel,rvec_in,rvec_out,Rinv,AN,AM,opt);


end