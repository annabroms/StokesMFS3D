%test sqrt of TWG and the effect of small negative eigenvalues of this
%matrix

clear; close all;
rng(5); 

sphere = 1;
tdesign = 0; %uniform discretisation of sphere
dense = 1; %build C densely? (sym(TG))
fmm = 0;
fast = 1;% different implementation of lanczos with reduced storage
plot_eigs = 1; 
use_sym = 0;


%% Set and create geometry
P = 1;
delta = 5;
E0 = [.4 .6 1]; %Type T in the mobility paper
sep = 0.125;
%sep = 0.2;

if sphere
    E0 = [1 1 1]; % this is a sphere
    sep = 0.9;
    sep = 0.3; 
end

% Set resolution (following ellipsoid_mobility_run)
Nv = 40;
Nv = 15; % need a much coarser resolution



rng(8);
[~, R, qvec, ~] = ellipsoid_cluster(E0,P,delta);
R{1} = eye(3); 
disp('Geometry created')


visualise = 1;
show_dir = 1; %visualise eigenvector/singvector corresponding to smallest eigval / singval?


a = 1.5; %Determines oversampling factor for the collocation points
a = 2; 

% Get discretization and normals
N1 = a*0.75*Nv;
N2 = a*Nv;
if sphere && tdesign
    Nt = 336;
    q = cell2mat(cellfun(@transpose, qvec, 'UniformOutput', false).');
    Rp = 0.5;
    [rvec_in,rvec_out,opt] = init_spheres(q,Rp,Nt,a);
    M = length(rvec_out)/P;
    n_out = rvec_out - kron(q,ones(M,1));
    w = 4*pi/M*ones(P*M,1);
    Rprecond = [];
       
else
    [rvec_in,rvec_out,~,~,~,w,n_out] = getEllipsoidGrids(E0,P,delta,N1,N2,sep,R,qvec); %Assign source and collocation points
    Rprecond = R;

end

M = length(rvec_out)/P;
N = length(rvec_in)/P;


% Visualise geometry
figure(10)
scatter3(rvec_in(:,1),rvec_in(:,2),rvec_in(:,3));
hold on
scatter3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3));
axis equal

% Noise
dW2 = randn(3*N*P,1);

%% block-diagonal preconditioner 
% Heavy for many particles!
[C,TW,A] = getSymTWS(rvec_out,rvec_in,n_out,w);
 %% How close is C to it's SPD counterpart? (keeping only positive eigenvalues)
[Ve,De] = eig(C);
d = diag(De); 

ind = find(d>0);
d_plus = zeros(3*N*P,1);
d_plus(ind) = d(ind);
ind = setdiff(1:3*N*P,ind);
d_minus = zeros(3*N*P,1);
d_minus(ind) = d(ind);  

C_plus_sqrt = Ve*diag(sqrt(d_plus));
C_minus_sqrt = Ve*diag(sqrt(d_minus));
C_sqrt = Ve*diag(sqrt(d)); 

%Sanity check
C_plus_rec = C_plus_sqrt*C_plus_sqrt';
C_minus_rec = C_minus_sqrt*C_minus_sqrt.'; %want the transpose here and not the complex conjugate
C_minus = Ve*diag(d_minus)*Ve';
C_plus = Ve*diag(d_plus)*Ve';

C_rec = C_plus_rec+C_minus_rec;
norm(C_rec-C)
figure()
subplot(1,6,1)
imagesc(log10(abs(C)));
title('Symmetric matrix TG')
colorbar
subplot(1,6,2)
imagesc(log10(abs(C_rec)));
colorbar
title('Recovered matrix, after taking sqrt and outer product')
subplot(1,6,3)
imagesc(log10(abs(C_rec-C)));
colorbar
subplot(1,6,4); 
imagesc(C_plus_rec);
colorbar
subplot(1,6,5); 
imagesc(log10(abs(C_minus_rec)));
colorbar
subplot(1,6,6)
imagesc(log10(abs(C_plus_rec-C)));
colorbar

%% Eigendecomp for single body
Ci = getSymTWS(rvec_out(1:M,:),rvec_in(1:N,:),n_out(1:M,:),w(1:M));
[Vi,Di] = eig(Ci);
di = diag(Di); 

m = {'+','d','.','s'};
c = {'b','r','k','c','g'};

tol = 1e-12; 
maxit = 500; 

vars.fmm = fmm; 
vars.eps = 1e-10; %for use in FMM 
tolvec = [1e-16 0];
tolvec = 1e-10;

dense_res1 = C_sqrt*dW2;
dense_res2 = sqrtm(C)*dW2; %should be the same?
%dense_res3 = sqrtm(C_rec)*dW2; %should be the same?
%dense_res4 = chol(C)*dW2; %should only be the same in expectation.
dense_res5 = sqrtm(A)*dW2;

%dense_res = C_sqrt*dW2; %should be ok with tdesign disc

% Loop over different truncation levels of the eigvals to check the effect
for i = 1:length(tolvec)

    % Build preconditioner as in Bao et al
    trunc = tolvec(i); 
    ind = find(real(di)>trunc);
    diff_set = setdiff(1:length(di),ind);  
    dsqrt = sqrt(di);
    dsqrt_inv = 1./dsqrt;
    dsqrt_inv(diff_set) = 0;
    Pi = diag(dsqrt_inv)*Vi';
    dsqrt_plus  = dsqrt;
    dsqrt_plus(diff_set) = 0;
    Pi_plus = Vi*diag(dsqrt_plus);

    % dsqrt_plus = sqrt(d); 
    % dsqrt_plus(diff_set) = 0;
    % Ci_plus = Ve*diag(dsqrt_plus);    
    
    PPplus = @(x) apply_block_diag(x, Pi_plus, P, Rprecond);

    if dense
        Pglob = kron(eye(P),Pi);
        if use_sym
            C_precond = Pglob*C*Pglob';
        else
            C_precond = Pglob*A*Pglob';
        end

        if plot_eigs
            e_precond = eig(C_precond); 
            e = eig(C); 
            %Eigenvalues are clustered at 1 with preconditioning! 
            figure(2)
            plot(real(e_precond),imag(e_precond),'+')
            hold on 
            plot(real(e),imag(e),'o')
        end
        CBC = @(x) C_precond*x;

        a = rand(N*P*3,1);
        pre = apply_block_diag(a, Pi_plus', P, Rprecond);
        post = CBC(pre);
        res1 = apply_block_diag(post, Pi_plus, P, Rprecond);
        res2 = C*a; %supposed to be the same as res1

    else
        % Build matrix via matvec
         %Mfun = @(x) getPrecondTG(x,P,rvec_in,rvec_out,n_out,w,Pi,vars);
         %[noise,iter_errv,iters] = lanczosSqrt_fast(Mfun,dW,tol_lanczos,maxit,PPplus);

        CBC = @(x) getPrecondTG(x,P,rvec_in,rvec_out,n_out,w,Pi,vars,Rprecond);
    end
    if fast
        [y,iter_err] = lanczosSqrt_fast(CBC,dW2,tol,maxit,PPplus);
        [y2,iter_err2] = arnoldiSqrt_fast(CBC,dW2,tol,maxit,PPplus);
    else
        [y,iter_err] = lanczosSqrt(CBC,dW2,tol,maxit,PPplus);
        %[y,iter_errv5] = lanczosSqrt(CBC,dW2,tol,maxit,[]);
    end
   % y = Pi_plus*y;% debug only
    y_hist(:,i) = y;    

    %Use the non-symmetrized version of A - does NOT work
    % A_precond = C*A*C';
    % [y,iter_errv6] = lanczosSqrt(A_precond,dW2,tol,maxit,C);
    % figure(11)
    % semilogy(iter_errv6,'Marker',m{mod(i,4)+1},'Color',c{mod(i,5)+1},'DisplayName',num2str(tolvec(i)));
    % hold on
    
end

lgd = legend('show','Location','best');
lgd.Title.String = 'Eigenvalue truncation (preconditioner)';

axis tight
grid on
xlabel('Iteration number')
ylabel('Estimated error using Lanczos')
