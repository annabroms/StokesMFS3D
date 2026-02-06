%test sqrt of TG with preconditioning for 2 bodies. This script in addition displays
%convergence with the RPY tensor and without preconditioning.

clear; close all;

P = 2; %number of bodies
delta = 10; %smallest particle particle distance 
if P == 2
    q = [0 0 0; 2+delta 0 0]; %center coordiante matrix for P particles, x,y,z: size P x 3
else
    q = [0 0 0]; 
end

%random configurations
%[q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
  
fmm = 0; %only activate if many particles (say, more than 40)

%Test first with very low resolution
Rp = 0.30;
N = 100; 
% Rp = 0.15;
% N = 50; 
a = 1.2; 
a = 2; %or play with SVD truncation level

%Rp = 0.68; %proxy radius
%N = 700; % approximate number of proxy sources on every particle


[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a);
M = size(rvec_out,1)/P; 
N = size(rvec_in,1)/P; 
n = repmat(rvec_out(1:M,:),P,1);

Tt = getTractionMat(rvec_in,rvec_out,n); 
G = stokes_SLP_mat(rvec_in, rvec_out);

A  = Tt'*G;

%% Check convergence with unpreconditioned Lanczos
%Do eigval decomp
B = (A+A')/2;
[V,D] = eig(B);
d = diag(D);
trunc = 1e-4; 
ind = find(real(d)>trunc);
diff_set = setdiff(1:length(d),ind);
d2 = d;
d2(diff_set) = 0;
inv_V = V\eye(size(V));
A2 = V*diag(d2)*inv_V;
A3 = V*diag(d2)*V';
% D = chol(real(A2)); 


% yexact = D*dW; 
S = generate_RPY_matrix(rvec_out,0.01);
dW1 = rand(size(S,1),1);
dW2 = rand(size(A,1),1);
maxit = 500; 
tol = 1e-8; 
[~,iter_errv1] = KrylovSqrtMsing(S,dW1,tol,maxit);
[y,iter_errv2] = KrylovSqrtMsing(A,dW2,tol,maxit);
[~,iter_errv3] = KrylovSqrtMsing(A+A',dW2,tol,maxit);
[~,iter_errv4] = KrylovSqrtMsing(A3,dW2,tol,maxit);
y2 = krylov_sqrt(A,dW2,20);
norm(y-y2,inf)/norm(y,inf)

%%
figure()
semilogy(iter_errv1,'+-');
hold on
semilogy(iter_errv2,'o-');
semilogy(iter_errv3,'*-');
semilogy(iter_errv4,'s--');
xlabel('Iteration number')
ylabel('Estimated error using Lanczos')
axis tight
legend('RPY tensor','TtG','Symmetrized TtG')

%% Block diagonal preconditioner
Tt = getTractionMat(rvec_in(1:N,:),rvec_out(1:M,:),rvec_out(1:M,:)); 
S = stokes_SLP_mat(rvec_in(1:N,:), rvec_out(1:M,:));
Ai = Tt'*S;
Bi = (Ai+Ai')/2;
%Bi = Ai; 
[Ve,De] = eig(Bi);
d = diag(De); 

tolvec = logspace(-15,-3,13);
m = {'+','d','.','s'};
c = {'b','r','k','c','g'};
% Loop over different tolerances to check the effect
for i = 1:length(tolvec)

    trunc = tolvec(i); 
    ind = find(real(d)>trunc);
    diff_set = setdiff(1:length(d),ind);
    %d2 = sqrt(d(ind));
    
    dsqrt = 1./sqrt(d);
    dsqrt(diff_set) = 0;
    Ci = diag(dsqrt)*Ve';
    
    dsqrt_plus = sqrt(d); 
    dsqrt_plus(diff_set) = 0;
    Ci_plus = Ve*diag(dsqrt_plus);
    
    
    %Ci = diag(d2)*Ve(:,ind)';
    
    C = kron(eye(P),Ci);
    Cplus = kron(eye(P),Ci_plus);
    B_precond = C*B*C';
    
    e_precond = eig(B_precond); 
    e = eig(B); 

    %% Display eigvals of preconditioned matrix
    %Eigenvalues are clustered at 1! 
    figure(2)
    plot(real(e_precond),imag(e_precond),'+')
    title('Eigvals of onebody preconditioned matrix')
    
    
    % semilogy(sort(e_precond));
    % hold on
    % semilogy(sort(e))
    %% Run lanczos with this truncation level in the preconditioner
    [y,iter_errv5] = KrylovSqrtMsing(B_precond,dW2,tol,maxit,Cplus);
    
    figure(10)
    semilogy(iter_errv5,'Marker',m{mod(i,4)+1},'Color',c{mod(i,5)+1},'DisplayName',num2str(tolvec(i)));
    hold on

    %Use the non-symmetrized version of A - does NOT work
    % A_precond = C*A*C';
    % [y,iter_errv6] = KrylovSqrtMsing(A_precond,dW2,tol,maxit,C);
    % figure(11)
    % semilogy(iter_errv6,'Marker',m{mod(i,4)+1},'Color',c{mod(i,5)+1},'DisplayName',num2str(tolvec(i)));
    % hold on
    
end

legend show

axis tight
grid on
xlabel('Iteration number')
ylabel('Estimated error using Lanczos')






