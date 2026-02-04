%test sqrt of traction computation

clear; close all;
rng(5); 

P = 500; %number of bodies
delta = 0.1; %smallest particle particle distance 
if P == 2
    q = [0 0 0; 2+delta 0 0]; %center coordiante matrix for P particles, x,y,z: size P x 3
elseif P == 1
    q = [0 0 0]; 
else
    %random configurations
    [q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
end


  
fmm = 1; %only activate if many particles (say, more than 40)

%Test first with very low resolution
Rp = 0.30;
N = 100; 
% Rp = 0.15;
% N = 50; 
a = 1.2; 
a = 2; %or play with SVD truncation level

dense = 0; %build full system matrices (costly)
plot_eigs = 0; 

%Rp = 0.68; %proxy radius
%N = 700; % approximate number of proxy sources on every particle


[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a);
M = size(rvec_out,1)/P; 
N = size(rvec_in,1)/P; 
n_out = repmat(rvec_out(1:M,:),P,1);

if dense
    T = stokes_DLP(rvec_out,rvec_in,n_out);
    S = generate_stokes_mat(rvec_in, rvec_out);
    A  = T*S;
    B = (A+A')/2;
end

dW2 = rand(3*N*P,1);

%% Block block-diagonal preconditioner 
Ti = stokes_DLP_mat(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:));
Si = generate_stokes_mat(rvec_in(1:N,:), rvec_out(1:M,:));
Ai = Ti*Si;
Bi = (Ai+Ai')/2;
[Ve,De] = eig(Bi);
d = diag(De); 

tolvec = logspace(-15,-3,13);
m = {'+','d','.','s'};
c = {'b','r','k','c','g'};

tol = 1e-8; 
maxit = 500; 

vars.fmm = fmm; 
vars.eps = 1e-10; %for use in FMM 
tolvec = 1e-12;

% Loop over different truncation levels of the eigvals to check the effect
for i = 1:length(tolvec)

    % Build preconditioner as in Bao et al
    trunc = tolvec(i); 
    ind = find(real(d)>trunc);
    diff_set = setdiff(1:length(d),ind);    
    dsqrt = 1./sqrt(d);
    dsqrt(diff_set) = 0;
    Ci = diag(dsqrt)*Ve';

    dsqrt_plus = sqrt(d); 
    dsqrt_plus(diff_set) = 0;
    Ci_plus = Ve*diag(dsqrt_plus);    
    Cplus = @(x) apply_block_diag(x, Ci_plus, P);

    if dense
        C = kron(eye(P),Ci);
        B_precond = C*B*C';
        if plot_eigs
            e_precond = eig(B_precond); 
            e = eig(B); 
            %Eigenvalues are clustered at 1 with preconditioning! 
            figure(2)
            plot(real(e_precond),imag(e_precond),'+')
            hold on 
            plot(real(e),imag(e),'o')
        end
        CBC = @(x) B_precond*x;

    else
        CBC = @(x) getPrecondTG(x,P,rvec_in,rvec_out,n_out,Ci,vars);
    end
    
    [y,iter_errv5] = KrylovSqrtMsing(CBC,dW2,tol,maxit,Cplus);
    
    figure(10)
    semilogy(iter_errv5,'Marker',m{mod(i,4)+1},'Color',c{mod(i,5)+1}, ...
        'DisplayName',sprintf('trunc tol = %.1e', tolvec(i)));
    hold on

    %Use the non-symmetrized version of A - does NOT work
    % A_precond = C*A*C';
    % [y,iter_errv6] = KrylovSqrtMsing(A_precond,dW2,tol,maxit,C);
    % figure(11)
    % semilogy(iter_errv6,'Marker',m{mod(i,4)+1},'Color',c{mod(i,5)+1},'DisplayName',num2str(tolvec(i)));
    % hold on
    
end

function y = apply_block_diag(x, B, P)
%APPLY_BLOCK_DIAG Apply block-diagonal operator with P identical blocks B.
% x is (3N*P) x 1, B is (3N) x (3N).
nb = size(B,1);
y = zeros(size(x));
for k = 1:P
    idx = (k-1)*nb + (1:nb);
    y(idx) = B * x(idx);
end
end

lgd = legend('show','Location','best');
lgd.Title.String = 'Eigenvalue truncation (preconditioner)';

axis tight
grid on
xlabel('Iteration number')
ylabel('Estimated error using Lanczos')
