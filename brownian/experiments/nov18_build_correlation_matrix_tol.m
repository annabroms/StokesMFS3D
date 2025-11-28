clear; close all; 
%Builds correlation matrix for <UU'> when the noise is given by the 
% fluctuating sources with correlation TG. This script checks the error in
% the correlation matrix <UU^T> relative to the known mobility matrix. 

% set particle coordinates
P=2;
delta = 1; 
if P==1
    q = [0 0 0];
else
    delta = 3; %initial particle particle distance
    q = [0 0 0; 2+delta 0 0];
end

% Set resolution
% Rp = 0.63; %fine resolution
% N = 700;

% Rp = 0.15; %Proxy sphere radius -- very coarse resolution
% N = 50;  % Number of proxy sources
% % 
Rp = 0.30;
N = 100; 

a = 2; %Determines oversampling factor for the collocation points

%Discretise particles 
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points
N = size(rvec_in,1)/P;
MM = size(rvec_out,1)/P; 

%Get all matrices needed to solve the system densely
B = getKmat(rvec_out,q);
K = getKmat(rvec_in,q);

n = repmat(rvec_out(1:MM,:),P,1);
if P == 1
    T = getTraction(rvec_in,rvec_out,n); 
else
    Tblock = getTraction(rvec_in(1:N,:),rvec_out(1:MM,:),rvec_out(1:MM,:));
    T = kron(eye(P),Tblock);
end
S = generate_stokes_mat(rvec_in, rvec_out);
visualise = 1; 
tol = 1e-13;

A = -4*pi/MM*T*S;

[Y,U]  = getPseudoFactors(S,tol,visualise);

R = K'*Y*(U'*B);
M = R\eye(P*6); 

[V,D] = eig(A);
d = diag(D);

if P == 2
    M_true = M;
else
    M_true = diag([1/6/pi*ones(1,3) 1/8/pi*ones(1,3)]);
    M-M_true
end


Gplus = Y*U';

figure()
semilogy(abs(d));

S_RPY  = generate_RPY_matrix(rvec_out,0.1);
[V_RPY,D_RPY] = eig(S_RPY); 
hold on
semilogy(diag(D_RPY));

BB = (A+A')/2;
Binv = BB\eye(size(BB));
% %C = chol(B);
R_b = (K'*Binv*K);
M_b = R_b\eye(6*P);
M_corr_c = M_b*K'*Binv*BB*Binv'*K*M_b';
Ainv = A\eye(size(A)); 
R_b = (K'*Ainv*K);
M_b = R_b\eye(6*P);
M_corr_d = M_b*K'*Ainv*BB*Ainv'*K*M_b';
[VS,DS] =  eig(BB);
ds = diag(DS); 


tolvec = logspace(-15,-1,20);
for i = 1:length(tolvec)
    i
    tol = tolvec(i);
    [YT,UT]  = getPseudoFactors(A,tol,0);
    TGplus = YT*UT';
    AA = TGplus-TGplus';
    aa(i) = norm(AA,inf); %this is the symmetry error in the pseudoinverse
    %eTG_plus = eig(TGplus); 


    %Approach 1a based in eigendecomp
   %  ind = find(real(d)>tol);
   %  diff_set = setdiff(1:length(d),ind);
   %  d2 = 1./d; 
   %  d2(diff_set) = 0;
   % % A2 = V*diag(d2)*inv(V);
   %  TG_inv = V*diag(d2)*(V\eye(size(V)));
   %  R_b = K'*TG_inv*K;
   %  M_b = R_b\eye(6*P);
   % 
   %  M_corr1a = M_b*K'*TG_inv*A*TG_inv'*K*M_b';

    % Do eigendecomp of symmetriced matrix instead
    ind = find(real(ds)>tol);
    diff_set = setdiff(1:length(ds),ind);
    invd2 = 1./ds; 
    invd2(diff_set) = 0;

    d2 = ds;
    d2(diff_set) = 0;
   % A2 = V*diag(d2)*inv(V);
    TG_inv = VS*diag(invd2)*VS';
    TG_reg = VS*diag(d2)*VS';
    R_b = (K'*TGplus*K);
    M_b = R_b\eye(6*P);


    M_corr1a = M_b*K'*TGplus*TG_reg*TGplus'*K*M_b';
    erra(i) = norm(M_true-M_corr1a,2)/norm(M_true,2);

    M_corr1b = M_b*K'*TGplus*BB*TGplus'*K*M_b';
    errb(i) = norm(M_true-M_corr1b,2)/norm(M_true,2);

    %Approach 1b based on Eig decomp 

    %M_corr1b = M_b*K'*TGplus*A*TGplus'*K*M_b';
   % errb(i) = norm(M_true-M_corr1b,2)/norm(M_true,2);

    R_b = K'*TG_inv*K;
    M_b = R_b\eye(6*P);
    M_corr1c = M_b*K'*TG_inv*TG_reg*TG_inv'*K*M_b';
    errc(i) = norm(M_true-M_corr1c,2)/norm(M_true,2);

   
    %Approach 2: 
    R_2 = (K'*Gplus*B);
    M_2 = R_2\eye(6*P);
    % M_corr2 = M_2*K'*TGplus*A*TGplus'*K*M_2';
    % err2(i) = norm(M_true-M_corr2,2)/norm(M_true,2);

    %Approach 3:
    % M_corr3 = M_2*K'*Gplus*S*TGplus*S'*Gplus'*K*M_2';
    % err3(i) = norm(M_true-M_corr3,2)/norm(M_true,2);

    %Approach 4 (standard solve of problem)
    M_corr4 = M_2*K'*Gplus*S*TG_inv*S'*Gplus'*K*M_2';
    err4(i) = norm(M_true-M_corr4,2)/norm(M_true,2);

    M_corr5 = M_2*K'*Gplus*S*Binv*S'*Gplus'*K*M_2';
    err5(i) = norm(M_true-M_corr4,2)/norm(M_true,2);

    

end
%%
figure()
loglog(tolvec,erra,'+-','DisplayName','Pseduoinv (SVD) TS solve -- eigendecomp sym(TS) noise');
hold on
loglog(tolvec,errb,'*-','DisplayName','Pseduoinv (SVD) TS solve -- TS noise (no regularisation)');
loglog(tolvec,errc,'.-','DisplayName','Pseduoinv (Eigen) sym(TS) solve -- eigendecomp sym(TS) noise');
%loglog(tolvec,err2,'*-');
%loglog(tolvec,err3,'o-');
loglog(tolvec,err4,'o-','DisplayName','Pseudoinv (SVD) S solve -- eigendecomp TS noise');
loglog(tolvec,err5,'o-','DisplayName','Pseudoinv (SVD) S solve -- eigendecomp sym(TS) noise');
axis tight
xlabel('TOL','Interpreter','latex')
ylabel('$\|M-UU^T\|_2/\|M\|_2$','interpreter','latex')
ylim([1e-11,1])
grid on
legend
set(legend,'interpreter','latex')
%Do the same with an eigenvalue factorization! 

figure()
loglog(tolvec,aa)
xlabel('TOL','Interpreter','latex')
ylabel('Symmetry error in SVD pseudoinverse of TS','Interpreter','latex')
axis tight
grid on 
