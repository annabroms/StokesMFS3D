%test sqrt of traction computation

clear; close all;

P = 1; %number of bodies
delta = 1; %smallest particle particle distance 
q = [0 0 0; 2+delta 0 0]; %center coordiante matrix for P particles, x,y,z: size P x 3
q = [0 0 0]; 

%random configurations
%[q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
  
fmm = 0; %only activate if many particles (say, more than 40)

%Test first with very low resolution
Rp = 0.30;
N = 100; 
a = 1.2; 
a = 2; %or play with SVD truncation level

Rp = 0.68; %proxy radius
N = 700; % approximate number of proxy sources on every particle


[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a);
M = size(rvec_out,1)/P; 
n = repmat(rvec_out(1:M,:),P,1);
T = getTraction(rvec_in,rvec_out,n); 
S = generate_stokes_mat(rvec_in, rvec_out);

A  = -T*S; 

%Do eigval decomp
[V,D] = eig(A);
d = diag(D);
tol = 1e-8; 
ind = find(real(d)>tol);
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
[y,iter_errv1] = KrylovSqrtMsing(S,dW1,tol);
[y,iter_errv2] = KrylovSqrtMsing(A,dW2,tol);
y2 = krylov_sqrt(A,dW2,20);
norm(y-y2,inf)/norm(y,inf)
%%
figure()
semilogy(iter_errv1,'+-');
hold on
semilogy(iter_errv2,'o-');
xlabel('Iteration number')
ylabel('Estimated error using Lanczos')
axis tight
legend('RPY tensor','TS')



