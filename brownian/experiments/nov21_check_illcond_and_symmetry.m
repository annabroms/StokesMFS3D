clear; 
close all;

%Test of the symmetry error in the matrix A = TS. We also check the
%eigenvalues and singular values

P=2;
delta = 1; 
if P==1
    q = [0 0 0];
else
    delta = 3; %initial particle particle distance
    q = [0 0 0; 2+delta 0 0];
end

[q,B] = grow_cluster(P,delta);

% Set resolution
Rp = 0.63; %fine resolution
N = 700;

% Rp = 0.15; %Proxy sphere radius -- very coarse resolution
% N = 50;  % Number of proxy sources
% % 
Rp = 0.30;
%Rp = 0.10; 
N = 100; 

tol = 1e-15;
visualise = 1; 


for a = 2 %; %[1 1.2 1.5 2 3] %Determines oversampling factor for the collocation points

    [rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points
    N = size(rvec_in,1)/P;
    MM = size(rvec_out,1)/P; 
    
    n = repmat(rvec_out(1:MM,:),P,1);
    %remember T is a block diagonal matrix
    if P==1
        T = getTraction(rvec_in,rvec_out,n); 
        Nmat = normal_outer_block(rvec_in(1:N,:), rvec_out(1:MM,:));
        T = T+Nmat;
    else
        Tblock = getTraction(rvec_in(1:N,:),rvec_out(1:MM,:),rvec_out(1:MM,:));
        Nmat = normal_outer_block(rvec_in(1:N,:), rvec_out(1:MM,:));
        %Tblock = Tblock+Nmat; 
        T = kron(eye(P),Tblock);
    end
    %T = getTraction(rvec_in,rvec_out,n);

    S = generate_stokes_mat(rvec_in, rvec_out);

    A = -4*pi/MM*T*S;

    [Y,U]  = getPseudoFactors(S,tol,4);
    [YA,UA]  = getPseudoFactors(A,tol,3);

    [V,D] = eig(A); 
    figure(1)
    plot(real(diag(D)),imag(diag(D)),'+-');
    hold on

    [V,D] = eig((A+A')/2); 
    figure(2)
    semilogy(diag(D),'+-');
    hold on

end

%%
xlabel('j','Interpreter','latex')
ylabel('$\lambda_j$','Interpreter','latex')

figure(1)
xlabel('Re($\lambda$)','Interpreter','latex')
ylabel('Im($\lambda$)','Interpreter','latex')

visualise = 1; 
tol = 1e-13;

legend('a = 1', 'a = 1.2','a = 1.5', 'a = 2', 'a = 3')
set(legend,'interpreter','latex')

