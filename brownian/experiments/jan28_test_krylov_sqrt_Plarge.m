%test sqrt of traction computation

clear; close all;
rng(5); 

P = 10; %number of bodies
delta = 1; %smallest particle particle distance 
if P == 2
    q = [0 0 0; 2+delta 0 0]; %center coordiante matrix for P particles, x,y,z: size P x 3
elseif P == 1
    q = [0 0 0]; 
else
    %random configurations
    [q,~] = grow_cluster(P,delta); %Every particle has at least one neigbour at distance delta
end


  
fmm = 0; %only activate if many particles (say, more than 40)

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
Ti = stokes_DLP(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:));
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
    C = kron(eye(P),Ci);
    Cplus = kron(eye(P),Ci_plus);

    if dense

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

lgd = legend('show','Location','best');
lgd.Title.String = 'Eigenvalue truncation (preconditioner)';

axis tight
grid on
xlabel('Iteration number')
ylabel('Estimated error using Lanczos')

xsrc = r(:,1); 
ysrc = r(:,2);
zsrc = r(:,3);
xtar = rout(:,1); 
ytar = rout(:,2); 
ztar = rout(:,3); 
nx = nn(:,1);
ny = nn(:,2); 
nz = nn(:,3); 
mu = 1; 

Nsrc = numel(xsrc);
Ntar = numel(xtar);
T = zeros(3*Nsrc,3*Ntar);

prefac = 6 / (8*pi*mu);

for j = 1:Ntar
    %Vector differences r = x - y_j
    rx = xtar(j) - xsrc;
    ry = ytar(j) - ysrc;
    rz = ztar(j) - zsrc;
    r2 = rx.^2 + ry.^2 + rz.^2;
    r5 = sqrt(r2).^5;

    %Dot with normal at target
    rdotn = rx*nx(j) + ry*ny(j) + rz*nz(j);

    %Compute 3x3 block per target
    Txx = prefac * rdotn .* (rx.*rx) ./ r5;
    Txy = prefac * rdotn .* (rx.*ry) ./ r5;
    Txz = prefac * rdotn .* (rx.*rz) ./ r5;
    Tyx = prefac * rdotn .* (ry.*rx) ./ r5;
    Tyy = prefac * rdotn .* (ry.*ry) ./ r5;
    Tyz = prefac * rdotn .* (ry.*rz) ./ r5;
    Tzx = prefac * rdotn .* (rz.*rx) ./ r5;
    Tzy = prefac * rdotn .* (rz.*ry) ./ r5;
    Tzz = prefac * rdotn .* (rz.*rz) ./ r5;

    %Assemble into block structure
    T(:, 3*j-2:3*j) = reshape( ...
    [Txx; Txy; Txz;  ...
     Tyx; Tyy; Tyz;  ...
     Tzx; Tzy; Tzz], ...
    3, 3*Nsrc).';
end


end