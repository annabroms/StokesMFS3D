clear; 
close all;

%Test of the symmetry error in the matrix TG (visualising its eigenvalues). We also check the
%eigenvalues and singular values of the symmetric part of TG and the
%corresponding singular vectors. Conclusion: TG, G and sym(TG) all have a
%joint null-space, namely in the normal direction of the proxy surface. 

P=1;
delta = 0.2; 
if P==1
    q = [0 0 0];
else
    q = [0 0 0; 2+delta 0 0];
    q = q+[1 1 1];

    [q,B] = grow_cluster(P,delta);
    q = q+rand(1,3); 
end



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


for a = [1 1.2 1.5 2 3] %Determines oversampling factor for the collocation points

    [rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points
    N = size(rvec_in,1)/P;
    MM = size(rvec_out,1)/P; 

    if P == 1
        n = rvec_out;
    else
        n = rvec_out - kron(q,ones(MM,1));
    end

    % Visualise geometry
    figure(10)
    scatter3(rvec_in(:,1),rvec_in(:,2),rvec_in(:,3));
    hold on
    scatter3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3));
    axis equal
    quiver3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3),n(:,1),n(:,2),n(:,3));
    
    % Attempt at removing the null space of Tt (per-body correction)
    if P==1
        Tt = getTractionMat(rvec_in,rvec_out,n); 
        Nmat = normal_outer_block(rvec_in(1:N,:), rvec_out(1:MM,:))';
        %Tt = Tt+Nmat;
    else
        Tt = getTractionMat(rvec_in,rvec_out,n); 
        Nmat_block = normal_outer_block(rvec_in(1:N,:), rvec_out(1:MM,:))';
        %Tt = Tt + kron(eye(P), Nmat_block);
    end
    %Tt = getTractionMat(rvec_in,rvec_out,n);

    G = stokes_SLP_mat(rvec_in, rvec_out);

    A = 4*pi/MM*Tt'*G; 

    %% Visualise all singvals
    [Y,U]  = getPseudoFactors(G,tol,4);
    [YA,UA]  = getPseudoFactors(A,tol,3);

    %% Visualise all eigvals
    [V,D] = eig(A); 
    figure(1)
    plot(real(diag(D)),imag(diag(D)),'+-');
    hold on

    [Vsym,Dsym] = eig((A+A')/2); 
    dsym = diag(Dsym);
    figure(2)
    semilogy(dsym,'+-');
    hold on

    %% Visualize smallest eigen/singular vectors
    n_modes = 1;
    if P == 2
        n_modes = 2;
    end

    [~,orda] = sort(abs(diag(D)),'ascend');
    va1 = V(:,orda(1));
    va2 = [];
    if n_modes == 2
        va2 = V(:,orda(2));
    end
    plot_quiver_field_pair(rvec_in, va1, va2, 'TG eigvecs (blue=min eig, red=2nd)');

    [~,ordas] = sort(abs(diag(Dsym)),'ascend');
    vas1 = Vsym(:,ordas(1));
    vas2 = [];
    if n_modes == 2
        vas2 = Vsym(:,ordas(2));
    end
    plot_quiver_field_pair(rvec_in, vas1, vas2, 'sym(TG) eigvecs (blue=min eig, red=2nd)');

    [Ug,Sg,Vg] = svd(G,'econ');
    vmin_g1 = Vg(:,end);
    vmin_g2 = [];
    if n_modes == 2
        vmin_g2 = Vg(:,end-1);
    end
    plot_quiver_field_pair(rvec_in, vmin_g1, vmin_g2, 'G right singvecs (blue=min sing, red=2nd)');

    [Uc,Sc,Vc] = svd(A,'econ');
    vmin_a1 = Vc(:,end);
    vmin_a2 = [];
    if n_modes == 2
        vmin_a2 = Vc(:,end-1);
    end
    plot_quiver_field_pair(rvec_in, vmin_a1, vmin_a2, 'TG right singvecs (blue=min sing, red=2nd)');
    

end

%%
figure(2)
xlabel('j','Interpreter','latex')
ylabel('$\lambda_j$','Interpreter','latex')
title('Eigvals of sym(TG)')

figure(1)
xlabel('Re($\lambda$)','Interpreter','latex')
ylabel('Im($\lambda$)','Interpreter','latex')
title('Eigvals of TG')

visualise = 1; 
tol = 1e-13;

legend('a = 1', 'a = 1.2','a = 1.5', 'a = 2', 'a = 3')
set(legend,'interpreter','latex')

function plot_quiver_field_pair(rvec, vec1, vec2, plot_title)
    if isempty(rvec) || isempty(vec1)
        return;
    end
    v1 = real(vec1(:));
    n = size(rvec,1);
    if numel(v1) ~= 3*n
        return;
    end
    v1 = reshape(v1,3,[])';
    figure();
    quiver3(rvec(:,1), rvec(:,2), rvec(:,3), v1(:,1), v1(:,2), v1(:,3), 'AutoScale','on','Color','b');
    hold on
    if ~isempty(vec2)
        v2 = real(vec2(:));
        if numel(v2) == 3*n
            v2 = reshape(v2,3,[])';
            quiver3(rvec(:,1), rvec(:,2), rvec(:,3), v2(:,1), v2(:,2), v2(:,3), 'AutoScale','on','Color','r');
        end
    end
    axis equal
    grid on
    title(plot_title)
end

fprintf(['Fig 1 and 2 display the eigvals of TG and sym(TG),\n Fig 4 and 3 the singvals of the same matrices\n ...' ...
    'There is clearly a nontrivial null space... now to be visualised'])


alignfigs; 
