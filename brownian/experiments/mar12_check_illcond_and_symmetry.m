clear; 
close all;

%Test of the symmetry error in the matrix TG (visualising its eigenvalues). We also check the
%eigenvalues and singular values of the symmetric part of TG and the
%corresponding singular vectors. This script runs the test on a small cluster of ellipsoids.

sphere = 0;
tdesign = 0; %uniform discretisation of sphere

P = 1;
delta = 1;
E0 = [.4 .6 1]; %Type T in the mobility paper
sep = 0.125;
%sep = 0.2;

if sphere
    E0 = [1 1 1]; % this is a sphere
    sep = 0.9;
    sep = 0.3; 
    %sep = 0.2; %test for a sphere
end

% Set resolution (following ellipsoid_mobility_run)
Nv = 40;
Nv = 15; % need a much coarser resolution



rng(8);
[~, R, qvec, ~] = ellipsoid_cluster(E0,P,delta);

disp('Geometry created')

tol = 1e-15;
visualise = 1;
show_dir = 1; %visualise eigenvector/singvector corresponding to smallest eigval / singval?


for a = [1.5 2] %Determines oversampling factor for the collocation points

    % Get discretization and normals
    N1 = a*0.75*Nv;
    N2 = a*Nv;
    if sphere && tdesign
        Nt = 336;
        [rvec_in,rvec_out,opt] = init_spheres(vertcat(qvec{:})',1-sep,Nt,a);
        M = length(rvec_out)/P;
        n = rvec_out - kron(vertcat(qvec{:})',ones(M,1));
           
    else
        [rvec_in,rvec_out,~,~,~,w,n] = getEllipsoidGrids(E0,P,delta,N1,N2,sep,R,qvec); %Assign source and collocation points
        M = length(rvec_out)/P;
    end
    

    % Visualise geometry
    figure(10)
    scatter3(rvec_in(:,1),rvec_in(:,2),rvec_in(:,3));
    hold on
    scatter3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3));
    axis equal
    quiver3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3),n(:,1),n(:,2),n(:,3));

    Tt = getTractionMat(rvec_in,rvec_out,n); 
    G = stokes_SLP_mat(rvec_in, rvec_out);
    [C,TW,A] = getSymTWS(rvec_out,rvec_in,n,w);

    %% Visualise all singvals
    [Y,U]  = getPseudoFactors(G,tol,5);
    [YA,UA]  = getPseudoFactors(A,tol,4);
    [YB,UB]  = getPseudoFactors(C,tol,3);

    %% Visualise all eigvals
    [V,D] = eig(A); 
    figure(1)
    plot(real(diag(D)),imag(diag(D)),'+-');
    hold on
   
    [Vsym,Dsym] = eig(C); 
    dsym = diag(Dsym);
    figure(2)
    subplot(1,2,1)
    semilogy(dsym,'+-');
    hold on
    subplot(1,2,2)
    semilogy(abs(dsym),'+-');
    hold on


    %% Visualize smallest eigen/singular vectors
    if show_dir
        n_modes = 1;
    
        % [~,orda] = sort(abs(diag(D)),'ascend');
        % va1 = V(:,orda(1));
        % va2 = [];
        % if n_modes == 2
        %     va2 = V(:,orda(2));
        % end
        % plot_quiver_field_pair(rvec_in, va1, va2, 'TG eigvecs (blue=min eig, red=2nd)');
        % 
        % [~,ordas] = sort(abs(diag(Dsym)),'ascend');
        % vas1 = Vsym(:,ordas(1));
        % vas2 = [];
        % if n_modes == 2
        %     vas2 = Vsym(:,ordas(2));
        % end
        % plot_quiver_field_pair(rvec_in, vas1, vas2, 'sym(TG) eigvecs (blue=min eig, red=2nd)');
        % 
        % [Ug,Sg,Vg] = svd(G,'econ');
        % vmin_g1 = Vg(:,end);
        % vmin_g2 = [];
        % if n_modes == 2
        %     vmin_g2 = Vg(:,end-1);
        % end
        % plot_quiver_field_pair(rvec_in, vmin_g1, vmin_g2, 'G right singvecs (blue=min sing, red=2nd)');
    
        [Uc,Sc,Vc] = svd(A,'econ');
        vmin_a1 = Vc(:,end);
        vmin_a2 = [];
        if n_modes == 2
            vmin_a2 = Vc(:,end-1);
        end
        plot_quiver_field_pair(rvec_in, vmin_a1, vmin_a2, 'TG right singvecs (blue=min sing, red=2nd)');
        hold on
        scatter3(rvec_in(:,1),rvec_in(:,2),rvec_in(:,3))
    end

end

%%
figure(2)
xlabel('j','Interpreter','latex')
ylabel('$\lambda_j$','Interpreter','latex')
sgtitle('Eigvals of sym(TG)')

figure(1)
xlabel('Re($\lambda$)','Interpreter','latex')
ylabel('Im($\lambda$)','Interpreter','latex')
title('Eigvals of TG')

legend('a = 2')
set(legend,'interpreter','latex')

fprintf(['Fig 1 and 2 display the eigvals of TG and sym(TG),\n Fig 4 and 3 the singvals of the same matrices, singvals of G in Fig 5\n ...' ...
    'There is clearly a nontrivial null space... now to be visualised'])

alignfigs; 

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

