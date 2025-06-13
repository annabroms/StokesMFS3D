function ellipsoid_mobility_run(P,delta,Nv,visualise,solve_res,read_name,save_name)
%ELLIPSOID_MOBILITY_RUN Runs a mobility problem for a cluster of ellipsoids.
%
%   ELLIPSOID_MOBILITY_RUN(P, delta, Nv, visualise, solve_res,read_name,save_name)
%   solves a Stokes mobility problem for a cluster of P ellipsoidal
%   particles using MFS. The ellipsoids are discretized with a
%   quasi-uniform grid based on Gauss–Legendre and trapezoidal nodes.
%
%   INPUTS:
%       P               - Number of ellipsoidal particles.
%       delta           - Smallest particle-particle distance
%       Nv              - Number of Gauss–Legendre nodes used in one of the
%                         directions around the particle.
%       visualise       - plot geometry
%       solve_res       - in addition, solve a resistance problem. 
%       read_name       - string: If nonempty, the particle configuration
%                         and force/torque vector are loaded from file with
%                         this name
%                         (rather than generated within the function).
%       save_name     - string: save to file with this name
%
%   NOTES:
%       - In the large-scale test shown in Example 4 / Fig. 1 of the
%         mobility MFS paper, we use P = 10000, Nv = 34 and solve_res = 0.
%         (A resistance problem with many particles requires many GMRES
%         iterations).
%       - Both the proxy and collocation surfaces are discretized using the
%         quasi-uniform ellipsoid grid described by Stein and Barnett
%         in [Stein2022: Quadrature by Fundamental solutions, ACOM]. Each ellipsoid 
%         is parameterized in Cartesian coordinates as
%             (a*sqrt(1 - t^2)*cos(s), b*sqrt(1 - t^2)*sin(s), c*t),
%         where (s, t) ∈ [0, 2π] × [-1, 1]. The parameter t is discretized
%         using Gauss–Legendre nodes, while s uses a
%         periodic trapezoidal rule.
%
%   To test, call without arguments
%
% Anna Broms, June 13, 2025

%% Set aspect ratio of ellipsoid and distance to others
E0 = [.5 .5 1]; %Type S in the mobility paper
%E0 = [1 1 1]; % Sphere
E0 = [.4 .6 1]; %Type T in the mobility paper


if nargin<1
   disp("Runs a basic demo")
   show_basic_test(); 
   return;
elseif nargin<2
    error("You must specify number of particles and discretisation");
elseif nargin <3
    delta = 1; 
    solve_res = 0;
elseif nargin<4
    visualise = 0;
    save_name = "test_ellipsoids";
    solve_res = 0; 
elseif nargin <5
    solve_res = 0;
elseif nargin<6
    save_name = "test_ellipsoids";
elseif nargin < 7
    save_name ="test_ellipsoids";
end

if ~isempty(read_name)
    read_from_file = 1; 
    read_file = sprintf('%s_P%u_dmin%1.2d_E0_%1.2d_%1.2d_%1.2d_N%u.mat',read_name,P,delta,E0(1),E0(2),E0(3),Nv);
else
    read_from_file = 0;
end
save_file = sprintf('%s_P%u_dmin%1.2d_E0_%1.2d_%1.2d_%1.2d_N%u.mat',save_name,P,delta,E0(1),E0(2),E0(3),Nv);


rng(8);

profile = 0; 
if profile
    memorygraph('start')
end
%

%For many particles, enable fmm
if P>40
    fmm = 1; 
else 
    fmm = 0;
end


if read_from_file
    f = load(read_file);
    qvec = f.qvec; %center coordinates
    R = f.R; %rotation matrices
    Fvec = f.Fvec; %vector of forces and torques
    for k = 1:P
        q(k,:) = qvec{k}';
    end
else
    Pmax = 500;
    if P<Pmax
        [E R qvec xnear] = ellipsoid_cluster(E0,P,delta);
        for k = 1:P
            q(k,:) = qvec{k}';
        end
    else
        %If the cluster is really large, it takes a very long time to generate
        %the geometry using `ellipsoid_cluster.m´. Then, instead, we first
        %generate a suspension of spherical shells with minimum separation
        %delta to at least one neighbour. Then, an ellipsoid is generated at
        %random within each such sphere. 
        [q,B] = grow_cluster(P,delta);
        
        Rz = @(t) [cos(t) -sin(t) 0;sin(t) cos(t) 0; 0 0 1];  % z-ax rot mat
        Ry = @(t) [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)];  % y-ax rot mat
        for k = 1:P
            u = randn(3,1); u=u/norm(u);  % rand on S2 defines (alpha, beta) Euler
            alpha = atan2(u(2),u(1));
            beta = acos(u(3));
            gamma = 2*pi*rand;
            R{k} = Rz(alpha) * Ry(beta) * Rz(gamma); 
            %R{k} = eye(3); 
            qvec{k} = q(k,:)';
            
        end
    end
    Fvec = rand(6*P,1);
end

%Now we have the center coordinates and orientations of everybody and can
%generate the proxy and collocation points. 
sep = 0.125; %separation between proxy and collocation surfaces (in normal direction) good choice for E0 = [.5 .5 1];
%sep = 0.13;
[rin,rout,~,~] = getEllipsoidGrids(E0,P,delta,0.75*Nv,Nv,sep,R,qvec);

disp('Grids computed...')

%% Solve mobility problem

opt = init_MFS();
opt.ellipsoid = 1; 

if solve_res
    Uref = rand(6*P,1);
    [Fvec, iter_r, lambda_norm_r, uerr_r] = solve_resistance(q,rin,rout,Uref, opt,R,E0);
    [rin,rout,~,~] = getEllipsoidGrids(E0,P,delta,0.75*Nv,Nv,sep*1.05,R,qvec); %avoid "inverse crimes"
    [Uvec, iter_m, lambda_norm, uerr] = solve_mobility(q,rin,rout,Fvec, opt,R,E0);
    fprintf("Two way error: %1.2e, resistance residual %1.2e, mobility residual %1.2e\n",...
        norm(Uref-Uvec,inf)/norm(Uref,inf),uerr_r,uerr)
else
    %only solve the mobility problem
    [Uvec, iter_m, lambda_norm, uerr] = solve_mobility(q,rin,rout,Fvec, opt,R,E0);
end




%% Visualise? 

if visualise 
    figure()
    %res = 30;
    res = 8;
    cross_mat = @(x) -[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
    for i = 1:P
        [x, y, z] = ellipsoid_matlab(0, 0, 0, E0(1), E0(2), E0(3), res);
        XX = qvec{i}'+(R{i}*[x(:) y(:) z(:)]')';


        Ui = Uvec((i-1)*6+1:i*6);
        
        for k = 1:length(XX)        
            Usurf(k) = norm(Ui(1:3)+cross_mat(XX(k,:)-qvec{i}')*Ui(4:6));
        end

        X = reshape(XX(:,1),res+1,res+1);
        Y = reshape(XX(:,2),res+1,res+1);
        Z = reshape(XX(:,3),res+1,res+1);
        Us = reshape(Usurf,res+1,res+1);

        surf(X,Y,Z,Us);
        hold on

       
    end
    shading interp
    axis equal
    camlight   
    hold on
    
    c = colorbar;
    ylabel(c,'Surface velocity magnitude from computed RBM','interpreter','latex','FontSize',16)
    axis tight
    axis off

end

if profile
    memorygraph('label','finished mobility')
    [bytes est_times cpu_times cpu_usages labelstrings labeltimes] = memorygraph('plot');
    maxRAM = max(bytes);
else
    labeltimes = [];
    labelstrings = {};
    bytes = [];
    maxRAM = [];
end

%end
if ~isempty(save_name)
    save(save_file,'iter_m','E0','Uvec','Fvec','labeltimes','labelstrings','sep','Nv',...
        'R','qvec','delta','maxRAM')
end

if profile
    memorygraph('done');
end


end

function show_basic_test()
    P = 7;
    Nv = 40;
    delta = 0.5; 
    visualise = 1;
    solve_res = 1; %in addition to mobility, solve a resistance problem and display error metrics
    read_name = [];
    %read_name = "test_ellipsoids";
    save_name = "test_ellipsoids";
    save_name = '';
    ellipsoid_mobility_run(P,delta,Nv,visualise,solve_res,read_name,save_name)
end



