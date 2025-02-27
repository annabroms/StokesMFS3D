function [rbase_in,rbase_out,body_inner,weights] = getGrids(Rg,opt)
%Returns an inner and outer grid for the sphere, weights is a diagonal
%vector of quadrature weights. Uses 1. a grid of fibonacci nodes, 2.
%spherical designs or 3.  quasi-uniform discretisation of an ellipsoid from
%the work of Stein and Barnett 2022. 2. is preferred as it gives a more
%uniform distribution.

if opt.fib
    if opt.fib2
        [x,y,z] = fibonacci_spiral2(ceil(opt.ls_glob*opt.fib_n),2*pi,1);
    else
        [x,y,z] = fibonacci_spiral(ceil(opt.ls_glob*opt.fib_n),2*pi,1);
    end
    [~,~,~,weights,~,~,~,~,~,~,~,~,~]=voronoiSphereArea2(x',y',z'); 

    rbase_out = [x' y' z']; 
    
    %just for debugging minimum distances
    % for k = 1:size(rbase_out,1)
    %     d(:,k) = vecnorm(rbase_out(k,:)-rbase_out,2,2);
    % end
    % d = d(:); 
    % d = d(d>0); 
    % sqrt(size(rbase_out,1))*min(d);



    %rotate grid with a given random orientation
    %v = [0.7315    0.9748    0.1534];
    v = [1 0 0];
    q = vecToQuat(v);
    R = quatToMat(q);
    rbase_out = (R*rbase_out')';
    if opt.fib2
        [x,y,z] = fibonacci_spiral2(opt.fib_n,2*pi,Rg);
    else
        [x,y,z] = fibonacci_spiral(opt.fib_n,2*pi,Rg);
    end
    
    rbase_in = [x' y' z']; 
    %v = [ 0.4019  0.7301 0.9551];
    v = [1 0 0];
    q = vecToQuat(v);
    R = quatToMat(q);
    rbase_in = (R*rbase_in')';
    %n_in  = [x' y' z']./Rg; 
    %n_in = (R*n_in')';

    body_inner = [];
elseif opt.design
    [X,w] = get_sphdesign(opt.des_n);
    rbase_in = Rg*X;


    [X,w] = get_sphdesign(ceil(opt.ls_glob*opt.des_n));
    %Check minimum distances
    % for k = 1:size(X,1)
    %     d(:,k) = vecnorm(X(k,:)-X,2,2);
    % end
    % d = d(:); 
    % d = d(d>0); 
    % sqrt(size(X,1))*min(d)
    rbase_out = X; 
    weights = w;
    body_inner = [];
else
    
    nph = opt.nph;
    nth = opt.nth;
    body_inner = ellipsoid(Rg,Rg,Rg);   % baseline object at the oridin, aligned
    body_inner = setupsurfquad(body_inner,[nph,nth]);
    
    body_outer = ellipsoid(1,1,1);   % baseline object at the oridin, aligned
    body_outer = setupsurfquad(body_outer,[round(opt.ls*nph),round(opt.ls*nth)]);
    
    rbase_in = body_inner.x';
    rbase_out = body_outer.x';
   
    [x,y,z] = fibonacci_spiral(opt.fib_n,2*pi,Rg);
    rbase_in = [x' y' z']; 
    v = [ 0.4019  0.7301 0.9551];
    q = vecToQuat(v);
    R = quatToMat(q);
    rbase_in = (R*rbase_in')';
    
    if opt.stokes > 1 %enables the use of a double layer formulation
       % n_in = body_inner.n;
        n_in = body_inner.nx'; %Barnett
    else
        n_in = [];
    end

end

if opt.plot
    figure()
    %KTHlred!80!KTHyellow
    red2 = 0.8*[231,51,57]+0.2*[251,186,0];
    scatter3(rbase_in(:,1),rbase_in(:,2),rbase_in(:,3),'filled','MarkerFaceColor','b')
    hold on
    scatter3(rbase_out(:,1),rbase_out(:,2),rbase_out(:,3),'MarkerEdgeColor','k');
   % [25,105,188]
end



% if opt.quad
%     WW = repmat(body_inner.w,1,3);
%     WW = WW';
%     WW = WW(:);
% else
%     WW = ones(1,3*size(rbase_in,1)); 
% end
% 
% Q = diag(WW);

end