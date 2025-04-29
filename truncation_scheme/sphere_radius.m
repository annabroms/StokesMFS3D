function [Rt,Rr] = sphere_radius(Rg,N,scale)

    if nargin<3
        scale = 1;
    end
    
    opt = init_MFS(N);
    [rin_c,~,~,~] =  getGrids(Rg,opt);
    

    R_blob = track_resistance_direct(rin_c,[0,0,0],scale,0);

    Rt = mean(diag(R_blob(1:3,1:3)))/(6*pi);
    Rr = (mean(diag(R_blob(4:6,4:6)))/(8*pi))^(1/3);
% 
%     Rh = max(Rt,Rf);
%     F = 1-Rh;

end