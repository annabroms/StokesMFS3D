function res = getFlow(tau_stokes,rin,rout,vars)
%getFlow 
%Compute velocity field given source
%strengths tau_stokes located in source points rin evaluated in target
%points rout. If vars.fmm FMM3D is used for the evaluation. 
if vars.fmm 
    nd = 1;
    srcinfo.nd = nd;
    
    ifppreg = 0;
    ifppregtarg = 1;
    
    srcinfo.sources = rin';
    srcinfo.stoklet = reshape(tau_stokes,3,[]);
    
    targ = rout';  
    eps = vars.eps; % was -6
    
    
    U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);    
    %U = st3ddir(srcinfo,targ,ifppregtarg); %Try to use this one

    res = U.pottarg(:);

    clear U srcinfo;
else
    targ = rout; 
    srcinfo.stoklet = reshape(tau_stokes,3,[]);
    U = SE0P_Stokeslet_direct_full_ext_mex(rin, srcinfo.stoklet', struct('eval_ext_x', targ));
       % instead for comparison
    U = U';
    res = 1/8/pi*U(:);

end




end