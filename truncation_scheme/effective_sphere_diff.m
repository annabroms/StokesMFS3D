function f = effective_sphere_diff(Rg,vars,mob,relative)
%EFFECTIVE_SPHERE_DIFF_MOB returns the difference between the qbx and multiblob
%mobility coefficients for a sphere


N = vars.N;
[Rt,Rr] = sphere_radius(Rg,N,vars.scale);


%f = [1-Rt; 1-Rf];
if ~mob
    f = [6*pi-6*pi*Rt; 8*pi-8*pi*Rr^3];
else
    f = [1/(6*pi*Rt)-1/6/pi; 1/(8*pi*Rr^3)-1/8/pi];
   % f = [1/(6*pi*Rt)-1/6/pi; 1/(8*pi*Rr^3)-1/8/pi; Rt-Rr];
    %f = [1/(6*pi*Rt)-1/6/pi]; %if only minimising the translational error
end

if relative
    if ~mob
        f = [6*pi-6*pi*Rt; 8*pi-8*pi*Rr^3]./[6*pi; 8*pi];
    else
        f = [1/(6*pi*Rt)-1/6/pi; 1/(8*pi*Rr^3)-1/8/pi; Rt-Rr]./[1/6/pi; 1/8/pi; 1];
        f = [1/(6*pi*Rt)-1/6/pi; 1/(8*pi*Rr^3)-1/8/pi]./[1/6/pi; 1/8/pi];
    end
end

end
    
