function [c,s]=rotmat(a,b)
%Help function in helsing_gmres_mv. Computes givens rotations
if  b==0
c=1;
s=0;
elseif abs(b)>abs(a)
temp=a/b;
s=1/sqrt(1+temp^2);
c=temp*s;
else
temp=b/a;
c=1/sqrt(1+temp^2);
s=temp*c;
end
end
